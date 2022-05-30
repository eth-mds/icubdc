
bin_labels <- function(breaks, unit, lower0 = TRUE) {

  x_labels <- sapply(1:(length(breaks)-1),
    function(x) paste0("[", breaks[x], "-", breaks[x+1], "]")
  )
  first_label <- paste0("< ", breaks[1])
  if (lower0) first_label <- paste0("[0-", breaks[1], "]")
  x_labels <- c(first_label, x_labels, paste0("> ", breaks[length(breaks)]))
  x_labels <- paste(x_labels, unit)

  return(x_labels)
}

reshape_frame <- function(frame, break_list) {

  for(col in names(break_list)) {

    breaks <- break_list[[col]]

    x <- frame[, col]
    x <- as.factor(.bincode(x, breaks = c(-Inf, sort(breaks), Inf)))
    if(length(breaks) == 1 & breaks[1] == 0.5) {
      levels(x) <- c("_no", "_yes")
    } else {
      levels(x) <- paste(" >", c(0, breaks))
    }
    frame[, col] <- x

  }

  return(frame)

}

extract_score <- function(model, p.threshold = 0.05) {

  score <- model$coefficients * as.integer(summary(logit)$coefficients[, 4] < p.threshold)
  score <- score[-1]
  score <- score / min(abs(score[score != 0]))
  score <- round(score)

  score

}

discretize <- function(x, breaks) {

  if (is.null(breaks)) return(x)
  return(.bincode(x, c(-Inf, breaks, Inf)))

}

repl_na <- function(x, val) replace(x, is.na(x), val)

carry_fwd <- function(x, time = 1L) {

  res <- x[length(x)]

  if(is.na(res)) {

    cand <- which(!is.na(x))
    if(length(cand) == 0) return(res)
    cand <- max(cand)
    if ((length(x) - cand) <= time) res <- x[cand]

  }

  res

}

srcwrap <- function(src) {

  if (length(src) > 1L) {
    return(vapply(src, srcwrap, character(1L)))
  }

  switch(src,
    mimic = "MIMIC-III",
    miiv = "MIMIC-IV",
    eicu = "eICU",
    hirid = "HiRID",
    aumc = "AUMC",
    mimic_demo = "MIMIC Demo",
    eicu_demo = "eICU Demo",
    stop("unknown data source")
  )
}

median_na <- function(x) {
  if(all(is.na(x))) return(x[1])
  return(median(x, na.rm = T))
}

make_plots <- function(res, save_plots, folder, plot_name, bottom = NULL,
                       width = 12, height = 5, dpi.res = 300, n.cols = length(res)) {

  wd <- getwd()

  ranges <- lapply(res, function(x) ggplot_build(x)$layout$panel_scales_y[[1]]$range$range)
  y_limits <- Reduce(function(x, y) c(min(x, y), max(x, y)), ranges)
  res <- lapply(res, function(x) x + ylim(y_limits))
  plot <- plot_grid(plotlist = res,
    ncol = n.cols, labels = c("A", "B", "C", "D")[1:length(res)])
  grid.arrange(arrangeGrob(plot, bottom = bottom))

  if(save_plots) {

    ggsave(file.path(wd, "paper", "figures", folder, paste0(plot_name, ".tif")),
      device = "tiff", width = width, height = height, dpi = dpi.res)

  }

}

carry_forward_to <- function(x, limits, cols = names(limits)) {

  assert_that(
    is_ts_tbl(x), has_no_gaps(x), is_intish(limits), all(limits >= 0),
    length(limits) == length(cols), all(cols %in% data_vars(x)),
    !"lag" %in% colnames(x)
  )

  names(limits) <- cols

  limits <- limits[limits > 0]

  if (length(limits) == 0L) {
    return(x)
  }

  on.exit(x[, lag := NULL])

  for (i in seq_along(limits)) {

    lmt <- limits[i]
    col <- sub("_raw$", "_imp", names(lmt))

    msg_ts("carrying `", names(lmt), "` forward by max ",
           format(as.difftime(lmt, units = time_unit(x))))

    x <- x[, lag := get(names(lmt))]
    x <- x[, c(col) := lag]

    for (j in seq_len(lmt)) {
      x <- x[, lag := data.table::shift(lag), by = c(id_vars(x))]
      x <- x[is.na(get(col)) & !is.na(lag), c(col) := lag]
    }

    gc()
  }

  x
}

impute_values <- function(x, values, cols = names(values), new_cols = cols) {

  assert_that(
    is.numeric(values), all(cols %in% data_vars(x)),
    all(cols %in% names(values)), length(cols) == length(new_cols)
  )

  values <- values[cols]

  msg_ts("imputing missing values")

  x <- x[, c(new_cols) := Map(repl_na, .SD, values), .SDcols = cols]

  x
}

augment <- function(x, fun, suffix,
                    cols = grep("_raw$", colnames(x), value = TRUE),
                    names = sub("_raw$", paste0("_", suffix), cols),
                    by = NULL, win = NULL, ..., info = NULL) {

  if (not_null(info)) {
    msg_ts("augmentation step ", info)
  }

  if (is.character(fun)) {

    win <- ceiling(win / as.double(interval(x), units = "hours"))
    fun <- switch(fun, min = roll::roll_min,  max = roll::roll_max,
                      mean = roll::roll_mean, var = roll::roll_var,
                      any  = roll::roll_any,
                  stop("unknown augmentation function ", fun, "()"))

    x <- x[, c(names) := lapply(.SD, fun, win, min_obs = 1L),
           .SDcols = cols, by = c(id_vars(x))]

  } else {

    assert_that(is.null(win))

    if (is.null(win) && is.null(by)) {

      x <- x[, c(names) := lapply(.SD, fun, ...), .SDcols = c(cols)]

    } else if (is.null(win)) {

      x <- x[, c(names) := lapply(.SD, fun, ...), .SDcols = c(cols),
             by = c(by)]

    }
  }

  x
}

w_value <- function(source, concept, dir, verbose = FALSE,
  upto = hours(24L), patient_ids = NULL, imp_val = NULL) {


  tbl <- load_concepts(concept, source, patient_ids = patient_ids, verbose = verbose)

  if(is.null(imp_val)) imp_val <- median(tbl[[concept]], na.rm = T)
  if(grepl("epineph", concept)) imp_val <- 0

  if(is_ts_tbl(tbl)) tbl <- tbl[get(index_var(tbl)) <= upto]

  if(is.null(patient_ids)) patient_ids <- unique(stay_windows(source)[[id_vars(tbl)]])

  pts <- data.table::data.table(patient_ids)
  data.table::setnames(pts, names(pts), id_vars(tbl))

  pts <- as_id_tbl(pts)

  w_fun <- ifelse(dir == "increasing", max_or_na, min_or_na)

  res <- tbl[, w_fun(get(concept)), by = eval(id_vars(tbl))]
  data.table::setnames(res, "V1", "w_val")

  res <- merge(pts, res, by = id_vars(res), all.x = T)

  if(length(imp_val) == 2) {
    res[is.na(w_val), "w_val"] <- runif(sum(is.na(res[["w_val"]])), imp_val)
  } else {
    res[is.na(w_val), "w_val"] <- imp_val
  }


  res

}

hypo <- function(source, patient_ids, hypo.threshold = 3.9, upto = hours(72L), verbose = FALSE) {

  hypo <- load_concepts("glu", source, patient_ids = patient_ids, verbose = verbose)

  hypo[, hg := (glu <= 18.016*hypo.threshold)]
  hypo <- hypo[hg == T & get(index_var(hypo)) <= upto, head(.SD, n = 1L), by = eval(id_vars(hypo))]

  hypo[, c(meta_vars(hypo), "hg"), with = FALSE]

}

fwrap <- function(feat) {
  if(feat == "alt") return("ALT")
  if(feat == "ast") return("AST")
  else return(feat)
}

diff_means <- function(d, ind, threshold) {
  use <- d[ind, ]
  mean(!is.na(use[w_val > threshold, "hg"])) - mean(!is.na(use[w_val <= threshold, "hg"]))
}

ate_sa <- function(d, ind, threshold) {
  use <- d[ind, ]

  H <- as.integer(!is.na(use[["hg"]]))
  L <- as.integer(use[["w_val"]] > threshold)
  S <- use[["w_sofa"]]
  S[is.na(S)] <- 0
  S[S > 12] <- 13

  assert_that(length(unique(S)) == 14L)
  tabular <- table(L, S)
  propensity <- tabular[1, ] / (tabular[2, ] + tabular[1, ])
  mean(H*(2*L-1) / (propensity[S+1]*(1-L) + (1-propensity[S+1])*L))

}

lwrap <- function(lab) {
  if(lab == "liver") return("Liver dysfunction")
  if(lab == "shock_state") return("Shock")
  if(lab == "ins_therapy") return("Insulin therapy")

  lab

}

compute_plot <- function(tbl, condition, legend.source, x.pos, y.pos,
  bins = NULL, unit = NULL, legend.rel.size = 1) {

    print(paste("Dataset:", srcwrap(unique(tbl[["dataset"]]))))
    print(paste("Conditioning on:", condition))
    x <- tbl[!is.na(lactbin) & !is.na(get(condition)) & !is.na(f_gluc)]
    print(paste("Measures used", nrow(x)))

    if (!is.null(bins)) {
      x[[condition]] <- .bincode(x[[condition]], breaks = c(-Inf, bins, Inf))
    }

    print(paste(
      "p-value of the CI test",
      gcm.test(X = as.matrix(x[["f_gluc"]]), Y = as.matrix(x[[condition]]),
        Z = as.matrix(x[["lact"]]), regr.method = "gam")[["p.value"]]
    ))

    df <- x[!is.na(lactbin) & !is.na(get(condition)),
      median(f_gluc, na.rm = TRUE), by = c("lactbin", condition)]


    df[[condition]] <- as.factor(df[[condition]])
    if (!is.null(bins)) {
      levels(df[[condition]]) <- bin_labels(bins, unit)
    } else {
      levels(df[[condition]]) <- c("no", "yes")
    }

    data.table::setnames(df, condition, "mid")

    p <- ggplot(df, aes(x = lactbin, y = V1, color = mid)) +
      geom_line(size = 3) + theme_bw(15) + ggtitle(srcwrap(unique(tbl[["dataset"]]))) +
      ylab("Glucose (mg/dL)") + xlab("Lactate (mmol/L)") +
      scale_x_continuous(labels=bin_labels(lactate_bins, NULL), breaks = c(1:(length(lactate_bins)+1))) +
      guides(color=guide_legend(title=lwrap(condition)))
    if(grepl(legend.source, unique(tbl[["dataset"]]))) {
      p <- p + theme(legend.position = c(x.pos, y.pos),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_text(size=rel(legend.rel.size)))
    } else {
      p <- p + theme(legend.position = "none")
    }

    p
}

sens_spec_table <- function(score, outcome, src) {
  value_set <- sort(unique(score))
  sens <- spec <- NULL
  for(thresh in value_set) {
    pred <- as.integer(score >= thresh)
    sens <- c(sens, paste0(round(100*sum(pred == 1 & outcome == 1) / sum(outcome == 1), 2), "%"))
    spec <- c(spec, paste0(round(100*sum(pred == 0 & outcome == 0) / sum(outcome == 0), 2), "%"))
  }

  res <- as.data.frame(cbind(sens, spec), stringsAsFactors = F)
  res <- cbind(threshold = value_set, res)
  names(res) <- c("threshold", paste0(c("sens", "spec"), "_", src))

  res
}

tgtl <- function(source, patient_ids = cohort(source),
                 verbose = FALSE) {

  lg <- load_concepts(c("lact", "gluc"), source,
                      patient_ids = patient_ids, verbose = verbose)

  hypo.threshold <- 18.016 * 3.9
  lactate.threshold <- 2

  # make the table regular
  lg <- fill_gaps(lg)

  lg[, hypo := as.integer(glu < hypo.threshold)]
  lg[, lactatemia := as.integer(lact > lactate.threshold)]

  lg[is.na(hypo), "hypo"] <- 0
  lg[is.na(lactatemia), "lactatemia"] <- 0

  lg[, lactatemia := cummax(lactatemia), by = eval(id_vars(lg))]
  lg[, hypo := cummax(hypo), by = eval(id_vars(lg))]

  lg[, coinc := (max(lactatemia)+max(hypo)), by = eval(id_vars(lg))]
  lg <- lg[coinc == 2]

  diff <- lg[, sum(lactatemia - hypo), by = eval(id_vars(lg))]
  diff <- diff[abs(V1) <= 48]

  p <- ggplot(diff) + theme_minimal(10) + ylab("Number of cases") + xlab(TeX("$t_G - t_L$")) +
    geom_histogram(aes(V1), breaks = seq(-48, 48, by = 6), closed = "left", fill = "white", color = "black") +
    ggtitle(srcwrap(source))
  print("---------------------")
  print(srcwrap(source))
  test_stat <- 2*(mean((diff$V1 + rnorm(nrow(diff), sd = 0.01)) > 0) - 0.5) * sqrt(nrow(diff))
  print(paste("p-value of the test is ", 1 - pnorm(test_stat)))
  print(paste("25%, 50% and 75% quantiles of t_G-t_L are",
    paste(quantile(diff$V1, probs = c(0.25, 0.5, 0.75)), collapse = ", ")))

  print(sprintf("Hyperlactatemia precedes hypoglycemia %.2f%% of the time", (diff$V1 > 0) / (diff$V1 != 0)))

  print("---------------------")
  return(p)
}

tg_tl <- function(source, patient_ids = cohort(source),
                 verbose = FALSE) {

  lg <- load_concepts(c("lact", "glu"), source,
                      patient_ids = patient_ids, verbose = verbose)

  hypo.threshold <- 18.016 * 3.9
  lactate.threshold <- 2

  # make the table regular
  lg <- fill_gaps(lg)

  lg[, hypo := as.integer(glu < hypo.threshold)]
  lg[, lactatemia := as.integer(lact > lactate.threshold)]

  lg[is.na(hypo), "hypo"] <- 0
  lg[is.na(lactatemia), "lactatemia"] <- 0

  lg[, lactatemia := cummax(lactatemia), by = eval(id_vars(lg))]
  lg[, hypo := cummax(hypo), by = eval(id_vars(lg))]

  lg[, coinc := (max(lactatemia)+max(hypo)), by = eval(id_vars(lg))]
  lg <- lg[coinc == 2]

  diff <- lg[, sum(lactatemia - hypo), by = eval(id_vars(lg))]
  diff <- diff[abs(V1) <= 48]

  diff[["V1"]]
}

tw_glucose <- function(source, patient_ids = cohort(source)) {

  x <- fill_gaps(load_concepts("glu", source, patient_ids = patient_ids, verbose = F))
  x[, glu := data.table::nafill(glu, "locf"), by = eval(id_var(x))]

  wins <- stay_windows(source)
  x <- merge(x, wins)
  x[get(index_var(x)) >= start & get(index_var(x)) <= end]

  x[, mean(glu, na.rm = T), by = eval(id_var(x))][["V1"]]
}

mean_glucose <- function(source, patient_ids = cohort(source)) {

  x <- fill_gaps(load_concepts("glu", source, patient_ids = patient_ids, verbose = F))
  #x[, glu := data.table::nafill(glu, "locf"), by = eval(id_var(x))]

  wins <- stay_windows(source)
  x <- merge(x, wins)
  x[get(index_var(x)) >= start & get(index_var(x)) <= end]

  x[, mean(glu, na.rm = T), by = eval(id_var(x))][["V1"]]
}

insulin_days <- function(source, patient_ids = cohort(source), upto = hours(10*24)) {

  wins <- stay_windows(source)
  patient_ids <- intersect(id_col(wins[!is.na(end)]), patient_ids)
  wins <- wins[get(id_var(wins)) %in% patient_ids]
  x <- load_concepts("ins", source, patient_ids = patient_ids, verbose = F)


  wins[, end := hours(24*round(end/24))] # round the stay to days
  wins[end > upto, "end"] <- upto
  wins[, num_days := as.integer(end/24)]

  x <- merge(x, wins, all = T)
  x <- x[get(index_var(x)) >= start & get(index_var(x)) < end]

  num_days <- sum(x[, length(unique(.bincode(get(index_var(x)), seq(-0.01, as.integer(upto)-24, 24)))), by = eval(id_var(x))][["V1"]])

  res <- rep(FALSE, sum(wins[["num_days"]]))
  res[1:num_days] <- TRUE

  res
}

is_alarm <- function(tim, pos, diff) {

  opt <- tim[pos]

  while (length(opt)) {

    now <- min(opt)
    cut <- now + diff
    opt <- opt[opt >= cut]

    pos[tim > now & tim < cut] <- NA
  }

  pos
}

not_null <- Negate(is.null)

reorder_cols <- function(...) data.table::setcolorder(...)

n_cores <- function(verbose = FALSE) {

  res <- as.integer(
    Sys.getenv("LSB_DJOB_NUMPROC", unset = parallel::detectCores() / 2L)
  )

  if (isTRUE(verbose)) {
    msg_ts("unsing ", res, " cores")
  }

  res
}

jobid <- function() {
  Sys.getenv(
    "LSB_JOBID",
    unset = substr(as.character(as.integer(Sys.time())), 1, 9)
  )
}

jobname <- function() {
  sub("\\[.+\\]$", "", Sys.getenv("LSB_JOBNAME", unset = "test"))
}

xtr_scalar <- function(type) function(x, i) vapply(x, `[[`, type, i)

xtr_num <- xtr_scalar(numeric(1L))
xtr_chr <- xtr_scalar(character(1L))

msg_ts <- function(...) {
  mem <- vapply(memuse::Sys.procmem(), as.character, character(1L))
  mem <- paste0(names(mem), ": ", mem, collapse = "; ")
  message(Sys.time(), " [",  mem, "] ", ...)
}

is_lsf <- function() !is.na(Sys.getenv("LSB_JOBID", unset = NA))

set_names <- function(x, val) `names<-`(x, val)

unique_name <- function() {
  paste0("run_", Sys.getenv("LSB_JOBINDEX", unset = rand_name()))
}

rand_name <- function(n = 1L, length = 12L, chars = c(letters, 0:9)) {

  res <- character()

  while (length(unique(res)) < n) {
    res <- replicate(n,
      paste(sample(chars, length, replace = TRUE), collapse = "")
    )
  }

  res
}

fmt4 <- function(x) format(x, trim = TRUE, digits = 4)

fmt_perc <- function(x, collapse = "/") {
  paste0(paste(fmt4(x * 100), collapse = collapse), "%")
}

rename <- function(x, new, old = names(x)) {

  assert_that(length(new) == length(old), all(old %in% names(x)),
              anyDuplicated(new) == 0, anyDuplicated(old) == 0)

  names(x)[match(old, names(x))] <- new

  x
}

update_controller <- function(x, ..., lst = NULL) {

  set <- c(list(...), lst)

  Map(assign, names(set), set, MoreArgs = list(envir = x))

  invisible(NULL)
}

call_do <- function(args, what) do.call(what, args)

create_dir <- function(dir) dir[dir.create(dir)]
