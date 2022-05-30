
patient_eval <- function(eval, delta, type = "circEWS") {
  
  if (!is.element("method", names(eval))) eval$method <- "RF"
  
  if (length(unique(eval[["tst_src"]])) > 1L | 
      length(unique(eval[["method"]])) > 1L ) {

    return(
      Reduce(
        rbind,
        lapply(split(eval, by = c("source", "method"), flatten = TRUE),
               patient_eval, delta, type)
      )
    )
    
  }
  
  
  if (length(delta) > 1L) {
    
    return(Reduce(rbind, lapply(delta, function(d) patient_eval(eval, d, type))))
    
  }
  
  convex_comb <- function(a, b, val1, val2, mid = 0.8) {
    
    assert_that(a >= mid, b <= mid)
    
    val2 + (val1 - val2) * (mid - b) / (a - b)
  }
  
  tpp <- "label"
  score_col <- "prediction"
  
  if (length(unique(eval[[score_col]])) < 50) {
    
    thresh <- sort(unique(eval[[score_col]]))
    
  } else {
    
    thresh <- unique(
      quantile(eval[, max(get(score_col)), by = "icu_stay"][["V1"]],
               prob = seq(0.01, 0.99, 0.01))
    )
    
  }
  
  res <- c(-1000, 1, 0, mean(eval[[tpp]]))
  res <- rbind(res, c(tail(thresh, n = 1L), 0, 1, 1))
  
  eval <- eval[, time_dbl := as.double(get(index_var(eval)))]
  lb_win <- as.double(delta, units = time_unit(eval))
  
  for (t in thresh[1:(length(thresh)-1)]) {
    
    eval[, pos := (get(score_col) > t)]
    eval[, alarm := is_alarm(time_dbl, pos, diff = lb_win),
         by = c(id_vars(eval))]
    
    dcs <- table(eval[!is.na(alarm)][[tpp]],
                 eval[!is.na(alarm)][["alarm"]])
    
    if (all(eval[!is.na(alarm)][["alarm"]])) {
      tn <- fn <- 0L
    } else {
      tn <- dcs["FALSE", "FALSE"]
      fn <- dcs["TRUE", "FALSE"]
    }
      
    fp <- dcs["FALSE", "TRUE"]
    tp <- dcs["TRUE", "TRUE"]
    
    if (is.element(type, c("circEWS", "honest"))) {
      sens <- eval[, list(is_hit = any(is_true(label) & is_true(alarm)), 
                          is_case = any(is_true(label))), by = "icu_stay"]
      sens <- mean(sens[is_true(is_case)][["is_hit"]])
    } else sens <- tp / (tp + fn)
    
    if (type == "honest") {
      spec <- eval[, list(is_falrm = any(is_false(label) & is_true(alarm)), 
                          is_control = all(is_false(label))), by = "icu_stay"]
      spec <- 1 - mean(spec[is_true(is_control)][["is_falrm"]])
    } else spec <- tn / (tn + fp)
    
    ppv  <- tp / (tp + fp)
    
    res <- rbind(res, c(t, sens, spec, ppv))
    
  }
  
  res <- as.data.frame(res)
  names(res) <- c("thresh", "sens", "spec", "ppv")
  res <- res[order(res$thresh), ]
  
  for (i in seq_len(nrow(res))) {
    if (res$sens[i] > 0.8 & res$sens[i+1] <= 0.8) {
      
      prec_80r <- convex_comb(
        res$sens[i], res$sens[i+1], res$ppv[i], res$ppv[i+1]
      )
    }
  }

  res <- cbind(res, Dataset = unique(eval[["tst_src"]]), 
               Delta = as.double(delta), Method = unique(eval[["method"]]), 
               prec80 = prec_80r)
  
  res
}

make_plot <- function(res) {

  if (length(unique(res$Dataset)) > 1L) {
    
    return(
      plot_grid(
        plotlist = lapply(split(res, res$Dataset), make_plot), ncol = 1L
      )
    )
    
  }
  
  integrate <- function(x, y) {
    assert_that(length(x) == length(y))
    
    len1 <- seq_len(length(x) - 1)
    len2 <- seq.int(2L, length(x))
    
    sum((x[len2] - x[len1]) * (y[len1] + y[len2]) / 2)
  }
  
  res <- data.table::as.data.table(res)
  
  res[, auroc := -integrate(1 - spec, sens), by = c("Dataset", "Delta", "Method")]
  res[, auprc := -integrate(sens, ppv), by = c("Dataset", "Delta", "Method")]
  
  res[, Delta_auroc := paste0(Delta, " (", round(auroc, 3), ")")]
  res[, Delta_auprc := paste0(Delta, " (", round(auprc, 3), ")")]
  res[, method_prec := paste0(Method, " (", 
                              round(prec_shift_module(prec80, res[1][["ppv"]], 
                                                      0.1), 3), ")")]
  
  roc <- ggplot(res, aes(x = 1-spec, y = sens, linetype = factor(Delta),
                         color = factor(Method))) +
    geom_line() +
    theme_bw() + ylim(c(0, 1)) + xlim(c(0,1)) +
    #facet_grid(rows = vars(srcwrap(Dataset))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotdash") +
    ggtitle(paste0("ROC on ", srcwrap(res$Dataset))) + 
    scale_color_brewer(palette="Dark2") +
    xlab("1 - Specificity") + ylab("Sensitivity") +
    theme(
      legend.position = c(0.7, 0.3),
      legend.box.background = element_rect(color = "black"),
      legend.background = element_blank()
    ) + labs(color = "Method", linetype = "Delta (hours)")
  
  prc <- ggplot(res, aes(x = sens, y = ppv, linetype = factor(Delta),
                         color = factor(method_prec))) +
    geom_line() +
    theme_bw() + ylim(c(0, 1)) +
    #facet_grid(rows = vars(srcwrap(Dataset))) +
    geom_hline(yintercept = res[1][["ppv"]], linetype = "dotdash") +
    ggtitle(paste0("PRC on ", srcwrap(res$Dataset), " with base PPV ", 
                   round(res[1][["ppv"]], 3))) + 
    scale_color_brewer(palette="Dark2") +
    ylab("PPV") + xlab("Sensitivity") +
    theme(
      legend.position = c(0.7, 0.7),
      legend.box.background = element_rect(color = "black"),
      legend.background = element_blank()
    ) + labs(color = "Method", linetype = "Delta (hours)")
  
  cowplot::plot_grid(roc, prc, ncol = 2L)
  
}

prec_shift_module <- function(prec, prev, target_prev) {
  s_1s <- prec / (1 - prec) * (1 - prev) / prev
  target_s_1s <- s_1s *  target_prev / (1 - target_prev)
  target_s_1s / (1 + target_s_1s)
}

select_job <- function(prefix = "test", jobid = NULL,
                       res_dir = results_dir()) {

  job_dir <- paste(prefix, jobid, sep = "_")

  cands <- list.dirs(res_dir)
  cands <- cands[grepl(job_dir, basename(cands))]

  job_dir <- cands[
    which.max(as.integer(xtr_chr(strsplit(basename(cands), "_"), 2L)))
  ]

  assert_that(length(job_dir) == 1L)

  job_dir
}

read_set <- function(run_name, dir = select_job()) {

  assert_that(is.string(run_name))

  if (!grep("\\.set", run_name)) {
    run_name <- paste0(run_name, ".set")
  }

  file <- file.path(dir, run_name)

  assert_that(file.exists(file))

  qs::qread(file)
}

read_sets <- function(run_names = NULL, job_prefix = "test",
                      dir = select_job(job_prefix)) {

  if (is.null(run_names)) {
    run_names <- list.files(dir, pattern = "\\.set")
  }

  lapply(run_names, read_set, dir)
}

read_pred <- function(run_name, tst_src, meta = list(run_name = run_name),
                      dir = select_job()) {

  assert_that(is.string(run_name), is.string(tst_src))

  file <- file.path(dir, paste0(run_name, "-", tst_src, ".prd"))

  assert_that(file.exists(file))

  res <- qs::qread(file)
  res <- res[, c("tst_src", names(meta)) := c(list(tst_src), meta)]

  res
}

select_sets <- function(sets, ...) {

  args <- list(...)
  hits <- rep(TRUE, length(sets))

  for (arg in names(args)) {
    temp <- lapply(sets, `[[`, arg)
    temp <- vapply(temp, `%in%`, logical(1L), args[[arg]])
    hits <- hits & temp
  }

  sets <- sets[hits]
  names(sets) <- xtr_chr(sets, "run_name")

  lapply(sets, `[`, names(args))
}

read_preds <- function(trn_src, tst_src = trn_src, ..., job_prefix = "test",
                       dir = select_job(job_prefix)) {

  ipply <- function(i, fun, x) lapply(x, fun, i)

  sets <- select_sets(read_sets(dir = dir), trn_src = trn_src, ...)

  prms <- expand.grid(run_name = names(sets), tst_src = tst_src,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  prms <- split(prms, seq_len(nrow(prms)))
  prms <- Map(c, prms, meta = lapply(sets, list), dir = list(dir))

  rbind_lst(lapply(prms, call_do, read_pred))
}

read_mod <- function(run_name, dir = select_job()) {

  assert_that(is.string(run_name))

  file <- file.path(dir, paste0(run_name, ".mod"))

  assert_that(file.exists(file))

  qs::qread(file)
}

read_mods <- function(..., job_prefix = "test", dir = select_job(job_prefix)) {

  set <- select_sets(read_sets(dir = dir), ...)
  res <- lapply(names(set), read_mod, dir)

  names(res) <- set

  res
}
