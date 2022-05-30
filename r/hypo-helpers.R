
liver_damage <- function(data_source = "mimic", id_type = "icustay",
                         patient_ids = NULL, verbose = FALSE){

  liv_dam <- load_concepts(
    c("bili", "alt", "ast"), src = data_source,
    id_type = id_type, patient_ids = patient_ids, verbose = verbose
  )
  liv_dam[, liver_damage := as.integer(bili > 2 | alt > 45 | ast > 45)]
  liv_dam <- liv_dam[, c(meta_vars(liv_dam), "liver_damage"), with = FALSE]

  liv_dam
}

shock <- function(data_source = "mimic", id_type = "icustay",
                  patient_ids = NULL, verbose = FALSE) {

  x <- load_concepts(c("dopa_rate", "norepi_rate", "dobu_rate", "epi_rate", "map"),
    data_source, id_type = id_type, patient_ids = patient_ids, verbose = verbose)

  if (data_source == "hirid") x[, dopa := NA]

  x <- x[, c(meta_vars(x), "dopa_rate", "norepi_rate", "dobu_rate", "epi_rate", "map"), with = FALSE]
  x[is.na(map), "map"] <- 100

  shock <- x[, as.integer(any(!is.na(dopa_rate), !is.na(dobu_rate), !is.na(epi_rate), !is.na(norepi_rate), map < 60)), by = eval(meta_vars(x))]
  shock <- data.table::setnames(shock, "V1", "shock")
  shock <- shock[shock == 1]

  return(shock)
}

mech_vent <- function(data_source = "mimic", id_type = "icustay",
                      patient_ids = NULL, verbose = FALSE) {

  mechv <- ricu:::sofa_vent(
    load_concepts("vent_start", data_source, id_type = id_type,
                  patient_ids = patient_ids, verbose = verbose),
    load_concepts("vent_end", data_source, id_type = id_type,
                  patient_ids = patient_ids, verbose = verbose),
    hours(6L), hours(2L), hours(1L)
  )

  mechv[, mech_vent := as.integer(vent)]
  mechv <- mechv[, c(meta_vars(mechv), "mech_vent"), with = FALSE]

  return(mechv)
}

load_ood <- function(data_source, concepts, id_type = "icustay", patient_ids = NULL) {

  res <- list()

  for(i in 1:length(concepts)) {

    f <- eval(parse(text = (concepts[i])))
    res[[i]] <- f(data_source, id_type = id_type, patient_ids = patient_ids)

  }

  if (length(res) == 1) return(res[[1]])

  res <- Reduce(function(x, y) merge(x, y, all = TRUE), res)

  return(res)

}

collect_hypo_cases <- function(tbl, max.hour = 240L) {
  # carry the hypo time backwards for 6 hours
  tbl <- slide(tbl, before = hours(0L), after = hours(5L), hypo_LA := max(hypo))

  # collect the relevant rows
  marks <- 6 * (1:(max.hour/6))
  rel_cols <- c(setdiff(names(tbl), meta_vars(tbl)))
  collect <- NULL
  for (mark in marks) {
    tmp <- tbl[get(index_var(tbl)) == mark]
    collect <- rbind(
      collect,
      as.matrix(tmp[, rel_cols, with = FALSE])
    )
  }

  collect <- data.table::data.table(collect)
  collect[, hypo := NULL]
  collect <- data.table::setnames(collect, "hypo_LA", "hypo")

  collect

}

glycemia_treatment <- function(data_source,
  vars = list(glu = list(time = 24L, imp_val = NA_real_),
    lact = list(time = 24L, imp_val = 1), ins = list(time = 12L, imp_val = 0),
    shock = list(time = 24L, imp_val = 0), mech_vent = list(time = 24L, imp_val = 0),
    liver_damage = list(time = 48L, imp_val = 0)), fill_na = FALSE,
    id_type = "icustay", patient_ids = NULL, hypo = TRUE, hypo.threshold = 3.9,
    verbose = FALSE) {

  dict <- ricu::get_config("concept-dict")

  in_dict <- intersect(names(vars), names(dict))
  out_dict <- setdiff(names(vars), in_dict)
  tbl1 <- load_concepts(in_dict, data_source, id_type = id_type,
                        patient_ids = patient_ids, verbose = verbose)

  if (length(out_dict) > 0L) {
    tbl2 <- load_ood(data_source, out_dict, id_type = id_type, patient_ids = patient_ids)
    tbl <- merge(tbl1, tbl2, all = TRUE)
  } else {
    tbl <- tbl1
  }

  if (data_source == "mimic" & is.element("ins", names(tbl))) tbl[ins == 0, "ins"] <- 2 # MIMIC carevue imputation
  # reorder the columns appropriately
  tbl <- tbl[, c(meta_vars(tbl), names(vars)), with = FALSE]

  # fill gaps
  tbl <- fill_gaps(tbl)
  tbl <- carry_values(tbl, vars)
  #fwd_times <- lapply(vars, function(x) x[1])
  #tbl <- slide_quo(tbl, before = hours(48L), substitute(Map(carry_fwd, .SD, fwd_times), list(fwd_times = fwd_times)))

  # fill NAs after the carry-forward if specified
  na_vals <- lapply(vars, function(x) x[["imp_val"]])
  if (fill_na) tbl <- tbl[, c(mget(meta_vars(tbl)), Map(repl_na, .SD, na_vals)), .SDcols = names(vars)]


  if(hypo) {
    # determine the first hypo time (glucose < 3.9 mmol/L)
    tbl[, hypo := as.integer(glu < hypo.threshold*18.016)]
    tbl[is.na(hypo), "hypo"] <- 0

    # delete everything after the first hypo time
    tbl <- tbl[, head(.SD, n = min(which(hypo == 1), length(hypo))), by = eval(id_vars(tbl))]
  }

  return(tbl)
}

lga <- function(g, liver, shock, ins, map) {

  if(all(is.na(g))) return(as.list(rep(NA_real_, 5)))

  upto <- min(which(!is.na(g)))
  range <- seq.int(1L, upto)

  list(
    g[upto], max_or_na(liver[range]), max_or_na(shock[range]),
    (sum(ins[range], na.rm = T) > 0)*1L, min_or_na(map[range])
  )
}

continuous_effects <- function(concepts, source, dir = "increasing", breaks,
  patient_ids = NULL, upto = hours(24L), soffa, sofa_breaks, condition = "feat", thresholds, imp_val) {

  hypo <- hypo(source, patient_ids = patient_ids)
  soffa[, sofa_adjust := (sofa - sofa_liver_comp)]
  soffa <- soffa[get(index_var(soffa)) <= upto]
  soffa[, w_sofa := max_or_na(sofa_adjust), by = eval(id_vars(soffa))]
  soffa <- unique(soffa[, c(id_vars(soffa), "w_sofa"), with = FALSE])
  soffa[, sofa_bins := .bincode(w_sofa, c(-Inf, sofa_breaks, Inf))]

  vals <- lapply(concepts,
    function(c) w_value(source, c, dir = dir, upto = upto,
                        patient_ids = patient_ids, imp_val = imp_val))

  vals <- lapply(vals, function(x) merge(x, hypo, by = id_vars(x), all.x = T))

  lapply(vals, function(x) x[, bins := .bincode(w_val, c(-Inf, breaks, Inf))])

  vals <- lapply(vals, function(x) merge(x, soffa, by = id_vars(x), all.x = T))
  vals <- lapply(vals, function(x) {
    x[!is.na(w_sofa) & !is.na(w_val)]
  })

  unadjusted_effect <- lapply(thresholds,
    function(t, boot_fun) {
      boot_out <- lapply(vals, function(x) boot(x, boot_fun, R = 100,
                                                threshold = t))
      boot_ci <- lapply(boot_out, function(boo) boot.ci(boo, type = "perc"))
      fin_ci <- lapply(boot_ci, function(ci) {c(ci$t0, ci$percent[4:5])})
      fin_ci[[1]]
    }, boot_fun = diff_means
  )

  unadjusted_effect <- data.frame(Reduce(rbind, unadjusted_effect))
  unadjusted_effect <- cbind(unadjusted_effect, thresholds, "unadjusted")

  adjusted_effect <- lapply(thresholds,
    function(t, boot_fun) {
      boot_out <- lapply(vals, function(x) boot(x, boot_fun, R = 100,
                                                threshold = t))
      boot_ci <- lapply(boot_out, function(boo) boot.ci(boo, type = "perc"))
      fin_ci <- lapply(boot_ci, function(ci) c(ci$t0, ci$percent[4:5]))
      fin_ci[[1]]
    }, boot_fun = ate_sa
  )

  adjusted_effect <- data.frame(Reduce(rbind, adjusted_effect))
  adjusted_effect <- cbind(adjusted_effect, thresholds, "SOFA adjusted")

  names(adjusted_effect) <- names(unadjusted_effect) <- c("y", "ymin", "ymax", "threshold", "Type")

  p <- ggplot(rbind(adjusted_effect, unadjusted_effect),
              aes(x = threshold, y = y, color = Type)) +
    geom_line(size = 3) + theme_bw(15) + ggtitle(srcwrap(source)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Type), alpha = 0.3,
                linetype = 0) +
    ylab("Estimated effect") + xlab("AST (IU/L)")
  if(grepl("mimic", source)) {
    p <- p + theme(legend.position = c(0.7, 0.1375),
                   legend.box.background = element_rect("black"),
      legend.title=element_text(size=rel(0.75)))
  } else {
    p <- p + theme(legend.position = "none")
  }

  p
}

RL_summary <- function(oos) {
  
  ncoh <- length(unique(id_col(oos[Q1 > 0])))
  smr <- oos[, 
             list(
               t0 = as.integer(head(get(index_var(oos))[Q1 > 0], n = 1L)),
               t1 = as.integer(head(get(index_var(oos))[Q1 > Q0 & Q1 > 0], n = 1L)),
               th = as.integer(tail(get(index_var(oos))[is_true(hypo)], n = 1L))
             ), 
             by = c(id_vars(oos))
  ]
  
  prop <- sum(is_true(smr$t0 == smr$t1)) / ncoh
  cat(round(100 * prop, 2), "% alarms unchanged \n")
  
  del <- smr[!is.na(t0) & !(t0 == t1)]
  
  cat("Of remaining", round(100 * (1-prop), 2), "% alarms:\n\n")
  
  fa <- round(100 * nrow(del[is.na(th)]) / nrow(del), 2)
  fa_save <- round(100 + del[is.na(th) & is.na(t1)] / nrow(del[is.na(th)]), 2)
  cat(fa, "% are false alarms, of which", fa_save, "% are saved by Pi_1\n")
  
  tim <- round(100 * nrow(del[!is.na(th) & (th - t0) <= 24L]) / nrow(del), 2)
  tim_ruin <- round(100 * del[!is.na(th) & (th - t0) <= 24L & is.na(t1)] / 
                      nrow(del[!is.na(th) & (th - t0) <= 24L]), 2)
  
  cat(tim, "% are timely alarms, of which", tim_ruin, "% are ruined by Pi_1\n")
  
  ear <- round(100 * nrow(del[!is.na(th) & (th - t0) > 24L]) / nrow(del), 2)
  ear_save <- round(100 * del[!is.na(th) & (th - t0) > 24L & (th - t1) <= 24L] / 
                      nrow(del[!is.na(th) & (th - t0) > 24L]), 2)
  
  cat(ear, "% are early alarms, of which", ear_save, "% are ruined by Pi_1\n")
}

