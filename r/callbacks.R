
map_beta <- function(beta) {

  assert_that(is.numeric(beta), is.scalar(beta))

  function (..., match_win = hours(2L), interval = NULL) {

    cnc <- c("map", "norepi_equiv")
    res <- ricu:::collect_dots(cnc, interval, ...)

    assert_that(ricu:::is_interval(match_win),
                match_win > ricu:::check_interval(res))

    on12 <- paste(meta_vars(res[[1L]]), "==", meta_vars(res[[2L]]))
    on21 <- paste(meta_vars(res[[2L]]), "==", meta_vars(res[[1L]]))

    res <- rbind(res[[1L]][res[[2L]], on = on12, roll = match_win],
                 res[[2L]][res[[1L]], on = on21, roll = match_win])

    res <- unique(res)

    res <- res[is.na(get(cnc[2L])), c(cnc[2L]) := 0]
    res <- res[!is.na(get(cnc[1L])), ]

    res <- res[,
      c(paste0("map_beta", beta)) := get(cnc[1L]) - beta * get(cnc[2L])
    ]

    res <- rm_cols(res, cnc, by_ref = TRUE)

    res
  }
}

ins_ifx_cb <- function(ins, interval) {
  if (id_var(ins) == "icustay_id") ins[ins == 0, "ins"] <- 2
  rename_cols(ins, "ins_ifx", "ins")
}

hypo_cb <- function(glu, ...) {

  onset_id <- function(x) replace(x, x, seq_len(sum(x)))

  glu_var <- data_var(glu)
  id_vars <- id_vars(glu)

  glu <- glu[, c("tmp") := is_true(get(glu_var) <= 3.9 * 18.016)]
  glu <- glu[, c("hypo") := onset_id(get("tmp")), by = c(id_vars)]
  glu <- glu[, c(glu_var, "tmp") := NULL]

  glu
}

hypo_epd <- function(hypo, min_dur = hours(6L), ...) {
  
  mrg <- function(x) cumsum(c(TRUE, x[-length(x)]))
  
  trm <- function(x, mx) {
    ind <- max(which(x > 0L))
    sft <- ind + mx
    if (sft <= length(x)) replace(x, sft, -x[ind]) else x
  }
  
  idv <- id_vars(hypo)
  idx <- index_var(hypo)
  
  assert_that(!any(c("diff", "hdif", "nhyp", "impu") %in% colnames(hypo)))
  
  hypo <- hypo[!is.na(hypo), diff := c(diff(get(idx)), Inf), by = c(idv)]
  hypo <- hypo[is_true(hypo > 0L), hdif := c(diff(get(idx)), Inf), by = c(idv)]
  hypo <- hypo[!is.na(hdif), nhyp := .N, by = c(idv)]
  
  hypo <- hypo[is_true(nhyp > 1L), hypo := mrg(hdif > min_dur | diff < hdif),
               by = c(idv)]
  hypo <- hypo[, c("diff", "hdif", "nhyp") := NULL]
  hypo <- hypo[, hypo_epi := data.table::nafill(hypo, fill = 0)]
  
  hypo
}

hypo_cnt_cb <- function(hypo_epi, ...) {
  
  hypo_epi <- fill_gaps(hypo_epi)
  hypo_epi <- hypo_epi[, lapply(.SD, nafill, type = "locf"), 
                       by = c(id_vars(hypo_epi))]
  hypo_epi <- replace_na(hypo_epi, 0)
  hypo_epi[, hypo_cnt := cummax(hypo_epi), by = c(id_vars(hypo_epi))]
  
  hypo_epi[, c(meta_vars(hypo_epi), "hypo_cnt"), with=FALSE]
  
}

dex_amount_callback <- function(...) {
  x <- list(...)[["dex"]]
  ivl <- list(...)[["interval"]]
  c_fct <- as.double(ivl, units = units(x$dur_var))
  # make into amount
  x[, dex := dex * as.double(dur_var) / c_fct]
  x[, distr := as.double(ricu:::re_time(dur_var, ivl)) + 1]
  x[, dex := dex / distr]
  x[, distr := NULL]
  x <- rename_cols(x, "dex_amount", "dex")
  
  expand(x, aggregate = "sum")
}

ts_to_win_2hours <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(120L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

ts_to_win_6hours <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

mimic_presc_cort <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(index_var(x)) := get(index_var(x)) + hours(9L)]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

hirid_pharma_win6 <- function(x, dur_var, group_var, ...) {
  
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  
  x[, c(dur_var) := max(get(index_var(x))) - min(get(index_var(x))) + hours(6L), 
    by = c(group_var)]
  x <- x[, head(.SD, n = 1L), by = c(group_var)]
  x[, c(dur_var) := `units<-`(get(dur_var), "mins")]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

hirid_pharma_win2 <- function(x, dur_var, group_var, ...) {
  
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  
  x[, c(dur_var) := max(get(index_var(x))) - min(get(index_var(x))) + hours(2L), 
    by = c(group_var)]
  x <- x[, head(.SD, n = 1L), by = c(group_var)]
  x[, c(dur_var) := `units<-`(get(dur_var), "mins")]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

aumc_cortico <- function(x, dur_var, ...) {
  
  x[, c(dur_var) := get(dur_var) + mins(360L)]
  
}