
#' @export
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

#' @export
ins_ifx_cb <- function(ins, interval) {
  if (id_var(ins) == "icustay_id") ins[ins == 0, "ins"] <- 2
  rename_cols(ins, "ins_ifx", "ins")
}

#' @export
hypo_cb <- function(glu, ...) {

  onset_id <- function(x) replace(x, x, seq_len(sum(x)))

  glu_var <- data_var(glu)
  id_vars <- id_vars(glu)

  glu <- glu[, c("tmp") := is_true(get(glu_var) <= 3.9 * 18.016)]
  glu <- glu[, c("hypo") := onset_id(get("tmp")), by = c(id_vars)]
  glu <- glu[, c(glu_var, "tmp") := NULL]

  glu
}

#' @export
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

#' @export
hypo_cnt_cb <- function(hypo_epi, ...) {
  
  hypo_epi <- fill_gaps(hypo_epi)
  hypo_epi <- hypo_epi[, lapply(.SD, nafill, type = "locf"), 
                       by = c(id_vars(hypo_epi))]
  hypo_epi <- replace_na(hypo_epi, 0)
  hypo_epi[, hypo_cnt := cummax(hypo_epi), by = c(id_vars(hypo_epi))]
  
  hypo_epi[, c(meta_vars(hypo_epi), "hypo_cnt"), with=FALSE]
  
}

#' @export
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

#' @export
ts_to_win_2hours <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(120L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

#' @export
ts_to_win_6hours <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

#' @export
mimic_presc_cort <- function(x, dur_var, ...) {
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(index_var(x)) := get(index_var(x)) + hours(9L)]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

#' @export
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

#' @export
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

#' @export
aumc_cortico <- function(x, dur_var, ...) {
  x[, c(dur_var) := get(dur_var) + mins(360L)]
  x
}

#' @export
sep3_info <- function(source, pids = NULL, keep_components = TRUE,
                         dat = NULL) {
  
  if (is.null(dat)) {
    dat <- load_concepts("sofa", source, patient_ids = pids,
                         keep_components = keep_components, verbose = FALSE)
  } else if (!is_ts_tbl(dat)) {
    dat <- data.table::copy(dat[["sofa"]])
  }
  
  stopifnot(is_ts_tbl(dat))
  
  if (grepl("eicu", source) | identical(source, "hirid")) {
    
    si <- load_concepts("susp_inf", source, abx_min_count = 2L,
                        id_type = "icustay", patient_ids = pids,
                        si_mode = "or", keep_components = keep_components,
                        verbose = FALSE)[susp_inf == TRUE]
  } else {
    
    si <- load_concepts("susp_inf", source, id_type = "icustay",
                        patient_ids = pids, keep_components = keep_components,
                        verbose = FALSE)
  }
  
  sep3_cmp(dat, si, si_window = "any", keep_components = keep_components,
           source = source)
}

sep3_cmp <- function (..., si_window = c("first", "last", "any"), 
                      delta_fun = delta_cummin, sofa_thresh = 2L, 
                      si_lwr = hours(48L), si_upr = hours(24L), 
                      keep_components = TRUE, interval = NULL, source = NULL) 
{
  cnc <- c("sofa", "susp_inf")
  res <- ricu:::collect_dots(cnc, interval, ...)
  assert_that(is.count(sofa_thresh), is.flag(keep_components), 
              ricu:::not_null(delta_fun))
  si_lwr <- ricu:::as_interval(si_lwr)
  si_upr <- ricu:::as_interval(si_upr)
  delta_fun <- ricu:::str_to_fun(delta_fun)
  si_window <- match.arg(si_window)
  sofa <- res[["sofa"]]
  susp <- res[["susp_inf"]]
  id <- id_vars(sofa)
  ind <- index_var(sofa)
  sus_cols <- setdiff(data_vars(susp), "susp_inf")
  sofa <- sofa[, `:=`(c("join_time1", "join_time2"), list(get(ind), 
                                                          get(ind)))]
  on.exit(rm_cols(sofa, c("join_time1", "join_time2"), by_ref = TRUE))
  susp <- susp[is_true(get("susp_inf")), ]
  susp <- susp[, `:=`(c("susp_inf"), NULL)]
  susp <- susp[, `:=`(c("si_lwr", "si_upr"), 
                      list(get(index_var(susp)) - si_lwr, 
                           get(index_var(susp)) + si_upr))
  ]
  if (si_window %in% c("first", "last")) {
    susp <- dt_gforce(susp, si_window, id)
  }
  join_clause <- c(id, "join_time1 >= si_lwr", "join_time2 <= si_upr")
  res <- sofa[susp, c(list(delta_sofa = delta_fun(get("sofa"))), 
                      mget(c(ind, sus_cols))), on = join_clause, by = .EACHI, 
              nomatch = NULL]
  res <- res[is_true(get("delta_sofa") >= sofa_thresh), ]
  cols_rm <- c("join_time1", "join_time2")
  if (!keep_components) {
    cols_rm <- c(cols_rm, "delta_sofa")
  }
  res <- rm_cols(res, cols_rm, by_ref = TRUE)
  res <- res[, head(.SD, n = 1L), by = c(id_vars(res))]
  res <- res[, `:=`(c("sep3"), TRUE)]

  if (!(grepl("eicu", source) | source == "hirid")) {
    
    res[, t_confirm := pmax(get(index_var(res)), abx_time, samp_time)]
    
  } else {
    
    # solve the sampling case
    res_samp <- res[!is.na(samp_time)]
    res_samp[, t_confirm := pmax(get(index_var(res_samp)), abx_time, samp_time,
                            na.rm = TRUE)]
    res_samp <- res_samp[, c(meta_vars(res_samp), "sep3", "t_confirm"),
                         with=FALSE]
    
    # go to abx case
    res_abx <- res
    res_abx <- res_abx[, c(id_var(res_abx), "abx_time"), with=F]
    add_abx <- load_concepts("abx", source, aggregate = "sum", verbose = FALSE)
    
    res_abx <- merge(add_abx, res_abx)
    res_abx <- res_abx[(get(index_var(res_abx)) > abx_time) |
                         (get(index_var(res_abx)) >= abx_time & abx > 1),
                       head(.SD, n = 1L),
                       by = c(id_vars(res_abx))]
    
    res_abx <- rename_cols(res_abx, "later_abx", index_var(res_abx))
    res_abx <- as_id_tbl(res_abx)
    res_abx <- res_abx[, c(id_vars(res_abx), "later_abx"), with = FALSE]
    
    res_abx <- merge(res[is.na(samp_time)], res_abx)
    res_abx[, t_confirm := pmax(get(index_var(res_abx)), abx_time, later_abx)]
    res_abx <- res_abx[, c(meta_vars(res_abx), "sep3", "t_confirm"), with=FALSE]
    
    res <- rbind(res_samp, res_abx) 
    
  }
  res <- res[, c(meta_vars(res), "sep3", "t_confirm"), with=FALSE]
  res <- rename_cols(res, "stay_id", id_var(res))
  res
}
