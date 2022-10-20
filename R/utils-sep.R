
#' Sepsis utilities
#' @param source Data source name
#' @param pids Patient IDs
#' @param keep_components Flag indicating whether to return individual score
#' components
#' @param dat Data with SOFA score (computed if `NULL`)
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
                      keep_components = TRUE, interval = NULL, source = NULL) {

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
    
    res[, c("t_confirm") := pmax(
      get(index_var(res)), get("abx_time"), get("samp_time")
    )]
    
  } else {
    
    # solve the sampling case
    res_samp <- res[!is.na(get("samp_time"))]
    res_samp[, c("t_confirm") := pmax(
      get(index_var(res_samp)), get("abx_time"), get("samp_time"), na.rm = TRUE
    )]
    res_samp <- res_samp[, c(meta_vars(res_samp), "sep3", "t_confirm"),
                         with=FALSE]
    
    # go to abx case
    res_abx <- res
    res_abx <- res_abx[, c(id_var(res_abx), "abx_time"), with=F]
    add_abx <- load_concepts("abx", source, aggregate = "sum", verbose = FALSE)
    
    res_abx <- merge(add_abx, res_abx)
    res_abx <- res_abx[
      (get(index_var(res_abx)) > get("abx_time")) |
        (get(index_var(res_abx)) >= get("abx_time") & get("abx") > 1),
      head(.SD, n = 1L),
      by = c(id_vars(res_abx))
    ]
    
    res_abx <- rename_cols(res_abx, "later_abx", index_var(res_abx))
    res_abx <- as_id_tbl(res_abx)
    res_abx <- res_abx[, c(id_vars(res_abx), "later_abx"), with = FALSE]
    
    res_abx <- merge(res[is.na(get("samp_time"))], res_abx)
    res_abx[, c("t_confirm") := pmax(
      get(index_var(res_abx)), get("abx_time"), get("later_abx")
    )]
    res_abx <- res_abx[, c(meta_vars(res_abx), "sep3", "t_confirm"), with=FALSE]
    
    res <- rbind(res_samp, res_abx) 
    
  }
  res <- res[, c(meta_vars(res), "sep3", "t_confirm"), with=FALSE]
  res <- rename_cols(res, "stay_id", id_var(res))
  res
}
