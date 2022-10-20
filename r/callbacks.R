
#' Callback functions for project-specific concept items
#' @param beta Vasopressor adjustement
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

#' @param ins Object containing insulin values
#' @param ... Further arguments are ignored
#' @rdname map_beta
#' @export
ins_ifx_cb <- function(ins, ...) {
  if (identical(id_vars(ins), "icustay_id")) {
    ins[get("ins") == 0, c("ins") := 2]
  } else if ("source" %in% id_vars(ins)) {
    ins[get("ins") == 0 & grepl("mimic", get("source")), c("ins") := 2]
  }
  rename_cols(ins, "ins_ifx", "ins")
}

#' @param glu Object containing glucose measurements
#' @rdname map_beta
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

#' @param hypo Object containing a hypoglycemia indicator
#' @param min_dur Minimal episode duration interval
#' @rdname map_beta
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

  hypo <- hypo[!is.na(get("hypo")), c("diff") := c(diff(get(idx)), Inf),
               by = c(idv)]
  hypo <- hypo[is_true(get("hypo") > 0L), c("hdif") := c(diff(get(idx)), Inf),
               by = c(idv)]
  hypo <- hypo[!is.na(get("hdif")), c("nhyp") := .N, by = c(idv)]

  hypo <- hypo[is_true(get("nhyp") > 1L),
    c("hypo") := mrg(get("hdif") > min_dur | get("diff") < get("hdif")),
    by = c(idv)
  ]
  hypo <- hypo[, c("diff", "hdif", "nhyp") := NULL]
  hypo <- hypo[, c("hypo_epi") := data.table::nafill(get("hypo"), fill = 0)]

  hypo
}

#' @param hypo_epi Object containing a hypoglycemia episode index
#' @rdname map_beta
#' @export
hypo_cnt_cb <- function(hypo_epi, ...) {

  hypo_epi <- fill_gaps(hypo_epi)
  hypo_epi <- hypo_epi[, lapply(.SD, nafill, type = "locf"),
                       by = c(id_vars(hypo_epi))]
  hypo_epi <- replace_na(hypo_epi, 0)
  hypo_epi[, c("hypo_cnt") := cummax(get("hypo_epi")),
           by = c(id_vars(hypo_epi))]

  hypo_epi[, c(meta_vars(hypo_epi), "hypo_cnt"), with=FALSE]
}

#' @param dex Object containing dextrose values
#' @rdname map_beta
#' @export
dex_amount_callback <- function(dex, interval, ...) {
  c_fct <- as.double(interval, units = units(dex$dur_var))
  # make into amount
  dex[, c("dex") := get("dex") * as.double(get("dur_var")) / c_fct]
  dex[, c("distr") := as.double(ricu:::re_time(get("dur_var"), interval)) + 1]
  dex[, c("dex") := get("dex") / get("distr")]
  dex[, c("distr") := NULL]
  dex <- rename_cols(dex, "dex_amount", "dex")
  
  expand(dex, aggregate = "sum")
}

#' @param interval,dur Duration
#' @rdname map_beta
#' @export
ts_to_win_hours <- function(dur) {

  assert_that(ricu:::is_interval(dur))

  function(x, dur_var, ...) {
    x[, c(list(...)$val_var) := NULL]
    x[, c(list(...)$val_var) := TRUE]
    x[, c(dur_var) := NULL]
    x[, c(dur_var) := dur]
    as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
  }
}

#' @rdname map_beta
#' @export
hirid_pharma_win <- function(dur) {

  assert_that(ricu:::is_interval(dur))

  function(x, dur_var, group_var, ...) {

    x[, c(list(...)$val_var) := NULL]
    x[, c(list(...)$val_var) := TRUE]
    x[, c(dur_var) := NULL]

    x[, c(dur_var) := max(get(index_var(x))) - min(get(index_var(x))) + dur,
      by = c(group_var)]
    x <- x[, head(.SD, n = 1L), by = c(group_var)]
    x[, c(dur_var) := `units<-`(get(dur_var), "mins")]
    as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
  }
}

#' @param x Data
#' @param val_var,dur_var Value and duration columns
#' @rdname map_beta
#' @export
mimic_presc_cort <- function(x, val_var, dur_var, ...) {
  x[, c(val_var) := NULL]
  x[, c(val_var) := TRUE]
  x[, c(index_var(x)) := get(index_var(x)) + hours(9L)]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(360L)]
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

#' @rdname map_beta
#' @export
aumc_cortico <- function(x, dur_var, ...) {
  x[, c(dur_var) := get(dur_var) + mins(360L)]
  x
}
