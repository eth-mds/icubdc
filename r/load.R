
set_settings <- function(x, ..., lst = NULL) {
  x <- data.table::setattr(x, "settings", c(get_settings(x), lst, list(...)))
  x
}

get_settings <- function(x) attr(x, "settings")

load_data <- function(src, regex = "_imp$", pids = NULL, ...) {

  dat <- data_read(src, regex = regex, pids = pids)
  cnc <- data_vars(dat)

  res <- get_response(src, ..., pids = pids)
  dat <- merge(dat, res, all.y = TRUE)

  dat <- reorder_cols(dat, c(meta_vars(dat), data_var(res), cnc))
  dat <- set_settings(dat, regex = regex, lst = get_settings(res))

  dat
}

get_response <- function(src, resp_fun = "first_hypo", ...) {

  res <- get(resp_fun, mode = "function")(src, ...)
  res <- hypo_clean(res)

  assert_that(is_ts_tbl(res), has_settings(res), "hypo" %in% data_vars(res))

  set_settings(res, resp_fun = resp_fun)
}

first_hypo <- function(src, max_onset = hours(120L), ...) {
  load_hypo(src, "hypo", max_onset = max_onset, ...)[hypo_epi == 1L, ]
}

second_hypo <- function(src, ...) load_hypo(src, "hypo", ...)[hypo_epi == 2L, ]

repeat_hypo <- function(src, ...) load_hypo(src, "hypo", ...)[hypo_epi > 1L, ]

any_hypo <- function(src, ...) load_hypo(src, "hypo", ...)

hypo_clean <- function(x) {
  x <- x[, c("hypo_imp", "hypo_epi", "start_time", "onset_time") := NULL]
  x
}

#' @param min_dur If two consecutive hypo episodes are not separate by more than
#' min_dur, they are merged into a single episode
#' @param max_dur If not terminated explicitly (by a non-hypo measurement),
#' a hypo episode is considered over after onset + max_dur
#' @param right,left Duration for which to carry back the label from a given
#' onset. The `left` value corresponds to the time point with which the label
#' is switched on and the `right` value switches the label off again
#' @param min_icu Time points preceding `min_icu` are removed in order to ensure
#' a certain data availability level.
#' @param max_onset For controls (hypo episodes with no onset), time-points
#' greater/equal to `max_onset` are discarded, while cases are only considered
#' if onset > max_onset.
load_hypo <- function(src, cols = "hypo",
                      min_dur = hours(6L), max_dur = hours(4L),
                      left = hours(24L), right = hours(1L),
                      min_icu = hours(6L), max_onset = hours(72L), ...) {

  res <- data_read(src, cols, regex = NULL, ...)
  res <- hypo_term(res, min_dur, max_dur)
  res <- hypo_augm(res)

  res <- res[, hypo := is_true(onset_time >= -left & onset_time <= -right)]
  res <- res[!is_true(onset_time > -right), ]

  res <- res[stay_time >= min_icu & start_time <= max_onset &
             !is_true(start_time - onset_time > max_onset), ]

  set_settings(res, min_dur = min_dur, max_dur = max_dur, left = left,
               right = right, min_icu = min_icu, max_onset = max_onset)
}

#' @return The returned `ts_tbl` contains a column `hypo` which has been
#' modified by merging adjacent hypo episodes and having not explicitly
#' terminated hypo episodes ended by a negative value.
hypo_term <- function(x, min_dur, max_dur) {

  mrg <- function(x) cumsum(c(TRUE, x[-length(x)]))

  trm <- function(x, mx) {
    ind <- max(which(x > 0L))
    sft <- ind + mx
    if (sft <= length(x)) replace(x, sft, -x[ind]) else x
  }

  idv <- id_vars(x)
  idx <- index_var(x)

  assert_that(!any(c("diff", "hdif", "nhyp", "impu") %in% colnames(x)))

  x <- x[!is.na(hypo), diff := c(diff(get(idx)), Inf), by = c(idv)]
  x <- x[is_true(hypo > 0L), hdif := c(diff(get(idx)), Inf), by = c(idv)]
  x <- x[!is.na(hdif), nhyp := .N, by = c(idv)]

  x <- x[is_true(nhyp > 1L), hypo := mrg(hdif > min_dur | diff < hdif),
             by = c(idv)]
  x <- x[, c("diff", "hdif", "nhyp") := NULL]

  x <- x[, impu := data.table::nafill(hypo, "locf"), by = c(idv)]
  x <- x[, impu := data.table::nafill(impu, fill = 0)]

  mxd <- as.integer(
    ceiling(max_dur / as.double(interval(x), units = units(max_dur)))
  )

  x <- x[impu > 0L, hypo := trm(hypo, mxd), by = c(idv, "impu")]
  x <- x[, impu := NULL]

  x
}

#' @return Constructed from a `hypo` column a returned by `hypo_term()`,
#' several columns are added to the passed `ts_tbl` and returned as such:
#' * `hypo_imp`: using an locf imputation scheme, hypo periods are marked by
#' ascending even numbers and the preceding non-hypo periods by odd integers.
#' * `hypo_epi`: hypo episodes are constructed from the `hypo_imp` column such
#' that a given hypo periods and the stretch leading up to it have assigned the
#' same integer.
#' * `start_time`: Time relative to the start of a given hypo episode.
#' * `onset_time`: Either `NA` if a given hypo episode does not contain a hypo
#' onset or the time relative to hypo onset.
hypo_augm <- function(x) {

  shift <- function(x) data.table::fifelse(x > 0L, x * 2L, x * -2L + 1L, 1L)
  timed <- function(tim, i, j) list(tim - tim[i], tim - tim[which(j)[1L]])

  idv <- id_vars(x)

  x <- x[, hypo_imp := data.table::fifelse(hypo == 0L, NA_integer_, hypo)]
  x <- x[, hypo_imp := data.table::nafill(hypo_imp, "locf"), by = c(idv)]
  x <- x[, hypo_imp := data.table::nafill(hypo_imp, fill = 0L)]
  x <- x[, hypo_imp := shift(hypo_imp)]

  x <- x[, hypo_epi := (hypo_imp + 1L) %/% 2L]

  x <- x[, c("start_time", "onset_time") := timed(
    stay_time, 1L, is_true(hypo > 0L)), by = c(idv, "hypo_epi")
  ]

  x
}

setup_data <- function(src, cfg = config("features"),
                       pids = cohort(src, type = "default"),
                       lwr_thresh = hours(-72L), upr_thresh = hours(336L)) {

  msg_ts("setting up ", src)

  dat <- load_concepts(c("hypo", names(cfg)), src, merge = FALSE,
                       patient_ids = pids)
  dat <- lapply(dat, function(x) {
    if (!is_win_tbl(x)) return(x) else {
      x <- expand(x)
      val_var <- setdiff(names(x), meta_vars(x))
      x[, c(val_var) := as.integer(get(val_var))]
      return(unique(x))
    }
  })
  tim <- vapply(dat, is_ts_tbl, logical(1L))

  msg_ts("merging concepts")

  if (all(tim)) {
    sta <- NULL
  } else {
    sta <- merge_lst(dat[!tim])
  }

  dat <- dat[tim]

  while(length(dat) > 1L) {
    dat[[1L]] <- merge(dat[[1L]], dat[[2L]], all = TRUE)
    dat[[2L]] <- NULL
  }

  dat <- dat[[1L]]

  lim <- collapse(dat)
  lim <- lim[start > 0, start := 0]
  lim <- lim[start < lwr_thresh, start := lwr_thresh]
  lim <- lim[end > upr_thresh, end := upr_thresh]

  msg_ts("making the time series regular")
  dat <- fill_gaps(dat, limits = lim)

  if (not_null(sta)) {
    dat <- merge(dat, sta, all.x = TRUE)
  }

  tim <- tim[names(tim) != "hypo"]

  if (!all(tim)) {
    sta <- xtr_chr(cfg[names(tim)[!tim]], "name")
    dat <- rename_cols(dat, paste0(sta, "_sta"), names(sta), by_ref = TRUE)
  }

  tim <- xtr_chr(cfg[names(tim)[tim]], "name")
  dat <- rename_cols(dat, paste0(tim, "_raw"), names(tim), by_ref = TRUE)

  dat <- carry_forward_to(dat, xtr_num(cfg[names(tim)], "carry_fwd"),
                          paste0(tim, "_raw"))

  imp <- xtr_num(cfg, "impute_val")
  nme <- xtr_chr(cfg, "name")
  imp <- set_names(imp, nme)
  nme <- paste0(nme, ifelse(names(cfg) %in% names(tim), "_imp", "_sta"))

  dat <- impute_values(dat, set_names(imp, nme), nme,
                       sub("_sta$", "_imp", nme))

  dat <- augment(dat, Negate(is.na), "ind", info = "indicators")
  col <- grep("_ind$", colnames(dat), value = TRUE)
  dat <- augment(dat, cumsum, cols = col, names = sub("_ind$", "_cnt", col),
                 by = id_vars(dat), info = "counts")

  for (win in c(4L, 8L, 16L)) {

    for (fun in c("min", "max", "mean", "var")) {

      col <- grep("_raw$", colnames(dat), value = TRUE)
      suf <- paste0("_", fun, win, "lbk")
      new <- sub("_raw$", suf, col)

      dat <- augment(dat, fun, cols = col, names = new, win = win,
        info = paste0(sub("s$", " ", format(hours(win))), fun, "() lookback")
      )
      dat <- impute_values(dat, set_names(imp, paste0(names(imp), suf)), new)

      gc()
    }

    col <- grep("_ind$", colnames(dat), value = TRUE)
    new <- sub("_ind$", paste0("_any", win, "lbk"), col)
    dat <- augment(dat, "any", cols = col, names = new, win = win,
      info = paste0(sub("s$", " ", format(hours(win))), "any() lookback")
    )

    gc()
  }

  data_write(dat, src)

  msg_ts("completed set up of ", src)

  invisible(NULL)
}

indicator_encoding <- function(dat, cfg) {
  
  encode <- function(sequ, is_inc, name, inds, dat, mat) {

    x <- dat[[name]]

    if (is.null(x)) {
      ival <- cut(sequ, sequ, right = !is_inc)
    } else {
      ival <- cut(x, sequ, right = !is_inc)
    }
    
    if (is_inc) {
      col_names <- sub(",\\d+(\\.\\d+)?\\)$", ",Inf)", levels(ival))
    } else {
      col_names <- sub("^\\(\\d+(\\.\\d+)?,", "(-Inf,", levels(ival))
    }

    res <- matrix(FALSE, nrow = nrow(mat), ncol = nlevels(ival),
                  dimnames = list(NULL, col_names))
    
    if (not_null(x)) {
      
      if (is_inc) {

        seq_fun <- function(i) seq.int(1L, i)
        res[!is.na(x) & x >= cfg[[name]][["upper"]], ] <- TRUE
        
      } else {

        seq_fun <- function(i) seq.int(i, ncol(res))
        res[!is.na(x) & x <= cfg[[name]][["lower"]], ] <- TRUE
      }
      
      lvls <- as.integer(ival)
      
      for (i in seq_len(nlevels(ival))) {
        res[!is.na(lvls) & lvls == i, seq_fun(i)] <- TRUE
      }
    }

    mat[, inds] <- res

    data.table::data.table(
      col_inds = inds,
      concept = rep(name, length(inds)),
      threshold = if (is_inc) sequ[-length(sequ)] else sequ[-1],
      right = rep(is_inc, length(inds))
    )
  }

  imp_vars <- paste0(xtr_chr(cfg, "name"), "_imp")

  assert_that(all(colnames(dat) %in% imp_vars))


  seqs <- Map(seq, xtr_num(cfg, "lower"), xtr_num(cfg, "upper"),
              xtr_num(cfg, "step"))
  incr <- xtr_chr(cfg, "direction") == "increasing"

  inds <- cumsum(lengths(seqs) - 1L)
  inds <- Map(seq.int, c(1L, inds[-length(inds)] + 1L), inds)

  mat <- bigmemory::filebacked.big.matrix(
    nrow(dat), sum(lengths(seqs) - 1L), type = "double",
    backingfile = "ind_enc.mat", descriptorfile = "ind_enc.desc",
    backingpath = create_dir(tempfile()), binarydescriptor = TRUE
  )

  res <- Map(encode, seqs, incr, imp_vars, inds,
             MoreArgs = list(dat = dat, mat = mat))

  list(data = mat, meta = data.table::rbindlist(res))
}
