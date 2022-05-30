
json_file <- function(name, dir, value = NULL, simplifyVector = TRUE,
  simplifyDataFrame = FALSE, simplifyMatrix = FALSE, ...) {
  
  assert_that(dir.exists(dir))
  
  file <- paste0(file.path(dir, name), ".json")
  
  if (!is.null(value)) {
    
    assert_that(is.list(value))
    jsonlite::write_json(value, file, ...)
    
  } else {
    
    if (!file.exists(file)) {
      stop("config file ", basename(file), " does not exists.")
    }
    
    jsonlite::read_json(file, simplifyVector = simplifyVector,
      simplifyDataFrame = simplifyDataFrame,
      simplifyMatrix = simplifyMatrix, ...)
  }
}

config <- function(name, value = NULL, ...) {
  json_file(name, here::here("config"), value, ...)
}

cohort <- function(src, type = "default") {

  # `default` and `all` refer to the same set of patients:
  # - `default` is used for `setup_data()`
  # - `all` is used in `train_predict()` and just uses everything

  if (identical(type, "all")) {
    return(NULL)
  }

  if (identical(type, "ins_avail") && !identical(src, "eicu")) {
    type <- "default"
  }

  res <- config("cohort")

  assert_that(has_name(res, src), has_name(res[[src]], type))

  res[[src]][[type]]
}

ensure_dir <- function(dir) {

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  dir
}

results_dir <- function(subdir = NULL) {

  dir <- here::here("results")

  if (not_null(subdir)) {
    dir <- file.path(dir, subdir)
  }

  ensure_dir(dir)
}

data_path <- function(src, v_ricu = packageVersion("ricu"),
                      v_hypo = data_version) {

  file <- paste0(src, "-", sub("\\.", "_", v_ricu), "-", v_hypo, ".fst")
  file <- file.path(ensure_dir(here::here("data")), file)

  set_settings(file, v_ricu = v_ricu, v_hypo = v_hypo)
}

data_read <- function(src, cols = NULL, regex = NULL, pids = NULL,
                      time_filter = NULL, ..., ts_tbl = TRUE) {

  fil <- data_path(src, ...)

  if (!file.exists(fil)) {
    stop("could not read `", src, "` data from ", fil)
  }

  if (is_lsf()) {
    file.copy(fil, tempdir())
    fil <- file.path(tempdir(), basename(fil))
  }

  fst <- fst::fst(fil)

  if (is.null(cols)) {
    cols <- colnames(fst)
  }

  if (not_null(regex)) {
    cols <- grep(regex, cols, value = TRUE)
  }

  cols <- unique(c("icu_stay", "stay_time", cols))

  if (is.null(pids) && is.null(time_filter)) {

    dat <- fst::read_fst(fil, cols, as.data.table = TRUE)

  } else {

    rows <- TRUE

    if (not_null(pids)) {
      rows <- rows & (fst[["icu_stay"]] %in% pids)
    }

    if (not_null(time_filter)) {
      rows <- rows & time_filter(fst[["stay_time"]])
    }

    dat <- fst[rows, cols]
  }

  if (isTRUE(ts_tbl)) {
    dat <- as_ts_tbl(dat, "icu_stay", "stay_time", interval = hours(1L),
                     by_ref = TRUE)
  }

  set_settings(dat, lst = get_settings(fil))
}

data_write <- function(x, src, ...) {

  assert_that(is_ts_tbl(x), all.equal(interval(x), hours(1L)))

  fil <- data_path(src, ...)

  if (is_lsf()) {
    dst <- dirname(fil)
    fil <- file.path(tempdir(), basename(fil))
  }

  msg_ts("writing data to ", fil)

  x <- rename_cols(x, c("icu_stay", "stay_time"), meta_vars(x), by_ref = TRUE)

  fst::write_fst(x, fil, compress = 100L)

  if (is_lsf()) {
    msg_ts("moving temp file to ", dst)
    file.copy(fil, dst)
  }

  invisible(NULL)
}
