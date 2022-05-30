
Sys.setenv(
  RICU_CONFIG_PATH = here::here("config"),
  RICU_SRC_LOAD = "mimic,mimic_demo,eicu,eicu_demo,hirid,aumc,miiv"
)

ncpu <- n_cores()

options(fst_threads = ncpu)

pkgs <- c("ricu", "jsonlite", "ranger", "data.table", "biglasso", "here",
          "qs", "optparse", "roll", "memuse", "precrec")

if (!all(vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE))) {
  stop("Packages ", paste0(pkgs, collapse = ", "),
       " are required in order to proceed.")
  if (!interactive()) q("no", status = 1, runLast = FALSE)
}

data.table::setDTthreads(ncpu)

library(ricu)

data_version <- "003"

if (is_lsf()) {
  RcppParallel::setThreadOptions(numThreads = ncpu)
}
