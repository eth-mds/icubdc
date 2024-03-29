
#' Misc utilities
#' @param src Data source name
#' @export
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
