#' @importFrom assertthat assert_that is.scalar is.count is.flag on_failure<-
is_intish <- function(x) {
  is.integer(x) || (is.numeric(x) && all(x == trunc(x)) && !is.na(x))
}

on_failure(is_intish) <- function(call, env) {
  paste0("Not all elements of ", deparse(call$x), " are integer valued.")
}

has_settings <- function(x) not_null(get_settings(x))

on_failure(has_settings) <- function(call, env) {
  paste0(deparse(call$x), " does not contain a settings attribute.")
}
