
delayedAssign("job_index",
  optparse::make_option(
    c("-i", "--index"),
    type = "integer",
    default = 1L,
    help = "job index [default = \"%default\"]",
    metavar = "index",
    dest = "ind"
  )
)

parse_args <- function(...) {
  parser <- optparse::OptionParser(option_list = list(...))
  tryCatch(optparse::parse_args(parser), error = function(e) {
    optparse::print_help(parser)
    q("no", status = 1, runLast = FALSE)
  })
}

check_index <- function(opt, ...) {

  opt <- as.integer(Sys.getenv("LSB_JOBINDEX", unset = opt$ind))

  assert_that(is.count(opt))

  if (...length() == 1L) {

    arg_opts <- as.list(..1)

    assert_that(opt <= length(arg_opts))

    setNames(arg_opts[[opt]], names(list(...)))

  } else {

    arg_opts <- expand.grid(..., KEEP.OUT.ATTRS = FALSE,
                            stringsAsFactors = FALSE)

    assert_that(opt <= nrow(arg_opts))

    as.list(arg_opts[opt, ])
  }
}
