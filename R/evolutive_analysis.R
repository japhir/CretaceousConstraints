##' Evolutive Harmonics Analysis Wrapper
##'
##' @param data Tibble TODO.
##' @param nest Tibble TODO.
##' @param x,y Column names.
##' @export
evolutive_analysis <- function(data, nest, x, y) {
  if (! "data.frame" %in% class(data)) {
    cli::cli_abort(c(
           "{.var data} must be a {.cls data.frame}",
           "x" = "You supplied a {.cls {class(data)}}"))
  }
  if (! all(nest %in% colnames(data))) {
    cli::cli_abort(c("{.var nest} must have columns that exist in {.var data}",
                     "i" = "{.var data} has column{?s} {.val {colnames(data)}}",
                     "x" = "You've supplied {.val {nest}}"))
  }

  data |>
    tidyr::nest(.by = tidyr::all_of(nest)) |>
    dplyr::mutate(
      eha = purrr::map(.data$data,
                ~ . |>
                  dplyr::select({{x}}, {{y}}) |>
                  astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
                  astrochron::eha(output = 3, genplot = FALSE, verbose = FALSE) |>
                  tidyr::pivot_longer(-freq) |>
                  dplyr::mutate(name = stringr::str_sub(name, 2, -1) |>
                                  readr::parse_double()) |>
                  dplyr::rename({{x}} := name))
    ) |>
    dplyr::select(-.data$data) |>
    tidyr::unnest(.data$eha)
}
