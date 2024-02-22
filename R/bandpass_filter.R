#' Bandpass filter
#'
##' @param data data.frame() to bandpass filter.
##' @param frequencies data.frame() with columns `target`, `flow`, and `fhigh`.
##'   Can contain multiple rows for multiple frequency filtering.
##' @param x Column name in `data` that holds the depth/age information.
##' @param y Column name in `data` that holds the proxy variable name.
##' @param ... Additional parameters, yet to be determined.
##' @param linterp If `TRUE`, linearly interpolate.
##' @param linterp_dt Resolution to linearly interpolate to.
##' @param window Window type for bandpass filter: 0 = rectangular, 1 = Gaussian, 2 = Cosine-tapered window (a.k.a. Tukey window).
##' @param add_depth Defaults to `FALSE`. If `TRUE`, find column `depth` in original and interpolate it back from age.
##' @examples
##' dat <- tibble::tibble(a = 1:10,
##'                       b = 11:20,
##'                       c = stats::rnorm(10),
##'                       d = sample(letters[1:3], 10, TRUE))
##' bandpass_filter(dat,
##'                 tibble::tibble(flow = 1, fhigh = 2, target = "group"),
##'                 x = a, y = c)
##' @export
##' @seealso [astrochron::bandpass], [astrochron::linterp]
bandpass_filter <- function(data, frequencies,
                            x, y, ...,
                            linterp = TRUE,
                            linterp_dt = NULL,
                            window = 0,
                            add_depth = FALSE) {
  if (! "data.frame" %in% class(data)) {
    cli::cli_abort(c(
           "{.var data} must be a {.cls data.frame}",
           "x" = "You supplied a {.cls {class(data)}}"))
  }
  if (! "data.frame" %in% class(frequencies)) {
    cli::cli_abort(c(
           "{.var frequencies} must be a {.cls data.frame}",
           "x" = "You've supplied a {.cls {class(frequencies)}}"))
  }
  if (! all(c("fhigh", "flow", "target") %in% colnames(frequencies))) {
    cli::cli_abort(c("{.var frequencies} must have columns `flow`, `fhigh`, and `target`.",
                     "i" = "{.var frequencies} has column{?s} {.val {colnames(frequencies)}}"))
  }

  out <- data |>
    dplyr::mutate(filt = list(frequencies)) |>
    tidyr::unnest("filt") |>
    tidyr::nest(.by = tidyr::all_of(colnames(frequencies)))

  if (linterp) {
    out <- out |>
      dplyr::mutate(
               lt = purrr::map(.data$data,
                               \(d) d |>
                                    dplyr::select({{x}}, {{y}}) |>
                                    astrochron::linterp(genplot = FALSE,
                                                        verbose = FALSE,
                                                        dt = linterp_dt)),
               bp = purrr::pmap(list(.data$lt, .data$flow, .data$fhigh),
                                \(d, l, h)
                                d |>
                                  astrochron::bandpass(flow = l, fhigh = h,
                                                       ...,
                                                       win = window,
                                                       demean = TRUE,
                                                       detrend = TRUE,
                                                       genplot = FALSE, verbose = FALSE) |>
                                dplyr::select(filter = {{y}}))) |>
      dplyr::select(-"data") |>
      tidyr::unnest(cols = c("lt", "bp"))
  } else {
    out <- out |>
      dplyr::mutate(
               bp = purrr::pmap(list(.data$data, .data$flow, .data$fhigh),
                                \(d, l, h)
                                d |>
                                astrochron::bandpass(flow = l, fhigh = h,
                                                     ...,
                                                     win = window,
                                                     demean = TRUE, detrend = TRUE,
                                                     genplot = FALSE, verbose = FALSE) |>
                                dplyr::select(filter = {{y}}))) |>
      dplyr::select(-"data") |>
      tidyr::unnest(cols = c("data", "bp"))
  }

  if (add_depth) {
    # interpolate depth back to new age scale
    out <- out |>
      dplyr::mutate(depth = stats::approx(data |> dplyr::pull({{x}}),
                                   data$depth,
                                   out |> dplyr::pull({{x}}))$y,
             .before = {{x}})
  }

  out
}


#' Nested bandpass filter
#'
#' @inheritParams bandpass_filter
#' @param nest character() Vector with column names in `data` to nest by.
nested_bandpass_filter <- function(data,
                                   frequencies,
                                   x, y,
                                   nest,
                                   ...,
                                   window = 0,
                                   add_depth = FALSE) {
  if (! all(nest %in% colnames(data))) {
    cli::cli_abort(c("{.var nest} must have columns that exist in {.var data}",
                     "i" = "{.var data} has column{?s} {.val {colnames(data)}}",
                     "x" = "You've supplied {.val {nest}}"))
  }

  data |>
    tidyr::nest(.by = tidyr::all_of(nest)) |>
    dplyr::mutate(bp = purrr::map(data,
               \(d) d |>
                    bandpass_filter(frequencies = frequencies,
                                    x = {{x}}, y = {{y}},
                                    ...,
                                    window = window,
                                    add_depth = add_depth))) |>
    tidyr::unnest("bp") |>
    dplyr::select(-"data")
}
