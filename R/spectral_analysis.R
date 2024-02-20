##' Spectral Analysis
##'
##' @param data Dataframe with input data.
##' @param x Column in `data`, i.e. `depth`, `height`, or `age`.
##' @param y Column in `data` with variable of interest.
##' @param method Spectral analysis method to apply: `"MTM"` or `"LOWSPEC"`.
##' @param ... Additional arguments to pass to spectral analysis method.
##' @return A [tibble::tibble()] with spectral outcome.
##' @examples
##' dat <- tibble::tibble(a = 1:10,
##'                       b = 11:20,
##'                       c = stats::rnorm(10),
##'                       d = sample(letters[1:3], 10, TRUE))
##' spectral_analysis(dat, x = a, y = c)
##' @seealso plot_spectrum
##' @export
spectral_analysis <- function(data, x, y, method = "MTM", ...) {
  supported_methods <- c("MTM", "LOWSPEC")
  ## if ("MTM" != method) {
  ##   cli::cli_abort(c("Only method MTM is currently supported.",
  ##                    "i" = "Other methods I may implement:",
  ##                    "*" = "FFT = periodogram",
  ##                    "*" = "Blackman-Tukey",
  ##                    "*" = "MTLS?"))
  ## }
  if (!method %in% supported_methods) {
    cli::cli_abort(c("Method {.q {method}} is not currently supported.",
                     "i" = "Supported methods are {.q {supported_methods}}."))
  }

  if (method == "MTM") {
    out <- data |>
      dplyr::select({{x}}, {{y}}) |>
      astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
      astrochron::mtm(output = 1, genplot = FALSE, verbose = FALSE, ...) |>
      tidyr::pivot_longer(c(.data$AR1_90_power,
                            .data$AR1_95_power,
                            .data$AR1_99_power),
                          names_to = c("AR1", ".width"),
                          names_pattern = "^(AR1)_(9[950])",
                          values_to = "ar1_power") |>
      dplyr::select(-.data$AR1) |>
      dplyr::mutate(.width = readr::parse_double(paste0(".", .data$.width))) |>
      dplyr::rename(frequency = .data$Frequency,
                    power = .data$Power,
                    harmonic_cl = .data$Harmonic_CL,
                    ar1_cl = .data$AR1_CL,
                    ar1_fit = .data$AR1_fit)
  } else if (method == "LOWSPEC") {
   out <- data |>
      dplyr::select({{x}}, {{y}}) |>
      astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
      astrochron::lowspec(output = 1, genplot = FALSE, verbose = FALSE, ...) |>
      tidyr::pivot_longer(c(.data$LOWSPEC_90_power,
                            .data$LOWSPEC_95_power,
                            .data$LOWSPEC_99_power),
                          names_to = c("LOWSPEC", ".width"),
                          names_pattern = "^(LOWSPEC)_(9[950])",
                          values_to = "lowspec_power") |>
      dplyr::select(-.data$LOWSPEC) |>
      dplyr::mutate(.width = readr::parse_double(paste0(".", .data$.width))) |>
      dplyr::rename(frequency = .data$Frequency,
                    power = .data$Prewhite_power,
                    harmonic_cl = .data$Harmonic_CL,
                    lowspec_cl = .data$LOWSPEC_CL,
                    lowspec_fit = .data$LOWSPEC_back)
  }
  return(out)
}

##' Nested Spectral Analysis
##'
##' @inheritParams spectral_analysis
##' @param nest Character vector to nest by.
##' @export
nested_spectral_analysis <- function(data, nest, ..., x, y) {
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
      mtm = purrr::map(.data$data,
                       spectral_analysis,
                       ...,
                       x = {{x}}, y = {{y}})
    ) |>
    dplyr::select(-.data$data) |>
    tidyr::unnest(.data$mtm)
}
