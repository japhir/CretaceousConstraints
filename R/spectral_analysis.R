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
spectral_analysis <- function(data, x = NULL, y = NULL, method = "MTM", ...) {
  supported_methods <- c("MTM", "MTM AR1", "LOWSPEC", "MTM PL", "MTM ML96")
  ## if ("MTM" != method) {
  ##   cli::cli_abort(c("Only method MTM is currently supported.",
  ##                    "i" = "Other methods I may implement:",
  ##                    "*" = "FFT = periodogram",
  ##                    "*" = "Blackman-Tukey",
  ##                    "*" = "MTLS?"))
  ## }
  # TODO: come up with uniform names for elements of the output so that I can always use the same simple plotting function?
  if (!method %in% supported_methods) {
    cli::cli_abort(c("Method {.q {method}} is not currently supported.",
                     "i" = "Supported methods are {.q {supported_methods}}."))
  }

  if (is.null(x) && is.null(y)) {
    data <- data[, 1:2]
  }

  # TODO switch from if else if else if else to switch?
  if (method == "MTM" || method == "MTM AR1") {
    out <- data |>
      astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
      astrochron::mtm(output = 1, genplot = FALSE, verbose = FALSE, ...) |>
      tidyr::pivot_longer(c(.data$AR1_90_power,
                            .data$AR1_95_power,
                            .data$AR1_99_power),
                          names_to = c("AR1", ".width"),
                          names_pattern = "^(AR1)_(9[950])",
                          values_to = "background_power") |>
      dplyr::select(-"AR1") |>
      dplyr::mutate(.width = readr::parse_double(paste0(".", .data$.width))) |>
      dplyr::rename(frequency = .data$Frequency,
                    power = .data$Power,
                    harmonic_cl = .data$Harmonic_CL,
                    background_cl = .data$AR1_CL,
                    background_fit = .data$AR1_fit)
  } else if (method == "LOWSPEC") {
   out <- data |>
      astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
      astrochron::lowspec(output = 1, genplot = FALSE, verbose = FALSE, ...) |>
      tidyr::pivot_longer(c(.data$LOWSPEC_90_power,
                            .data$LOWSPEC_95_power,
                            .data$LOWSPEC_99_power),
                          names_to = c("LOWSPEC", ".width"),
                          names_pattern = "^(LOWSPEC)_(9[950])",
                          values_to = "background_power") |>
      dplyr::select(-"LOWSPEC") |>
      dplyr::mutate(.width = readr::parse_double(paste0(".", .data$.width))) |>
      dplyr::rename(frequency = .data$Frequency,
                    power = .data$Prewhite_power,
                    harmonic_cl = .data$Harmonic_CL,
                    background_cl = .data$LOWSPEC_CL,
                    background_fit = .data$LOWSPEC_back)
  } else if (method == "MTM PL") {
    out <- data |>
      astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
      astrochron::mtmPL(output = 1, genplot = FALSE, verbose = FALSE, ...) |>
      tidyr::pivot_longer(c(.data$PowerLaw_90_power,
                            .data$PowerLaw_95_power,
                            .data$PowerLaw_99_power),
                          names_to = c("PL", ".width"),
                          names_pattern = "^(PowerLaw)_(9[950])",
                          values_to = "background_power") |>
      dplyr::select(-"PL") |>
      dplyr::mutate(.width = readr::parse_double(paste0(".", .data$.width))) |>
      dplyr::rename(frequency = .data$Frequency,
                    power = .data$Power,
                    harmonic_cl = .data$Harmonic_CL,
                    background_cl = .data$PowerLaw_CL,
                    background_fit = .data$PowerLaw_fit)
  } else if (method == "MTM ML96") {
    out <- data |>
      astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
      astrochron::mtmML96(output = 1, genplot = FALSE, verbose = FALSE, ...) |>
      tidyr::pivot_longer(c(.data$AR1_90_power,
                            .data$AR1_95_power,
                            .data$AR1_99_power),
                          names_to = c("AR1", ".width"),
                          names_pattern = "^(AR1)_(9[950])",
                          values_to = "background_power") |>
      dplyr::select(-"AR1") |>
      dplyr::mutate(.width = readr::parse_double(paste0(".", .data$.width))) |>
      dplyr::rename(frequency = .data$Frequency,
                    power = .data$Power,
                    harmonic_cl = .data$Harmonic_CL,
                    background_cl = .data$AR1_CL,
                    background_fit = .data$AR1_fit)
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
    dplyr::select(-"data") |>
    tidyr::unnest("mtm")
}
