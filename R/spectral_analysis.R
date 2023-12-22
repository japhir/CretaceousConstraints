spectral_analysis <- function(data, x, y) {
  data |>
    select({{x}}, {{y}}) |>
    astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
    astrochron::mtm(output = 1, genplot = FALSE, verbose = FALSE) |>
    pivot_longer(c(AR1_90_power, AR1_95_power, AR1_99_power),
                 names_to = c("AR1", ".width"), names_pattern = "^(AR1)_(9[950])",
                 values_to = "ar1_power") |>
    select(-AR1) |>
    mutate(.width = parse_double(paste0(".", .width))) |>
    rename(frequency = Frequency, power = Power, harmonic_cl = Harmonic_CL,
           ar1_cl = AR1_CL, ar1_fit = AR1_fit, .width = .width)
}

nested_spectral_analysis <- function(data, nest, x, y) {
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
    nest(.by = all_of(nest)) |>
    mutate(
      mtm = map(data, spectral_analysis, x = {{x}}, y = {{y}})
    ) |>
    select(-data) |>
    unnest(mtm)
}
