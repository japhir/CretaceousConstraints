taner_filter <- function(data, frequencies, x, y, ...,
                         ## linterp = TRUE,
                         linterp_dt = NULL,
                         roll = 1e3,
                         add_depth = FALSE) {
  # linterp by default!
  out <- data |>
    mutate(filt = list(frequencies)) |>
    unnest("filt") |>
    nest(.by = all_of(colnames(frequencies))) |>
    mutate(lt = purrr::map(data, \(d) d |> dplyr::select({{x}}, {{y}}) |>
                                        astrochron::linterp(genplot = FALSE,
                                                            verbose = FALSE,
                                                            dt = linterp_dt)),
           bp = purrr::pmap(list(lt, flow, fhigh),
                            \(d, l, h)
                            d |>
                              astrochron::taner(flow = l, fhigh = h,
                                                roll = roll,
                                                ...) |>
                              dplyr::select(filter = {{y}}))) |>
    select(-"data") |>
    unnest(cols = c("lt", "bp"))

  # interpolate depth back to new age scale
  if (add_depth) {
    out <- out |>
      dplyr::mutate(depth = stats::approx(data |> dplyr::pull({{x}}),
                                          data$depth,
                                          out |> dplyr::pull({{x}}))$y,
                    .before = {{x}})
  }

  out
}

nested_taner_filter <- function(data, frequencies, x, y, nest, ...) {
    if (! all(nest %in% colnames(data))) {
    cli::cli_abort(c("{.var nest} must have columns that exist in {.var data}",
                     "i" = "{.var data} has column{?s} {.val {colnames(data)}}",
                     "x" = "You've supplied {.val {nest}}"))
  }

  data |>
    tidyr::nest(.by = tidyr::all_of(nest)) |>
    dplyr::mutate(bp = purrr::map(data,
               \(d) d |>
                    taner_filter(frequencies = frequencies,
                                    x = {{x}}, y = {{y}},
                                    ...))) |>
    tidyr::unnest("bp") |>
    dplyr::select(-"data")
}
