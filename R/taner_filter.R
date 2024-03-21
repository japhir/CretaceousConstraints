taner_filter <- function(data, frequencies, x, y, ...,
                         ## linterp = TRUE,
                         ## linterp_dt = NULL,
                         roll = 1e3) {
  data |>
    mutate(filt = list(frequencies)) |>
    unnest("filt") |>
    nest(.by = all_of(colnames(frequencies))) |>
    mutate(lt = purrr::map(data, \(d) d |> dplyr::select({{x}}, {{y}}) |>
                                      astrochron::linterp(genplot = FALSE, verbose = FALSE)),
           bp = purrr::pmap(list(lt, flow, fhigh),
                            \(d, l, h)
                            astrochron::taner(d, flow = l, fhigh = h,
                                              roll = roll,
                                              ...))) |>
    select(-"lt", -"data") |>
    unnest(cols = c("bp"))
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
