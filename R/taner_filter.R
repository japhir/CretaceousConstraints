taner_filter <- function(data, frequencies, x = NULL, y = NULL, ...,
                         ## linterp = TRUE,
                         linterp_dt = NULL,
                         roll = 1e3,
                         add_depth = FALSE) {
  # linterp by default!
  out <- data |>
    mutate(filt = list(frequencies)) |>
    unnest("filt") |>
    nest(.by = all_of(colnames(frequencies))) |>
    mutate(lt = purrr::map(data, \(d) {
      if (is.null(x) & is.null(y)) {
        d <- d[, 1:2]
      } else {
        d <- d |>
          dplyr::select({{x}}, {{y}})
      }
      d |>
        astrochron::linterp(genplot = FALSE,
                            verbose = FALSE,
                            dt = linterp_dt)}),
      bp = purrr::pmap(list(lt, flow, fhigh),
                       \(d, l, h) {

                       d <- d |>
                         astrochron::taner(flow = l, fhigh = h,
                                           roll = roll,
                                           ...)

                       if (is.null(x) && is.null(y)) {
                         d <- rename(d, filter = 2)
                       } else {
                         d <- d |> dplyr::select(filter = {{y}})
                       }
                       d |> pull(2)
                       })
      ) |>
  unnest(cols = c("lt", "bp"))


  # interpolate depth back to new age scale
  if (add_depth) {
    out <- out |>
      dplyr::mutate(depth = stats::approx(data |> dplyr::pull({{x}}),
                                          data$depth,
                                          out |> dplyr::pull({{x}}))$y,
                    .before = {{x}})
  }

  out |>
    select(-"data")
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
                                 ...) ## |>
                      ## select(-nest) # otherwise it'll get duplicated
               )) |>
    tidyr::unnest("bp") |>
    dplyr::select(-"data")
}
