#' Bandpass filter
#'
##' @param data data.frame()
##' @param frequencies data.frame() with columns `target`, `flow`, and `fhigh`.
##'   Can contain multiple rows for multiple frequency filtering.
##' @param x Column name in `data` that holds the depth/age information.
##' @param y Column name in `data` that holds the proxy variable name.
##' @param add_depth Logical(1) Add depth back in if x is age?
##' @export
bandpass_filter <- function(data, frequencies, x, y, add_depth = FALSE) {
  if (! "data.frame" %in% class(data)) {
    cli::cli_abort(c(
           "{.var data} must be a {.cls data.frame}",
           "x" = "You supplied a {.cls {class(data)}}"))
  }
  if (! "data.frame" %in% class(frequencies)) {
    cli::cli_abort(c(
           "{.var freqs} must be a {.cls data.frame}",
           "x" = "You've supplied a {.cls {class(freqs)}}"))
  }
  if (! all(c("fhigh", "flow", "target") %in% colnames(frequencies))) {
    cli::cli_abort(c("{.var freqs} must have columns `flow` and `fhigh`",
                     "i" = "{.var freqs} has column{?s} {.val {colnames(freqs)}}"))
  }

  out <- data |>
    mutate(filt = list(frequencies)) |>
    unnest(filt) |>
    nest(.by = all_of(colnames(frequencies)))

  out <- out |>
    mutate(
      lt = map(data, \(d) d |>
                          select({{x}}, {{y}}) |>
                          astrochron::linterp(genplot = FALSE, verbose = FALSE)),
      bp = pmap(list(lt, flow, fhigh), \(d, l, h)
                d |>
                astrochron::bandpass(flow = l, fhigh = h, win = 0,
                                     genplot = FALSE, verbose = FALSE) |>
                select(filter = {{y}}))) |>
    select(-data) |>
    unnest(cols = c(lt, bp))

  if (add_depth) {
    # interpolate depth back to new age scale
    out <- out |>
      mutate(depth = approx(data |> pull({{x}}),
                            data$depth,
                            out |> pull({{x}}))$y,
             .before = {{x}})
  }

  out
}


#' Nested bandpass filter
#'
#' @inheritParams bandpass_filter
#' @param nest character() Vector with column names in `data` to nest by.
nested_bandpass_filter <- function(data, frequencies, x, y, nest, add_depth = FALSE) {
  if (! all(nest %in% colnames(data))) {
    cli::cli_abort(c("{.var nest} must have columns that exist in {.var data}",
                     "i" = "{.var data} has column{?s} {.val {colnames(data)}}",
                     "x" = "You've supplied {.val {nest}}"))
  }

  data |>
    nest(.by = all_of(nest)) |>
    mutate(bp = purrr::map(data,
               \(d) d |>
                    bandpass_filter(frequencies = frequencies,
                                    x = {{x}}, y = {{y}},
                                    add_depth = add_depth))) |>
    unnest(bp) |>
    select(-data)
}
