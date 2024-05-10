#' Taner filter wrapper
#' @inheritParams bandpass_filter
#' @param roll Roll-off rate, in dB/octave.  Typical values are \eqn{10^3} to \eqn{10^12}, but can be larger.
#' @param ... Additional parameters, passed on to [astrochron::taner()].
#' @seealso [astrochron::taner()] for the original implementation, [astrochron::linterp] for linearly interpolating.
#' @returns A [tibble][tibble::tibble-package] with the target frequencies and the associated filter output.
#' @examples
#' dat <- tibble::tibble(a = 1:10,
#'                       b = 11:20,
#'                       c = stats::rnorm(10),
#'                       d = sample(letters[1:3], 10, TRUE), depth = 11:20)
#' # TODO: fix example so it actually filters something
#' taner_filter(dat,
#'              tibble::tibble(flow = 1, fhigh = 2, target = "group"),
#'              x = a, y = c,
#'              genplot = FALSE, verbose = FALSE, add_depth = TRUE)
#' @export
taner_filter <- function(data, frequencies, x = NULL, y = NULL, ...,
                         ## linterp = TRUE,
                         linterp_dt = NULL,
                         roll = 1e3,
                         add_depth = FALSE) {
  null_tracker <- 0L
  if (rlang::quo_is_null(enquo(x))) {
    null_tracker <- null_tracker + 1
  }
  if (rlang::quo_is_null(enquo(y))) {
    null_tracker <- null_tracker + 1
  }
  if (null_tracker == 1L) {
    stop("Must provide either both 'x' and 'y', or neither.")
  }

  # linterp by default!
  out <- data |>
    mutate(filt = list(frequencies)) |>
    unnest("filt") |>
    nest(.by = tidyr::all_of(colnames(frequencies))) |>
    mutate(lt = purrr::map(data, \(d) {
      if (null_tracker == 2L) {
        d <- d[, 1:2]
      } else {
        d <- d |>
          dplyr::select({{x}}, {{y}})
      }
      d |>
        astrochron::linterp(genplot = FALSE,
                            verbose = FALSE,
                            dt = linterp_dt)
    })) |>
    mutate(bp = purrr::pmap(list(lt, flow, fhigh),
                            \(d, l, h) {
                              d <- astrochron::taner(d, flow = l, fhigh = h,
                                                     roll = roll,
                                                     ...)
                              if (null_tracker == 2L) {
                                d <- dplyr::rename(d, filter = 2)
                              } else {
                                d <- dplyr::rename(d, filter = {{y}})
                              }
                              # remove depth so unnesting works. make sure it's a named tibble still
                              d |> dplyr::select("filter")
                            })
           ) |>
    unnest(cols = c("lt", "bp"))


  # interpolate depth back to new age scale
  if (add_depth) {
    if (null_tracker > 0L) {
      message("Cannot add depth when 'x' and/or 'y' are NULL.")
    } else {
      out <- out |>
        select(-"data") |>
        dplyr::mutate(depth = stats::approx(x = dplyr::pull(data, {{x}}),
                                            y = data$depth,
                                            xout = dplyr::pull(out, {{x}}))$y,
                      .before = {{x}})
    }
  }

  out |>
    select(-tidyselect::any_of("data"))
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
                                 genplot = FALSE, verbose = FALSE, ...) ## |>
                      ## select(-nest) # otherwise it'll get duplicated
               )) |>
    tidyr::unnest("bp") |>
    dplyr::select(-"data")
}
