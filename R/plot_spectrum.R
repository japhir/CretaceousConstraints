##' Plot Spectrum
##'
##' Plot the result of [spectral_analysis()] or [nested_spectral_analysis()].
##'
##' @param spec A [tibble::tibble()] with the output of [spectral_analysis()]
##'   or [nested_spectral_analysis()].
##' @param group Character(1) Column name to colour lines and ribbons by.
##' @param logx,logy If `TRUE`, plot the axis on a log10 scale. Defaults to
##'   `FALSE`.
##' @param domain Domain of spectral analysis: `"depth"` or `"time"`.
##' @param confidence If `TRUE`, plot the AR1 or LOWSPEC confidence levels.
##' @param periods Periods of interest to plot along the top axis. Defaults to
##'   `NULL` if we're in the depth domain, and to the short and long
##'   eccentricity cycles if we're in the time domain.
##' @return A [ggplot2::ggplot()] object.
##' @examples
##' dat <- tibble::tibble(a = 1:10,
##'                       b = 11:20,
##'                       c = stats::rnorm(10),
##'                       d = sample(letters[1:3], 10, TRUE))
##' plot_spectrum(spectral_analysis(dat, x = a, y = c))
##' @export
plot_spectrum <- function(spec,
                          group = "none",
                          logx = FALSE,
                          logy = FALSE,
                          domain = "depth",
                          confidence = FALSE,
                          periods = NULL) {
  if (domain == "depth") {
    xlab <- "Frequency (cycles" ~ m^{-1} * ")"
    sec_xlab <- "Period (m)"
  } else if (domain == "time") {
    xlab <- "Frequency (cycles" ~ kyr^{-1} * ")"
    sec_xlab <- "Period (kyr)"
  }

  if (is.null(periods)) {
    if (domain == "time") {
      periods <- c(405, 131, 124, 99, 95)
    } else if (domain == "depth") {
      sec_xlab <- ""
    }
  }

  if (confidence) {
      if (!"background_fit" %in% colnames(spec)) {
        cli::cli_warn("{.q background_fit} not found in {colnames(spec)}.")
      }
  }

  if ("none" %in% group) {
    pl <- spec |>
      ggplot2::ggplot(ggplot2::aes(x = .data$frequency, y = .data$power))
    if (confidence) {
      pl <- pl +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$background_fit,
                                          ymax = .data$background_power,
                                          group = .data$.width),
                             outline.type = "lower",
                             alpha = .1)
    }
  } else {
    pl <- spec |>
      ggplot2::ggplot(ggplot2::aes(x = .data$frequency,
                                   y = .data$power,
                                   colour = .data[[group]]))
    if (confidence) {
      if ((spec[[group]] |> unique() |> length()) == 1) {
        pl <- pl +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$background_fit,
                                            ymax = .data$background_power,
                                            fill = .data[[group]],
                                            group = .data$.width,
                                            ),
                               outline.type = "lower",
                               alpha = .1)
      } else {
        pl <- pl +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$background_fit,
                                            ymax = .data$background_power,
                                            fill = .data[[group]],
                                            group = paste(.data$.width, .data[[group]]),
                                            ),
                               outline.type = "lower",
                               alpha = .1)
      }
    }
  }
    ##   } else if ("lowspec_fit" %in% colnames(spec)) {
    ##     pl <- pl +
    ##       # this is a bit silly, cannot pass group = paste(group, .width)
    ##       # so have to repeat myself
    ##       ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
    ##                                         ymax = .data$lowspec_power,
    ##                                         fill = .data[[group]],
    ##                                         linetype = NA),
    ##                    alpha = .1,
    ##                    data = \(x) x |> filter(.width == .9)) +
    ##       ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
    ##                                         ymax = .data$lowspec_power,
    ##                                         fill = .data[[group]],
    ##                                         linetype = NA),
    ##                            alpha = .1,
    ##                            data = \(x) x |> filter(.width == .95)) +
    ##       ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
    ##                                         ymax = .data$lowspec_power,
    ##                                         fill = .data[[group]],
    ##                                         linetype = NA),
    ##                            alpha = .1,
    ##                            data = \(x) x |> filter(.width == .99))
    ##   }
    ## }

  if (!logx){
    pl <- pl +
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::sec_axis(breaks = periods,
                                     trans = \(x) 1 / x,
                                     name = sec_xlab))
  } else {
    pl <- pl +
      ggplot2::scale_x_continuous(transform = "log10",
                                  guide = guide_axis_logticks(long = 2, mid = 1, short = 0.5),
                                  sec.axis = ggplot2::sec_axis(breaks = periods,
                                                               trans = \(x) 1 / x,
                                                               name = sec_xlab))
  }

  if (logy) {
    pl <- pl +
      ggplot2::scale_y_continuous(transform = "log10",
                                  guide = guide_axis_logticks(long = 2, mid = 1, short = 0.5))
  }

  pl <- pl +
    ggplot2::geom_line() +
    ggplot2::labs(x = xlab,
                  y = "Spectral power (-)")

  return(pl)
}
