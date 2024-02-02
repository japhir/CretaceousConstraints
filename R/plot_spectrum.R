##' Plot Spectrum
##'
##' Plot the result of [spectral_analysis()] or [nested_spectral_analysis()].
##'
##' @param spec Tibble with output of [spectral_analysis()] or [nested_spectral_analysis()].
##' @param group Character(1). Group column to colour by.
##' @param logx,logy Logical(1). Make the x- and/or y-axis log10.
##' @param confidence Logical(1). Plot the AR1 or LOWSPEC confidence levels.
##' @param periods Periods of interest to plot along the top axis.
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
      periods <- c(405, 132.5, 124, 99.7, 95)
    } else if (domain == "depth") {
      sec_xlab <- ""
    }
  }

  if ("none" %in% group) {
    pl <- spec |>
      ggplot2::ggplot(ggplot2::aes(x = .data$frequency, y = .data$power))
    if (confidence) {
      if ("ar1_fit" %in% colnames(spec)) {
        pl <- pl +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ar1_fit,
                                            ymax = .data$ar1_power,
                                            linetype = NA, group = .data$.width),
                               alpha = .1)
      } else if ("lowspec_fit" %in% colnames(spec)) {
        pl <- pl +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
                                            ymax = .data$lowspec_power,
                                            linetype = NA,
                                            group = .data$.width),
                               alpha = .1)
      }
    }
  } else {
    pl <- spec |>
      ggplot2::ggplot(ggplot2::aes(x = .data$frequency,
                                   y = .data$power,
                                   colour = .data[[group]]))
    if (confidence) {
      if ("ar1_fit" %in% colnames(spec)) {
        pl <- pl +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ar1_fit,
                                            ymax = .data$ar1_power,
                                            fill = .data[[group]],
                                            linetype = NA,
                                            group = paste(c(.data[[group]],
                                                            .data$.width))),
                               alpha = .1)
      } else if ("lowspec_fit" %in% colnames(spec)) {
        pl <- pl +
          # this is a bit silly, cannot pass group = paste(group, .width)
          # so have to repeat myself
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
                                            ymax = .data$lowspec_power,
                                            fill = .data[[group]],
                                            linetype = NA),
                       alpha = .1,
                       data = \(x) x |> filter(.width == .9)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
                                            ymax = .data$lowspec_power,
                                            fill = .data[[group]],
                                            linetype = NA),
                               alpha = .1,
                               data = \(x) x |> filter(.width == .95)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lowspec_fit,
                                            ymax = .data$lowspec_power,
                                            fill = .data[[group]],
                                            linetype = NA),
                               alpha = .1,
                               data = \(x) x |> filter(.width == .99))
      }
    }
  }

  if (logx && logy) {
    pl <- pl + ggplot2::annotation_logticks() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_log10(
        sec.axis = ggplot2::sec_axis(breaks = periods,
                                     trans = \(x) 1 / x,
                                     name = sec_xlab))
  } else if (logx){
    pl <- pl + ggplot2::annotation_logticks(sides = "b") +
      ggplot2::scale_x_log10(
        sec.axis = ggplot2::sec_axis(breaks = periods,
                                     trans = \(x) 1 / x,
                                     name = sec_xlab))
  } else if (logy) {
    pl <- pl +
      ggplot2::annotation_logticks(sides = "l") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::sec_axis(breaks = periods,
                                     trans = \(x) 1 / x,
                                     name = sec_xlab))
  } else {# neither
    pl <- pl +
      ggplot2::scale_x_continuous(
        sec.axis = ggplot2::sec_axis(breaks = periods,
                                     trans = \(x) 1 / x,
                                     name = sec_xlab))
  }

  pl <- pl +
    ggplot2::geom_line() +
    ggplot2::labs(x = xlab,
                  y = "Spectral power (-)")

  return(pl)
}
