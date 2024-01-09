plot_spectrum <- function(spec, group = "none", logx = FALSE, logy = FALSE, ar1 = FALSE) {

  if ("none" %in% group) {
    pl <- spec |>
      ggplot(aes(x = frequency, y = power))
    if (ar1) {
      pl <- pl +
        geom_ribbon(aes(ymin = ar1_fit, ymax = ar1_power,
                        linetype = NA, group = .width),
                    alpha = .1)
    }
  } else {
    pl <- spec |>
      ggplot(aes(x = frequency, y = power, colour = .data[[group]]))
    if (ar1) {
      pl <- pl +
        geom_ribbon(aes(ymin = ar1_fit, ymax = ar1_power, fill = .data[[group]],
                        linetype = NA, group = paste(c(group, .width))),
                    alpha = .1)
      }
  }

  periods <- c(405, 132.5, 124, 99.7, 95)
  if (logx && logy) {
    pl <- pl + annotation_logticks() +
      scale_y_log10() +
      scale_x_log10(sec.axis = sec_axis(breaks = periods,
                                        trans = \(x) 1 / x, name = "Period (m)"))
  } else if (logx){
      pl <- pl + annotation_logticks(sides = "b") +
        scale_x_log10(sec.axis = sec_axis(breaks = periods,
                                          trans = \(x) 1 / x, name = "Period (m)"))
  } else if (logy) {
    pl <- pl + annotation_logticks(sides = "l") + scale_y_log10() +
      scale_x_continuous(sec.axis = sec_axis(breaks = periods,
                                             trans = \(x) 1 / x, name = "Period (m)"))
  } else {# neither
    pl <- pl +
      scale_x_continuous(sec.axis = sec_axis(breaks = periods,
                                             trans = \(x) 1 / x, name = "Period (m)"))
  }

  pl <- pl +
    geom_line() +
    labs(x = "Frequency (cycles/m)", y = "Spectral power (-)")
  return(pl)
}
