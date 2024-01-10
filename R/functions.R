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
  if (logx & logy) {
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

evolutive_analysis <- function(data, nest, x, y) {
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
      eha = map(data,
                ~ . |>
                  select({{x}}, {{y}}) |>
                  astrochron::linterp(genplot = FALSE, verbose = FALSE) |>
                  astrochron::eha(output = 3, genplot = FALSE, verbose = FALSE) |>
                  pivot_longer(-freq) |>
                  mutate(name = str_sub(name, 2, -1) |> parse_double()) |>
                  rename({{x}} := name))
    ) |>
    select(-data) |>
    unnest(eha)
}

#' Bandpass filter
#'
#' @param data data.frame()
#' @param frequencies data.frame() with columns `target`, `flow`, and `fhigh`.
#'   Can contain multiple rows for multiple frequency filtering.
#' @param x Column name in `data` that holds the depth/age information.
#' @param y Column name in `data` that holds the proxy variable name.
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
                          arrange({{x}}) |>
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

#' Construct eccentricity from filtered record
#'
#' @param data data.frame() Result of [bandpass_filter()]. Must contain columns
#'   `depth`, `value`, `filter`, and `target`.
#' @param sign numeric(0) with the desired phase relationship for this proxy. -1 for d13C and Lstar, +1 for MS!
# #' @param x, y Column names in `data` that contain depth and values.
# #' @param weights Numeric vector of weights to apply for each target.
construct_eccentricity <- function(data, sign = 1, weights = c("405" = 1, "100" = 1),
                                   id_cols = c(depth, age, value), f = filter
                                   ) {
  data |>
    pivot_wider(id_cols = {{id_cols}},
                names_from = target,
                values_from = {{f}}
                ) |>
    mutate(ecc = sign * weights[[1]] * scale(`405 kyr`)[, 1] +
             sign * weights[[2]] * scale(`100 kyr`)[, 1])
}

get_rmcd <- function(data, rmcd = "dat/ODP208_1267_rmcd.csv") {
  rmcd <- readr::read_csv(rmcd) |>
    separate(label, into = c("sitehole", "coretype", "Sec"),
               sep = "-", remove = FALSE) |>
    separate(sitehole, into = c("Site", "H"), sep = -1) |>
    separate(coretype, into = c("Core", "T"), sep = -1) |>
    # we do not rename the interval, may not be the same as in the data!
    # rename the CC sections into 7, the naming convention in the MS data
    mutate(Sec = ifelse(Sec == "7", "7", Sec),
           Sec = ifelse(Sec == "cc", "C", Sec)) |>
    mutate(diff = depth_rmcd - depth_mbsf, .after = depth_rmcd) |>
    mutate(diff2 = depth_rmcd2 - depth_mbsf2, .after = depth_rmcd2) |>
    mutate(row = 1:n())

  # the right side of the splice table only
  rmcd2 <- rmcd |>
    select(label, link, label2, interval2, depth_mbsf2, depth_rmcd2, diff2, row) |>
    separate(label2, into = c("sitehole", "coretype", "Sec"),
               sep = "-", remove = FALSE) |>
    separate(sitehole, into = c("Site", "H"), sep = -1) |>
    separate(coretype, into = c("Core", "T"), sep = -1) |>
    # we do not rename the interval, may not be the same as in the data!
    mutate(Sec = ifelse(Sec == "7", "7", Sec),
           Sec = ifelse(Sec == "cc", "C", Sec))

  out <- data |>
    tidylog::left_join(rmcd |>
                     # make the types the same
                     mutate(across(c(Site, Core), parse_double)) |>
                     # do NOT match by section, only by core!
                     rename(section = Sec) |>
                     select(top = label, to = label2,
                            Site, H, Core, T, section, interval,
                            depth_mbsf, depth_rmcd, diff, row)) |>
    # add the right-hand side of the splice table
    tidylog::left_join(rmcd2 |>
                       mutate(across(c(Site, Core), parse_double)) |>
                       rename(section2 = Sec) |>
                       select(from = label, bot = label2,
                              Site, H, Core, T, section2, interval2,
                              depth_mbsf2, depth_rmcd2, diff2, row2 = row)) |>
  mutate(my_rmcd = case_when(
  (Sec <= section) | ((Sec == section) & (`Top (cm)` <= interval)) ~
    `Depth (mbsf)` + diff,
  (Sec >= section2) | ((Sec == section2) & (`Top (cm)` >= interval2)) ~
    `Depth (mbsf)` + diff2,
  TRUE ~ NA_real_)) |>
    mutate(on_splice = (Sec < section | ((Sec == section) &
                                         (`Top (cm)` <= interval))) &
             (Sec > section2 | ((Sec == section2) &
                                `Top (cm)` >= interval2))) |>
    mutate(on_splice = ifelse(is.na(on_splice), FALSE, on_splice))

  return(out)
}
