#' Wrap Age Model
#'
#' @param data data.frame(). The data must contain columns `depth`, `ecc`, and
#'   `Ma405`.
#' @param agemodel data.frame() The agemodel must contain columns `n`,
#'   `strat_bot`, and `age`.
#' @param astronomical_solution data.frame() The astronomical solution to
#'   calculate RMSD scores against. Must have column `sln`, `age`, and `scl`.
#' @param ... Other arguments, currently nothing just to force explicit argument names.
#' @param tiepoint_uncertainty Numeric vector of tiepoint uncertainty
#'   variations to try out in same units as data$depth.
#' @param proxy_phase Factor by which to multiply the scaled and filtered
#'   signal to make it in-phase with the orbital forcing.
#' @param target_periods Vector of periods to target. Defaults to 405 and 110
#'   kyr.
#' @param bandpass_window Window type for bandpass filter: 0 = rectangular, 1 = Gaussian, 2 = Cosine-tapered window (a.k.a. Tukey window).
#' @param frequency_fraction Fraction of frequency to adjust upper and lower
#'   filter limits by.
#' @param eccentricity_weights Weights for the 405 kyr and 100 kyr filter when
#'   creating the 'eccentricity' curve.
#' @param kpg_age The age of the K/Pg boundary to start with.
#' @param age_slider The age in kyr to slide the record by prior to matching to the solution.
#' @param max_error_range In case the tiepoint_uncertainty is large, this
#'   limits which values are used to find the best value.
#' @param RMSD_threshold If the difference in RMSD is less than this value,
#'   just use the error = 0 value.
#' @param genplot logical(1) Whether or not to generate a plot of the RMSD
#'   scores for each tiepoint uncertainty.
#' @param output character(1) Which output do you want the function to return?
#'   See details.
#' @return Depends on the output parameter. Either a dataframe, a ggplot
#'   object, or a list of dataframes.
#' @details Output must be one of `"default"`, `"details"`, `"plot"`,
#'   `"matched"`, or `"full"`.
#' @export
# TODO: add example!
wrap_age_model <- function(data,
                           agemodel,
                           astronomical_solution,
                           ...,
                           tiepoint_uncertainty = seq(-4, 4, .5),

                           proxy_phase = 1,
                           target_periods = c("405 kyr" = 405, "100 kyr" = 110),
                           frequency_fraction = 0.3,
                           bandpass_window = 0,
                           eccentricity_weights = c(1, 1),
                           age_slider = 200,
                           kpg_age = 65.9e3,
                           max_error_range = 2.5,
                           RMSD_threshold = 0.005,
                           genplot = FALSE,
                           output = "default") {

  data_input_columns <- c("depth", "value")
  if (!all(data_input_columns %in% colnames(data))) {
    cli::cli_abort(c(
      "{.var data} must contain columns {.and {.q {data_input_columns}}}",
      "x" = "You've supplied a dataframe with columns {.and {.q {colnames(data)}}}"))
  }
  if ("sol" %in% colnames(data)) {
    if (length(unique(data$sol)) > 1) {
      cli::cli_abort(c(
        "{.var data} must contain only one unique astronomical solution in column {.var sol}.",
        "x" = "{.var data} contains {.and {unique(data$sol)}}"
      ))
    }
  }
  if (nrow(data) == 0) {
    cli::cli_abort(c(
      "{.var data} must contain at least one row.",
      "x" = "{.var data} contains 0 rows."
    ))
  }

  valid_outputs <- c("default", "details", "matched", "plot", "full")
  # validate user inputs
  if (!output %in% valid_outputs) {
    cli::cli_abort(c(
      "{.var output} must be one of {.or {.q {valid_outputs}}}",
      "x" = "You've supplied {.q {output}}"))
  }

  if (length(unique(astronomical_solution$sln)) > 1) {
    cli::cli_abort(c(
      "{.var astronomical_solution} must contain only one unique astronomical solution in column {.var sln}.",
      "x" = "{.var astronomical_solution} contains {.and {unique(astronomical_solution$sln)}}"
    ))
  }

  solution_timestep <- astronomical_solution$age |>
    diff() |>
    unique()

  if (length(solution_timestep) > 1) {
    if (dplyr::near(solution_timestep[[1]], solution_timestep[[2]])) {
      # typically the AS has a step of 0.4
      # calculating diffs between this results in some double precision errors
      cli::cli_warn("The diff in the {.var astronomical_solution} timestep is not unique, but within machine precision.")
      solution_timestep <- solution_timestep |> round(2) |> unique()
    } else {
      cli::cli_abort(c("{.var astronomical_solution} must have a simple linear timestep!"))
    }
  }

  if (output %in% c("plot", "full") && !genplot) {
    cli::cli_inform("{.var output} is either {.q plot} or {.q full}, setting {.var genplot} = {.q TRUE}")
    genplot <- TRUE
  }

  if (!"Ma405" %in% colnames(data)) {
    data <- data |>
      dplyr::mutate(Ma405 = findInterval(.data$depth, agemodel$strat_bot))
  }
  # this is kind of hard
  # while the above was correct, I want to /interpolate/, so should always
  # include one tiepoint lower and one higher.
  # if the tiepoint isn't in the agemodel (i.e. for d13C)
  # it'll skip it anyway.
  # Note that we should NOT include the +1, I checked.
  tiepoints <- c(min(data$Ma405) - 1, unique(data$Ma405)) |>
    sort()
  # for Sopelana, the 14the tiepoint is NOT in the agemodel
  tiepoints <- tiepoints[tiepoints %in% agemodel$n]

  if (!0. %in% tiepoint_uncertainty) {
    cli::cli_abort(c(
      "{.q 0} must be present in {.var tiepoint_uncertainty}",
      "i" = "You have provided {.var tiepoint_uncertainty} = {tiepoint_uncertainty}."
    ))
  }

  if (length(target_periods) != 2) {
    cli::cli_abort(c(
      "Length of {.var target_periods} must be 2.",
      "i" = "You have provided {.var target_periods} = {target_periods}."
    ))
  }
  # end input validation

  # convert input frequency to tibble with desired frequency ranges
  my_filt_age <- tibble::tibble(target = c("405 kyr", "100 kyr"),
                        p = target_periods) |>
    dplyr::mutate(f = 1 / .data$p,
           range = frequency_fraction * .data$f,
           flow = .data$f - .data$range,
           fhigh = .data$f + .data$range,
           ref = "This study")

  # first calculate reference RMSD vs. age model with 0 error on the tiepoints
  tmp <- data |>
    dplyr::mutate(
             age_floating = purrr::map_dbl(.data$depth,
                                           ~ Hmisc::approxExtrap(
                                                      agemodel$strat_bot,
                                                      agemodel$age_floating,
                                                      xout = .x)$y),
                  .after = .data$depth)

  # data median timestep in kyr
  data_dt <- median(diff(tmp$age_floating))

  data_dt_sd <- sd(diff(tmp$age_floating))
  if (data_dt_sd > 0.0012) {
    cli::cli_warn(c(
           "The sd of the timestep after applying the age model is {data_dt_sd}.",
           "This is larger than 0.0012 (Zumaia < 109.26 m = 0.001149)"
         ))
  }

  # default resolution, 0.4 kyr -> 1.6 kyr
  ## linterp_dt <- solution_timestep * 4

  # round interpolation dt to nearest 0.4 interval
  ## linterp_dt <- solution_timestep * round(data_dt / solution_timestep)
  # for Zumaia depth < 109.26 data_dt = 2.7 ->
  # this results in 7, so 7 * 0.4 = 2.8
  linterp_dt <- solution_timestep * ceiling(data_dt / solution_timestep)

  flt <- tmp |>
    bandpass_filter(frequencies = my_filt_age,
                    x = .data$age_floating, y = .data$value,
                    window = bandpass_window,
                    linterp_dt = linterp_dt,
                    add_depth = TRUE)

  ## cli::cli_inform(c(
  ##        "nrow(tmp) = {nrow(tmp)}",
  ##        "nrow(flt) = {nrow(flt)}"
  ##      ))

  ecc <- flt |>
    # NOTE: this assumes that the frequencies tibble had column target
    # with names `405 kyr` and `100 kyr`!
    construct_eccentricity(id_cols = c(.data$depth,
                                       .data$age_floating,
                                             .data$value),
                           f = .data$filter,
                           # I'm now forcing sign = 1 here!
                           sign = 1,
                           weights = eccentricity_weights)

  # ensure that 0 is part of this!
  age_error <- seq(0, age_slider, linterp_dt)
  age_error <- c(-rev(age_error[-1]), age_error)

  # then optimize tiepoint uncertainty
  esd <- ecc |>
    dplyr::mutate(kpg_age = kpg_age, .after = .data$age_floating) |>
    dplyr::mutate(age_slider = list(age_error), .after = .data$kpg_age) |>
    tidyr::unnest(.data$age_slider) |>
    dplyr::mutate(age = .data$kpg_age +
                    .data$age_slider + .data$age_floating,
                  .after = .data$age_slider) |>
    # linearly interpolate the astronomical solution eccentricity
    dplyr::mutate(
             ecc_sln = stats::approx(astronomical_solution$age,
                                     astronomical_solution$scl,
                                     xout = .data$age)$y) |>
    # calculate SD between ecc and ecc_sol
    dplyr::mutate(SD = (proxy_phase * .data$ecc - .data$ecc_sln)^2,
                  .by = age_slider)

  # summarize into RMSD
  smy <- esd |>
    dplyr::summarize(RMSD = sqrt(mean(.data$SD)),
                     .by = age_slider)

  optimal_age_slider <- smy |>
    # get only the best one
    dplyr::filter(RMSD == min(RMSD)) |>
    pull(age_slider)

  if (length(optimal_age_slider) > 1) {
    cli::cli_abort(c(
           "Somehow {.var optimal_age_slider} has more than one element...",
           "{optimal_age_slider}"
         ))
  }

  # initialize the whole output dataframe
  # this will make looping less slow
  the_best_summary <- agemodel |>
    dplyr::filter(.data$n %in% tiepoints) |>
    dplyr::select(tidyr::all_of(c("n", "strat_bot", "age_floating"))) |>
    dplyr::mutate(kpg_age = kpg_age,
                  age_slider = optimal_age_slider,
                  tie_err = NA_real_, # what's the best tiepoint error in m?
                  RMSD_cum = NA_real_) # best RMSD score for full record

  the_best <- the_best_summary |>
    dplyr::select(-.data$tie_err, -.data$RMSD_cum) |>
    dplyr::mutate(
      optimal = list(tibble::tibble(error = tiepoint_uncertainty,
                                    RMSD_tie = NA_real_))
    ) |>
    tidyr::unnest(cols = c(.data$optimal))

  full_record <- esd |>
    dplyr::filter(.data$age_slider == optimal_age_slider)

  for (tiepoint in tiepoints) {
    if (!tiepoint %in% agemodel$n) {
      next
    } else {
      for (error in tiepoint_uncertainty) {

        if (!(tiepoint - 1) %in% the_best_summary$n) {
          am <- the_best_summary |>
            dplyr::mutate(tie = dplyr::case_when(
              .data$n == tiepoint ~ error,
              .default = NA_real_)) |>
            dplyr::mutate(depth = dplyr::case_when(
              !is.na(.data$tie) ~ .data$strat_bot + .data$tie,
              .default = .data$strat_bot))

        } else {
          am <- the_best_summary |>
            dplyr::mutate(tie = dplyr::case_when(
              .data$n == tiepoint ~ error,
              .data$n < tiepoint ~ tie_err,
              .default = NA_real_)) |>
            dplyr::mutate(depth = dplyr::case_when(
              !is.na(.data$tie) ~ .data$strat_bot + .data$tie,
              .default = .data$strat_bot))
        }

        # calculate age from age model for each uncertainty
        tmp <- data |>
          dplyr::mutate(age_floating = purrr::map_dbl(.data$depth,
                               ~ Hmisc::approxExtrap(
                                 # note that am |> pull(depth) clocks in at 266 µs, below at 1.25 µs!
                                 am$depth,
                                 am$age_floating,
                                 xout = .x)$y),
                 .after = .data$depth)

        flt <- tmp |>
          bandpass_filter(frequencies = my_filt_age,
                          x = .data$age_floating, y = .data$value,
                          window = bandpass_window,
                          linterp_dt = linterp_dt,
                          add_depth = TRUE)

        ecc <- flt |>
          # NOTE: this assumes that the frequencies tibble had column target
          # with names `405 kyr` and `100 kyr`!
          construct_eccentricity(id_cols = c(.data$depth,
                                             .data$age_floating,
                                             .data$value),
                                 f = .data$filter,
                                 # I'm now forcing sign = 1 here!
                                 sign = 1,
                                 weights = eccentricity_weights)

        esd <- ecc |>
          # get these from the no-tiepoint uncertainty run above
          dplyr::mutate(kpg_age = kpg_age, age_slider = optimal_age_slider) |>
          dplyr::mutate(age = .data$kpg_age + .data$age_slider + .data$age_floating) |>
          # linearly interpolate the astronomical solution eccentricity
          dplyr::mutate(
            ecc_sln = stats::approx(astronomical_solution$age,
                                    astronomical_solution$scl,
                                    xout = .data$age)$y) |>
          # calculate SD between ecc and ecc_sol
          dplyr::mutate(SD = (proxy_phase * .data$ecc - .data$ecc_sln)^2)

        # summarize into RMSD
        smy <- esd |>
          dplyr::summarize(RMSD = sqrt(mean(.data$SD)))

        # save the results to our tibble
        # NOTE: not using tidyverse syntax here because this is way faster than
        # e.g. pull.
        the_best$RMSD_tie[the_best$n == tiepoint &
                            the_best$error == error] <- smy$RMSD

        # before esd goes out scope!
        if (output %in% c("full", "matched")) {
          if (tiepoint == tail(tiepoints, 1) &&
                error == tail(tiepoint_uncertainty, 1))
            full_record <- esd
        }

      } # end error loop

      # find the best-scoring tiepoint error
      tb <- the_best[the_best$n == tiepoint, ]

      # limit the error options for choosing the best value
      # we do this here so that we can still compute it for a very wide range
      # so we can illustrate subsequent dips!
      tb <- tb |> dplyr::filter(abs(error) <= max_error_range)

      # THIS is where we decide which strategy to use!
      bst <- min(tb$RMSD_tie)

      # if the difference between the best and 0 error is small
      if (RMSD_threshold != 0) {
        if (bst > tb$RMSD_tie[tb$error == 0] - RMSD_threshold) {
          # just use the tiepoint from the field.
          bst <- tb$RMSD_tie[tb$error == 0]
        }
      }
      tb <- tb[tb$RMSD_tie == bst, ]

      # append the best result
      the_best_summary$tie_err[the_best_summary$n == tiepoint] <- tb$error
      the_best_summary$RMSD_cum[the_best_summary$n == tiepoint] <- tb$RMSD_tie

      # DEBUG
      ## if (8 == tiepoint) {
      ##   browser()
      ## }

    } # end if
  } # end tiepoint loop

  if (genplot) {
    # new plot against depth
    pl <- the_best |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$strat_bot + .data$error, y = .data$RMSD_tie,
        colour = factor(.data$n))) +
      ggplot2::labs(x = "Depth (m)", y = "RMSD") +
      ggplot2::scale_x_reverse() +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$strat_bot),
        colour = "gray",
        data = the_best_summary) +
      ggplot2::geom_text(ggplot2::aes(x = .data$strat_bot + .data$tie_err,
                                      y = .data$RMSD_cum, label = .data$n),
                nudge_y = -5 * RMSD_threshold,
                data = the_best_summary) +
      ggplot2::annotate("segment", x = 0, xend = 0,
                        y = 1.2, yend = 1.2 + RMSD_threshold) +
      ggplot2::geom_line(ggplot2::aes(group = .data$n)) +
      ggplot2::geom_line(ggplot2::aes(group = .data$n), linewidth = 1.2,
                data = \(x) x |> dplyr::filter(abs(error) <= max_error_range)) +
      ggplot2::geom_point(ggplot2::aes(x = .data$strat_bot + .data$tie_err,
                                       y = .data$RMSD_cum),
                 colour = "red", size = 3,
                 data = the_best_summary)
    if (output != "full") {
      print(pl)
    }
  }

  if (output == "default") {
    return(the_best_summary)
  }
  if (output == "details") {
    return(the_best)
  }
  if (output == "matched") {
    return(full_record)
  }
  if (output == "plot") {
    return(pl)
  }
  if (output == "full") {
    return(list(summary = the_best_summary,
                details = the_best,
                matched = full_record,
                plot = pl))
  }
} # end function definition
