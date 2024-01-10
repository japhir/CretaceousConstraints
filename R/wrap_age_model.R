#' Wrap Age Model
#'
#' @param data data.frame(). The data must contain columns `depth`, `ecc`, and `Ma405`.
#' @param agemodel data.frame() The agemodel must contain columns `sol`, `n`, `strat_bot`, and `age`.
#' @param astronomical_solution data.frame() The astronomical solution to calculate RMSD scores against. Must have column `sln`, `age`, and `scl`.
#' @param tiepoint_uncertainty Numeric vector of tiepoint uncertainty variations to try out in same units as data$depth.
#' @param proxy_phase Factor by which to multiply the scaled and filtered signal to make it in-phase with the orbital forcing.
#' @param max_error_range In case the tiepoint_uncertainty is large, this limits which values are used to find the best value.
#' @param RMSD_threshold If the difference in RMSD is less than this value, just use the error = 0 value.
#' @param genplot logical(1) Whether or not to generate a plot of the RMSD scores for each tiepoint uncertainty.
#' @param output character(1) Which output do you want the function to return? See details.
#'
#' @details
#' Output must be one of `"default"`, `"details"`, `"plot"`, or `"full"`.
wrap_age_model <- function(data,
                           agemodel,
                           astronomical_solution = sln,
                           tiepoint_uncertainty = seq(-4, 4, .5),
                           proxy_phase = 1,
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

  valid_outputs <- c("default", "details", "plot", "full")
  # validate user inputs
  if (!output %in% valid_outputs) {
    cli::cli_abort(c(
      "{.var output} must be one of {.or {.q {valid_outputs}}}",
      "x" = "You've supplied {.q {output}}"))
  }

  if (length(unique(agemodel$sol)) > 1) {
    cli::cli_abort(c(
      "{.var agemodel} must contain only one unique astronomical solution in column {.var sol}.",
      "x" = "{.var agemodel} contains {.and {unique(agemodel$sol)}}"
    ))
  }

  if (length(unique(astronomical_solution$sln)) > 1) {
    cli::cli_abort(c(
      "{.var astronomical_solution} must contain only one unique astronomical solution in column {.var sln}.",
      "x" = "{.var astronomical_solution} contains {.and {unique(astronomical_solution$sln)}}"
    ))
  }

  # can't check this in data, sometimes not passed
  if (!(astronomical_solution$sln[[1]] == agemodel$sol[[1]])) {
    cli::cli_abort(c(
      "{.var sol} must be the same in {.var astronomical_solution} and {.var agemodel}",
      "i" = "{.var sol} in {.var agemodel} = {agemodel$sol[[1]]}",
      "i" = "{.var sln} in {.var astronomical_solution} = {astronomical_solution$sln[[1]]}"
    ))
  }

  if (output == "plot" & !genplot) {
    cli::cli_inform("{.var output} = {.q plot}, setting {.var genplot} = {.q TRUE}")
    genplot <- TRUE
  }

  if (!"Ma405" %in% colnames(data)) {
    data <- data |>
      mutate(Ma405 = findInterval(depth, agemodel$strat_bot))
    # should be correct, see
    # [[file:~/SurfDrive/Postdoc1/prj/2023-05-19_cretaceous_constraints/cretaceous_constraints.org::*Did
    # we number the Ma405 correctly?][Did we number the Ma405 correctly?]]
  }
  # this is kind of hard
  # while the above was correct, I want to /interpolate/, so should always
  # include one tiepoint lower and one higher.
  # if the tiepoint isn't in the agemodel (i.e. for d13C)
  # it'll skip it anyway.
  tiepoints <- c(min(data$Ma405) - 1,
                 unique(data$Ma405)#,
                 ## max(data$Ma405) + 1
                 ) |> sort()

  if (!0. %in% tiepoint_uncertainty) {
    cli::cli_abort(c(
      "{.q 0} must be present in {.var tiepoint_uncertainty}",
      "i" = "You have provided {.var tiepoint_uncertainty} = {tiepoint_uncertainty}."
    ))
  }
  # end input validation

  # filter the astronomical solution
  astronomical_solution_filters <- astronomical_solution |>
    bandpass_filter(frequencies = my_filt_age |>
                      filter(target != "23 kyr"),
                    x = age, y = ecc) |>
    mutate(filter = scale(filter)[, 1], .by = target) |>
    pivot_wider(id_cols = age,
                values_from = filter,
                names_from = target)

  # make sure to fully initialize the whole output dataframe
  # this will make looping less slow
  the_best_summary <- agemodel |>
    filter(n %in% tiepoints) |>
    select(all_of(c("sol", "n", "strat_bot", "age"))) |>
    mutate(tie_err = NA_real_, # what's the best tiepoint error in m?
           RMSD_cum = NA_real_,
           RMSD_cum_405 = NA_real_,
           RMSD_cum_100 = NA_real_) # best RMSD score for full record

  the_best <- the_best_summary |>
    select(-tie_err, -RMSD_cum, -RMSD_cum_405, -RMSD_cum_100) |>
    mutate(
      optimal = list(tibble(error = tiepoint_uncertainty,
                            RMSD_tie = NA_real_,
                            RMSD_tie_405 = NA_real_,
                            RMSD_tie_100 = NA_real_))
    ) |>
    unnest(cols = c(optimal))

  full_record <- NA_real_

  for (tiepoint in tiepoints) {
    if (!tiepoint %in% agemodel$n) {
      next
    } else {
      for (error in tiepoint_uncertainty) {

        if (!(tiepoint - 1) %in% the_best_summary$n) {
          am <- the_best_summary |>
            mutate(tie = case_when(n == tiepoint ~ error,
                                   .default = NA_real_)) |>
            mutate(depth = case_when(!is.na(tie) ~ strat_bot + tie,
                                     .default = strat_bot))

        } else {
          am <- the_best_summary |>
            mutate(tie = case_when(n == tiepoint ~ error,
                                   n < tiepoint ~ tie_err,
                                   .default = NA_real_)) |>
            mutate(depth = case_when(!is.na(tie) ~ strat_bot + tie,
                                     .default = strat_bot))
        }

        # calculate age from age model for each uncertainty
        tmp <- data |>
          mutate(age = map_dbl(depth,
                               ~ Hmisc::approxExtrap(
                                 # note that am |> pull(depth) clocks in at 266 µs, below at 1.25 µs!
                                 am$depth,
                                 am$age,
                                 xout = .x)$y),
                 .after = depth)

        flt <- tmp |>
          # TODO: make frequencies a parameter!!! For now, inheriting
          # my_filt_age from global namespace
          bandpass_filter(frequencies = my_filt_age |>
                            filter(target != "23 kyr"),
                          x = age, y = value, add_depth = TRUE)

        ecc <- flt |>
          # TODO: pass in different weights, previously comb.
          construct_eccentricity(id_cols = c(depth, age, value),
                                 f = filter,
                                 # I'm now forcing sign = 1 here!
                                 sign = 1, weights = c(1, 1)) |>
          # note that scaling occurs after
          mutate(`405 kyr` = scale(`405 kyr`)[, 1],
                 `100 kyr` = scale(`100 kyr`)[, 1])

        esd <- ecc |>
          # linearly interpolate the astronomical solution eccentricity
          mutate(
            ecc_sln_405 = approx(astronomical_solution_filters$age,
                                 astronomical_solution_filters$`405 kyr`,
                                 xout = age)$y,
            ecc_sln_100 = approx(astronomical_solution_filters$age,
                                 astronomical_solution_filters$`100 kyr`,
                                 xout = age)$y,
            ecc_sln = approx(astronomical_solution$age,
                             astronomical_solution$scl,
                             xout = age)$y) |>
          # calculate SD between ecc and ecc_sol
          mutate(
            SD_405 = (proxy_phase * `405 kyr` - ecc_sln_405)^2,
            SD_100 = (proxy_phase * `100 kyr` - ecc_sln_100)^2,
            SD = (proxy_phase * ecc - ecc_sln)^2)

        # DEBUG
        ## if (8 == tiepoint) {
        ##   browser()
        ## }

        ## # DEBUG FIGS
        ## (
        ##   esd |>
        ##   ggplot(aes(x = age * 1e-3,
        ##              y = value)) +
        ##   scale_x_reverse() +
        ##   geom_vline(xintercept = am$age * 1e-3) +
        ##     geom_line(colour = "skyblue") +
        ##     labs(title = "MS")
        ##   ## geom_line(aes(y = ecc_sln), colour = "purple")
        ## ) / (
        ##   esd |>
        ##   ggplot(aes(x = age * 1e-3,
        ##              y = proxy_phase * ecc)) +
        ##   scale_x_reverse() +
        ##   geom_vline(xintercept = am$age * 1e-3) +
        ##   geom_line(aes(y = ecc_sln), colour = "purple") +
        ##   geom_line(colour = "cyan", linewidth = 1.4)
        ## ) / (
        ## # DEBUG plot 405 kyr filters
        ##   esd |>
        ##   ggplot(aes(x = age * 1e-3,
        ##              y = proxy_phase * `405 kyr`)) +
        ##   scale_x_reverse() +
        ##   geom_vline(xintercept = am$age * 1e-3) +
        ##   geom_line(aes(y = ecc_sln_405), colour = "purple") +
        ##   geom_line(colour = "cyan", linewidth = 1.4)
        ## ) / (
        ## # DEBUG plot 100 kyr filters
        ## esd |>
        ##   ggplot(aes(x = age * 1e-3,
        ##              y = proxy_phase * `100 kyr`)) +
        ##   scale_x_reverse() +
        ##   geom_vline(xintercept = am$age * 1e-3) +
        ##   geom_line(aes(y = ecc_sln_100), colour = "purple") +
        ##   geom_line(colour = "cyan", linewidth = 1.4)
        ## ) +
        ##   plot_layout(axes = "collect") & labs(x = "Age (Ma)", y = "Filtered proxy")

        # annotate which depth tiepoints were originally closest to the am.
        ## geom_point(data = \(x)
        ##            x |>
        ##              filter(
        ##                near(depth, am$depth[[1]],
        ##                     tol = 3) |
        ##                  near(depth, am$depth[[2]],
        ##                       tol = 1) |
        ##                  near(depth, am$depth[[3]],
        ##                       tol = 0.5) |
        ##                  near(depth, am$depth[[4]],
        ##                       tol = 1) |
        ##                  near(depth, am$depth[[5]],
        ##                       tol = 1) |
        ##                  near(depth, am$depth[[6]],
        ##                       tol = 1) |
        ##                  near(depth, am$depth[[7]],
        ##                       tol = 1) |
        ##                  near(depth, am$depth[[8]],
        ##                       tol = 1) |
        ##                  near(depth, am$depth[[9]],
        ##                       tol = 0.5) |
        ##                  near(depth, am$depth[[10]],
        ##                       tol = .2) |
        ##                  near(depth, am$depth[[11]],
        ##                       tol = .2) |
        ##                  near(depth, am$depth[[12]],
        ##                       tol = .2) |
        ##                  near(depth, am$depth[[13]],
        ##                       tol = .2) |
        ##                  near(depth, am$depth[[14]],
        ##                       tol = .2) #|
        ##                     ),
        ##            colour = "red")

        # summarize into RMSD
        smy <- esd |>
          summarize(
            RMSD_405 = sqrt(mean(SD_405)),
            RMSD_100 = sqrt(mean(SD_100)),
            RMSD = sqrt(mean(SD))
          )

        # save the results to our tibble
        # NOTE: not using tidyverse syntax here because this is way faster than
        # e.g. pull.
        the_best$RMSD_tie[the_best$n == tiepoint &
                            the_best$error == error] <- smy$RMSD
        the_best$RMSD_tie_405[the_best$n == tiepoint &
                                the_best$error == error] <- smy$RMSD_405
        the_best$RMSD_tie_100[the_best$n == tiepoint &
                                the_best$error == error] <- smy$RMSD_100

        # before esd goes out scop!
        if (output == "full") {
          if (tiepoint == last(tiepoints) &&
                error == last(tiepoint_uncertainty))
            full_record <- esd
        }

      } # end error loop

      # find the best-scoring tiepoint error
      tb <- the_best[the_best$n == tiepoint, ]

      # limit the error options for choosing the best value
      # we do this here so that we can still compute it for a very wide range
      # so we can illustrate subsequent dips!
      tb <- tb |> filter(abs(error) <= max_error_range)

      # THIS is where we decide which strategy to use!
      bst <- min(tb$RMSD_tie)

      # if the difference between the best and 0 error is small
      if (bst > tb$RMSD_tie[tb$error == 0] - RMSD_threshold) {
        # just use the tiepoint from the field.
        bst <- tb$RMSD_tie[tb$error == 0]
      }
      tb <- tb[tb$RMSD_tie == bst, ]

      # append the best result
      the_best_summary$tie_err[the_best_summary$n == tiepoint] <- tb$error
      the_best_summary$RMSD_cum[the_best_summary$n == tiepoint] <- tb$RMSD_tie
      the_best_summary$RMSD_cum_405[the_best_summary$n == tiepoint] <- tb$RMSD_tie_405
      the_best_summary$RMSD_cum_100[the_best_summary$n == tiepoint] <- tb$RMSD_tie_100

      # DEBUG
      ## if (8 == tiepoint) {
      ##   browser()
      ## }

    } # end if
  } # end tiepoint loop

  if (genplot) {
    # new plot against depth
    pl <- the_best |>
      ggplot(aes(x = strat_bot + error, y = RMSD_tie,
                 colour = factor(n))) +
      labs(x = "Depth (m)", y = "RMSD") +
      scale_x_reverse() +
      geom_vline(aes(xintercept = strat_bot), colour = "gray",
                 data = the_best_summary) +
      geom_text(aes(x = strat_bot + tie_err, y = RMSD_cum, label = n),
                nudge_y = -5 * RMSD_threshold,
                data = the_best_summary) +
      annotate("segment", x = 0, xend = 0,
               y = 1.2, yend = 1.2 + RMSD_threshold) +
      geom_line(aes(group = n)) +
      geom_line(aes(group = n), linewidth = 1.2,
                data = \(x) x |> filter(abs(error) <= max_error_range)) +
      geom_point(aes(x = strat_bot + tie_err, y = RMSD_cum),
                 colour = "red", size = 3,
                 data = the_best_summary)
    print(pl)
  }

  if (output == "default") {
    return(the_best_summary)
  }
  if (output == "details") {
    return(the_best)
  }
  if (output == "plot") {
    return(pl)
  }
  if (output == "full") {
    return(list(summary = the_best_summary, details = the_best,
                full = full_record, plot = pl))
  }
} # end function definition
