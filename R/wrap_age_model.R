#' Wrap Age Model
#'
#' @param data data.frame(). The data must contain columns `depth`, `ecc`, and `Ma405`.
#' @param agemodel data.frame() The agemodel must contain columns `sol`, `n`, `strat_bot`, and `age`.
#' @param astronomical_solution data.frame() The astronomical solution to calculate RMSD scores against. Must have column `sln`, `age`, and `scl`.
#' @param tiepoint_uncertainty Numeric vector of tiepoint uncertainty variations to try out in same units as data$depth.
#' @param genplot logical(1) Whether or not to generate a plot of the RMSD scores for each tiepoint uncertainty.
#' @param output character(1) Which output do you want the function to return? See details.
#'
#' @details
#' Output must be one of `"default"`, `"details"`, `"plot"`, or `"full"`.
wrap_age_model <- function(data,
                           agemodel,
                           astronomical_solution = sln,
                           tiepoint_uncertainty = seq(-4, 4, .5),
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

  # end input validation

  # make sure to fully initialize the whole output dataframe
  # this will make looping less slow
  the_best_summary <- agemodel |>
    filter(n %in% tiepoints) |>
    select(all_of(c("sol", "n", "strat_bot", "age"))) |>
    mutate(tie_err = NA_real_, # what's the best tiepoint error in m?
           RMSD_cum = NA_real_) # best RMSD score for full record

  the_best <- the_best_summary |>
    select(-tie_err, -RMSD_cum) |>
    mutate(
      optimal = list(tibble(error = tiepoint_uncertainty,
                            RMSD_tie = NA_real_))
    ) |>
    unnest(cols = c(optimal))

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
          # TODO: make frequencies a parameter!!!
          bandpass_filter(frequencies = my_filt_age |> filter(target != "23 kyr"),
                          x = age, y = value, add_depth = TRUE)

        ecc <- flt |>
          # TODO: sign should be passed based on proxy, -1 for d13C and Lstar, +1 for MS
          # TODO: pass in different weights, previously comb.
          construct_eccentricity(id_cols = c(depth, age, value), f = filter,
                                 sign = 1, weights = c(1, 1))

        esd <- ecc |>
          # linearly interpolate the astronomical solution eccentricity
          mutate(ecc_sln = approx(astronomical_solution$age,
                                  astronomical_solution$scl,
                                  xout = age)$y) |>
          # calculate SD between ecc and ecc_sol
          mutate(SD = (ecc - ecc_sln)^2) #|>

        # summarize into RMSD
        smy <- esd |>
          summarize(RMSD = sqrt(mean(SD)))

        # save the results to our tibble
        the_best$RMSD_tie[the_best$n == tiepoint & the_best$error == error] <- smy$RMSD

      } # end error loop

      tb <- the_best[the_best$n == tiepoint, ]
      bst <- min(tb$RMSD_tie)
      tb <- tb[tb$RMSD_tie == bst, ]

      # append the best result
      the_best_summary$tie_err[the_best_summary$n == tiepoint] <- tb$error
        # this clocks in at 81.6 µs, while dplyr alternatives below at 830 µs
#& the_best$RMSD_tie == min(the_best$RMSD_tie), "error"]##  |>
        ## filter(n == tiepoint) |>
        ## filter(RMSD_tie == min(RMSD_tie)) |>
        ## pull(error)

      the_best_summary$RMSD_cum[the_best_summary$n == tiepoint] <- tb$RMSD_tie
        ## the_best[the_best$n == tiepoint & RMSD_tie == min(RMSD_tie), "RMSD_tie"]##  |>
        ## filter(n == tiepoint) |>
        ## filter(RMSD_tie == min(RMSD_tie)) |>
        ## pull(RMSD_tie)

    } # end if
  } # end tiepoint loop

  if (genplot) {
    # faceted plot
    ## pl <- the_best |>
    ##   ggplot(aes(x = error, y = RMSD_tie, group = n)) +
    ##   facet_wrap(facets = vars(n)) +
    ##   geom_line() +
    ##   geom_point(aes(x = tie_err, y = RMSD_cum), colour = "red", size = 3,
    ##              data = the_best_summary)
    # new plot against depth
    pl <- the_best |>
      ggplot(aes(x = strat_bot + error, y = RMSD_tie)) +
      labs(x = "Depth (m)", y = "RMSD") +
      scale_x_reverse() +
      geom_vline(aes(xintercept = strat_bot), colour = "gray",
               data = the_best_summary) +
      geom_line(aes(group = n)) +
      geom_point(aes(x = strat_bot + tie_err, y = RMSD_cum),
                 colour = "red",
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
    return(list(summary = the_best_summary, details = the_best, plot = pl))
  }
} # end function definition
