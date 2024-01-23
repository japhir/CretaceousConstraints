#' Construct eccentricity from filtered record
#'
#' @param data data.frame() Result of [bandpass_filter()]. Must contain columns
#'   `depth`, `value`, `filter`, and `target`.
#' @param sign numeric(0) with the desired phase relationship for this proxy.
#-1 for d13C and Lstar, +1 for MS! #' @param x, y Column names in `data` that
#contain depth and values. #' @param weights Numeric vector of weights to apply
#for each target.
#' @param weights A vector of 2 numbers, the factors by which to multiply the
#'   405 kyr and 100 kyr filtered record.
#' @param id_cols The columns to maintain when pivoting wider.
#' @param f The column in data that contains the filter values.
#' @export
construct_eccentricity <- function(data, sign = 1, weights = c("405" = 1, "100" = 1),
                                   id_cols = c(.data$depth, .data$age, .data$value), f = .data$filter
                                   ) {
  data |>
    tidyr::pivot_wider(id_cols = {{id_cols}},
                       names_from = .data$target,
                       values_from = {{f}}
                ) |>
    dplyr::mutate(ecc = sign * weights[[1]] *
                    scale(.data$`405 kyr`)[, 1] +
                    sign * weights[[2]] *
                    scale(.data$`100 kyr`)[, 1])
}
