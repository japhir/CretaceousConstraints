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
