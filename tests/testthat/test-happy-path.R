test_that("bandpass_filter has a happy path", {
  set.seed(123)
  dat <- tibble::tibble(a = 1:10,
                      b = 11:20,
                      c = stats::rnorm(10),
                      d = sample(letters[1:3], 10, TRUE))
  expect_snapshot(bandpass_filter(dat,
                  tibble::tibble(flow = 1, fhigh = 2, target = "group"),
                  x = a, y = c,
                  window = 0))
})
