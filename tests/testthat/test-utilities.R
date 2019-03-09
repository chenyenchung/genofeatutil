context("test-utilities")

test_that("mtry series", {
  df <- data.frame("x1" = c(1:10),
                   "x2" = c(1:10),
                   "x3" = c(1:10),
                   "x4" = c(1:10),
                   "x5" = c(1:10),
                   "x6" = c(1:10),
                   "x7" = c(1:10),
                   "x8" = c(1:10),
                   "x9" = c(1:10),
                   "yf" = factor(c(1:10)),
                   "yn" = c(1:10)
                   )
  expect_equal(make_mtry_series(x = df[ , 1:9], y = df[ , 10]), c(1, 3, 6, 9))
  expect_equal(make_mtry_series(x = df[ , 1:4], y = df[ , 11]), c(1, 2, 4))
})
