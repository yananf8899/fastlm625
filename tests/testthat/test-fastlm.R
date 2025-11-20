test_that("fastlm_qr works", {
  data(mtcars)
  X <- as.matrix(mtcars[, c("wt", "hp")])
  y <- mtcars$mpg
  fit <- fastlm_qr(X, y)

  expect_true(is.numeric(fit$coefficients))
  expect_true(fit$r.squared > 0 && fit$r.squared < 1)
})

test_that("fastlm_chol works", {
  data(mtcars)
  X <- as.matrix(mtcars[, c("wt", "hp")])
  y <- mtcars$mpg
  fit <- fastlm_chol(X, y)

  expect_true(is.numeric(fit$coefficients))
  expect_true(fit$r.squared > 0 && fit$r.squared < 1)
})

test_that("fastlm matches lm()", {
  data(mtcars)
  fit_lm <- lm(mpg ~ wt + hp, data = mtcars)
  fit_fast <- fastlm(mpg ~ wt + hp, data = mtcars)

  expect_equal(as.numeric(coef(fit_lm)), as.numeric(fit_fast$coefficients), tolerance = 1e-10)
  expect_equal(summary(fit_lm)$r.squared, fit_fast$r.squared, tolerance = 1e-10)
})

test_that("QR and Cholesky give same results", {
  data(mtcars)
  X <- as.matrix(mtcars[, c("wt", "hp")])
  y <- mtcars$mpg

  fit_qr <- fastlm_qr(X, y)
  fit_chol <- fastlm_chol(X, y)

  expect_equal(fit_qr$coefficients, fit_chol$coefficients, tolerance = 1e-10)
  expect_equal(fit_qr$r.squared, fit_chol$r.squared, tolerance = 1e-10)
})

test_that("print works", {
  data(mtcars)
  fit <- fastlm(mpg ~ wt + hp, data = mtcars)
  expect_output(print(fit), "Coefficients")
})

test_that("predict works", {
  data(mtcars)
  fit <- fastlm(mpg ~ wt + hp, data = mtcars)
  pred <- predict(fit)

  expect_equal(pred, fit$fitted.values)
  expect_equal(length(pred), nrow(mtcars))
})
