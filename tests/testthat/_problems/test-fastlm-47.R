# Extracted from test-fastlm.R:47

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "fastlm625", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(mtcars)
fit <- fastlm(mpg ~ wt + hp, data = mtcars)
expect_output(print(fit), "Coefficients")
expect_output(summary(fit), "R-squared")
