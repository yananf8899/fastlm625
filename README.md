# fastlm625

Fast linear regression using QR and Cholesky decomposition.

## Installation

```r
devtools::install_github("yananf8899/fastlm625")
```

## Usage

```r
library(fastlm625)

# Formula interface
data(mtcars)
fit <- fastlm(mpg ~ wt + hp, data = mtcars)
summary(fit)
predict(fit)[1:5]

# Matrix interface - QR decomposition
X <- as.matrix(mtcars[, c("wt", "hp")])
y <- mtcars$mpg
fit_qr <- fastlm_qr(X, y)
print(fit_qr$coefficients)

# Cholesky decomposition (faster)
fit_chol <- fastlm_chol(X, y)
print(fit_chol$coefficients)
```

## Features

- QR decomposition (numerically stable)
- Cholesky decomposition (faster)
- Formula interface like base `lm()`
- S3 methods: print, summary, predict
- Full regression statistics

## Performance

About 1.5-2x faster than base R `lm()` for medium to large datasets.
Cholesky is usually fastest when the design matrix is well-conditioned