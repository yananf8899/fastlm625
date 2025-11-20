#' Fast Linear Regression using QR Decomposition
#'
#' Fits a linear regression model using QR decomposition.
#'
#' @param X Numeric matrix of predictors.
#' @param y Numeric vector of response values.
#' @param intercept Logical; whether to include an intercept (default = TRUE).
#'
#' @return A list of class \code{"fastlm"} containing:
#' \item{coefficients}{Estimated regression coefficients}
#' \item{fitted.values}{Fitted (predicted) values}
#' \item{residuals}{Residuals}
#' \item{sigma}{Estimated residual standard error}
#' \item{df.residual}{Residual degrees of freedom}
#' \item{r.squared}{R-squared}
#' \item{adj.r.squared}{Adjusted R-squared}
#' \item{vcov}{Variance-covariance matrix of coefficients}
#' \item{se}{Standard errors of coefficients}
#' @examples
#' data(mtcars)
#' X <- as.matrix(mtcars[, c("wt", "hp")])
#' y <- mtcars$mpg
#' fit <- fastlm_qr(X, y)
#' fit$coefficients
#' @export
fastlm_qr <- function(X, y, intercept = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- length(y)
  
  if (intercept) X <- cbind(Intercept = 1, X)
  p <- ncol(X)
  
  qr_obj <- qr(X)
  coef <- qr.coef(qr_obj, y)
  fitted <- qr.fitted(qr_obj, y)
  resid <- qr.resid(qr_obj, y)
  
  df <- n - p
  rss <- sum(resid^2)
  tss <- sum((y - mean(y))^2)
  r_squared <- 1 - rss / tss
  sigma <- sqrt(rss / df)
  
  R <- qr.R(qr_obj)
  R_inv <- backsolve(R, diag(p))
  vcov <- sigma^2 * tcrossprod(R_inv)
  se <- sqrt(diag(vcov))
  
  if (!is.null(colnames(X))) {
    names(coef) <- colnames(X)
    names(se) <- colnames(X)
  }
  
  result <- list(
    coefficients = coef,
    fitted.values = fitted,
    residuals = resid,
    sigma = sigma,
    df.residual = df,
    r.squared = r_squared,
    adj.r.squared = 1 - (1 - r_squared) * (n - 1) / df,
    vcov = vcov,
    se = se,
    call = match.call()
  )
  class(result) <- "fastlm"
  result
}

#' Fast Linear Regression using Cholesky Decomposition
#'
#' Fits a linear regression model using Cholesky decomposition
#' (usually faster than QR).
#'
#' @inheritParams fastlm_qr
#'
#' @examples
#' data(mtcars)
#' X <- as.matrix(mtcars[, c("wt", "hp")])
#' y <- mtcars$mpg
#' fit <- fastlm_chol(X, y)
#' fit$coefficients
#' @export
fastlm_chol <- function(X, y, intercept = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- length(y)
  
  if (intercept) X <- cbind(Intercept = 1, X)
  p <- ncol(X)
  
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  
  chol_XtX <- chol(XtX)
  coef <- backsolve(chol_XtX, forwardsolve(t(chol_XtX), Xty))
  coef <- as.vector(coef)
  
  fitted <- as.vector(X %*% coef)
  resid <- y - fitted
  
  df <- n - p
  rss <- sum(resid^2)
  tss <- sum((y - mean(y))^2)
  r_squared <- 1 - rss / tss
  sigma <- sqrt(rss / df)
  
  chol_inv <- backsolve(chol_XtX, diag(p))
  vcov <- sigma^2 * tcrossprod(chol_inv)
  se <- sqrt(diag(vcov))
  
  if (!is.null(colnames(X))) {
    names(coef) <- colnames(X)
    names(se) <- colnames(X)
  }
  
  result <- list(
    coefficients = coef,
    fitted.values = fitted,
    residuals = resid,
    sigma = sigma,
    df.residual = df,
    r.squared = r_squared,
    adj.r.squared = 1 - (1 - r_squared) * (n - 1) / df,
    vcov = vcov,
    se = se,
    call = match.call()
  )
  class(result) <- "fastlm"
  result
}

#' Fast Linear Regression using Formula Interface
#'
#' Fits a linear regression model from a formula and data frame.
#'
#' @param formula A formula like \code{y ~ x1 + x2}.
#' @param data A data frame.
#' @param method Character; either \code{"qr"} or \code{"chol"} (default \code{"qr"}).
#'
#' @examples
#' data(mtcars)
#' fit <- fastlm(mpg ~ wt + hp, data = mtcars)
#' summary(fit)
#' @return A \code{fastlm} object.
#' @export
fastlm <- function(formula, data, method = c("qr", "chol")) {
  method <- match.arg(method)
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0) {
    X <- X[, -intercept_col, drop = FALSE]
    intercept <- TRUE
  } else {
    intercept <- FALSE
  }
  
  if (method == "qr") {
    result <- fastlm_qr(X, y, intercept = intercept)
  } else {
    result <- fastlm_chol(X, y, intercept = intercept)
  }
  
  result$call <- match.call()
  result$formula <- formula
  result
}

#' Print method for fastlm objects
#'
#' @param x A \code{fastlm} model.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.fastlm <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nR-squared:", round(x$r.squared, 4), "\n\n")
  invisible(x)
}

#' Summary method for fastlm objects
#'
#' @param object A \code{fastlm} model.
#' @param ... Additional arguments (ignored).
#'
#' @return The input \code{fastlm} object invisibly.
#' @export
summary.fastlm <- function(object, ...) {
  t_stats <- object$coefficients / object$se
  p_values <- 2 * pt(abs(t_stats), df = object$df.residual, lower.tail = FALSE)
  
  coef_table <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `t value` = t_stats,
    `Pr(>|t|)` = p_values
  )
  
  cat("\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  printCoefmat(coef_table, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nResidual standard error:", round(object$sigma, 4),
      "on", object$df.residual, "degrees of freedom")
  cat("\nR-squared:", round(object$r.squared, 4),
      "Adjusted R-squared:", round(object$adj.r.squared, 4), "\n\n")
  invisible(object)
}

#' Print method for summary.fastlm objects
#'
#' @param x A summary of \code{fastlm} model.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.summary.fastlm <- function(x, ...) {
  invisible(x)
}

#' Predict method for fastlm objects
#'
#' @param object A \code{fastlm} model.
#' @param newdata Optional new data. If \code{NULL}, returns fitted values.
#' @param ... Additional arguments (ignored).
#'
#' @return A numeric vector of predictions.
#' @export
predict.fastlm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }
  
  if (!is.null(object$formula)) {
    terms_obj <- delete.response(terms(object$formula))
    mf <- model.frame(terms_obj, newdata)
    X <- model.matrix(terms_obj, mf)
  } else {
    if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
    X <- cbind(1, newdata)
  }
  
  as.vector(X %*% object$coefficients)
}
