# S3 methods for "nbreg" objects (project description, Section 4.2).

#' Methods for `nbreg` fits
#'
#' Standard extractor and summary methods for objects returned by [nbreg()].
#'
#' @param object,x An object of class `"nbreg"` (or `"summary.nbreg"` for
#'   the print method of a summary).
#' @param type For `residuals()`: `"pearson"` (default),
#'   \eqn{(y_i - \hat\mu_i)/\sqrt{\hat\mu_i + \hat\mu_i^2/\hat\theta}}, or
#'   `"response"`, \eqn{y_i - \hat\mu_i}. For `predict()`: `"response"`
#'   (default) for \eqn{\hat\mu}, or `"link"` for the linear predictor
#'   \eqn{\hat\eta = \log\hat\mu}.
#' @param newdata Optional data frame for `predict()`. If omitted, predictions
#'   are for the data used in the fit.
#' @param offset For `predict()` with `newdata`: an optional offset vector
#'   (in addition to any `offset()` term in the formula).
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `coef()` the coefficient vector; `vcov()` its covariance matrix;
#'   `logLik()` an object of class `"logLik"` with `df` = number of
#'   coefficients + 1 (for \eqn{\theta}) and `nobs`; `fitted()` the fitted means;
#'   `residuals()` a residual vector; `nobs()` the number of observations;
#'   `predict()` a vector of predictions; `summary()` an object of class
#'   `"summary.nbreg"` containing the coefficient table.
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(300))
#' d$y <- rnbinom(300, size = 2, mu = exp(0.5 + 0.8 * d$x))
#' fit <- nbreg(y ~ x, data = d)
#' fit
#' summary(fit)
#' coef(fit)
#' AIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#' @name nbreg-methods
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Negative binomial (NB2) regression, log link, fitted by IRLS + Newton-Raphson\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat(sprintf("\nTheta: %s   Log-likelihood: %s on %d df\n",
              format(x$theta, digits = digits),
              format(x$loglik, digits = digits + 2L),
              length(x$coefficients) + 1L))
  if (!x$converged) cat("Warning: the algorithm did not converge.\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- est / se
  # Wald tests: z = estimate / SE, two-sided p-value from N(0, 1)
  coef_table <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                      `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  ll <- stats::logLik(object)
  structure(list(call = object$call, coefficients = coef_table,
                 theta = object$theta, se_theta = object$se_theta,
                 loglik = object$loglik, df = attr(ll, "df"),
                 aic = stats::AIC(object), bic = stats::BIC(object),
                 iter = object$iter, converged = object$converged,
                 pearson = stats::residuals(object, type = "pearson")),
            class = "summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCoefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits)
  cat(sprintf("\nTheta: %s  (SE %s)\n", format(x$theta, digits = digits),
              format(x$se_theta, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d df   AIC: %s   BIC: %s\n",
              format(x$loglik, digits = digits + 2L), x$df,
              format(x$aic, digits = digits + 2L),
              format(x$bic, digits = digits + 2L)))
  cat(sprintf("Outer iterations (IRLS + Newton-Raphson): %d%s\n", x$iter,
              if (x$converged) "" else "  (NOT converged)"))
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  # df counts the regression coefficients plus theta
  structure(object$loglik, df = length(object$coefficients) + 1L,
            nobs = length(object$y), class = "logLik")
}

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  r <- object$y - mu
  # Pearson: divide by the NB2 standard deviation sqrt(mu + mu^2 / theta)
  if (type == "pearson") r <- r / sqrt(mu + mu^2 / object$theta)
  r
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) length(object$y)

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("response", "link"),
                          offset = NULL, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    # Rebuild the model matrix for the new data with the same terms and
    # factor levels as the fit
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                             xlev = object$xlevels)
    X <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    off <- stats::model.offset(mf)
    if (is.null(off)) off <- rep(0, nrow(X))
    if (!is.null(offset)) off <- off + offset
    eta <- drop(X %*% object$coefficients) + off
  }
  if (type == "link") eta else exp(eta)
}
