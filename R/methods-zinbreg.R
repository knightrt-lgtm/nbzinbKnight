# S3 methods for "zinbreg" objects (project description, Section 4.2).

#' Methods for `zinbreg` fits
#'
#' Standard extractor and summary methods for objects returned by
#' [zinbreg()]. Coefficient vectors and covariance matrices use the labels
#' `count_*` and `zero_*`, as `pscl::zeroinfl()` does.
#'
#' @param object,x An object of class `"zinbreg"` (or `"summary.zinbreg"`
#'   for the print method of a summary).
#' @param model For `coef()` and `vcov()`: `"full"` (default, count and zero
#'   parts together), `"count"`, or `"zero"`.
#' @param type For `residuals()`: `"pearson"` (default),
#'   \eqn{(y_i - \hat E Y_i)/\sqrt{\widehat{\mathrm{Var}}(Y_i)}}, or
#'   `"response"`, \eqn{y_i - \hat E Y_i}, where
#'   \eqn{E Y_i = (1-\pi_i)\mu_i} and
#'   \eqn{\mathrm{Var}(Y_i) = (1-\pi_i)\mu_i\{1 + \mu_i(\pi_i + 1/\theta)\}}.
#'   For `predict()`: `"response"` (default) for the mean
#'   \eqn{(1-\hat\pi)\hat\mu}, `"count"` for the count-model mean
#'   \eqn{\hat\mu}, or `"zero"` for the zero-inflation probability
#'   \eqn{\hat\pi}.
#' @param newdata Optional data frame for `predict()`. If omitted, predictions
#'   are for the data used in the fit.
#' @param offset For `predict()` with `newdata`: an optional count-model
#'   offset vector (in addition to any `offset()` term in the formula).
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `coef()` a named coefficient vector; `vcov()` its covariance
#'   matrix; `logLik()` an object of class `"logLik"` with `df` = number of
#'   count and zero coefficients + 1 (for \eqn{\theta}) and `nobs`; `fitted()`
#'   the fitted means \eqn{(1-\hat\pi)\hat\mu}; `residuals()` a residual vector;
#'   `nobs()` the number of observations; `predict()` a vector of predictions;
#'   `summary()` an object of class `"summary.zinbreg"` with separate
#'   coefficient tables for the count and zero parts.
#' @examples
#' set.seed(1)
#' n <- 500
#' d <- data.frame(x = rnorm(n))
#' d$y <- ifelse(rbinom(n, 1, plogis(-1 + 0.5 * d$x)) == 1, 0,
#'               rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
#' fit <- zinbreg(y ~ x | x, data = d)
#' fit
#' summary(fit)
#' coef(fit, model = "zero")
#' AIC(fit)
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "zero")
#' @name zinbreg-methods
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Zero-inflated negative binomial regression, fitted by EM with IRLS M-steps\n\n")
  cat("Count model coefficients (negbin with log link):\n")
  print.default(format(x$coefficients$count, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat(sprintf("\nTheta: %s   Log-likelihood: %s on %d df\n",
              format(x$theta, digits = digits),
              format(x$loglik, digits = digits + 2L),
              attr(stats::logLik(x), "df")))
  if (!x$converged) cat("Warning: the EM algorithm did not converge.\n")
  invisible(x)
}

# Wald coefficient table: estimate, SE, z = estimate / SE, two-sided p-value
#' @noRd
wald_table <- function(est, se) {
  z <- est / se
  cbind(Estimate = est, `Std. Error` = se, `z value` = z,
        `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  p <- length(object$coefficients$count)
  count_tab <- wald_table(object$coefficients$count, se[seq_len(p)])
  # Log(theta) row in the count table, as zeroinfl() reports it
  theta_row <- wald_table(log(object$theta), object$se_logtheta)
  rownames(theta_row) <- "Log(theta)"
  count_tab <- rbind(count_tab, theta_row)
  zero_tab <- wald_table(object$coefficients$zero, se[-seq_len(p)])
  ll <- stats::logLik(object)
  structure(list(call = object$call,
                 coefficients = list(count = count_tab, zero = zero_tab),
                 theta = object$theta, se_theta = object$se_theta,
                 loglik = object$loglik, df = attr(ll, "df"),
                 aic = stats::AIC(object), bic = stats::BIC(object),
                 iter = object$iter, converged = object$converged,
                 pearson = stats::residuals(object, type = "pearson")),
            class = "summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCount model coefficients (negbin with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits)
  cat(sprintf("\nTheta: %s  (SE %s)\n", format(x$theta, digits = digits),
              format(x$se_theta, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d df   AIC: %s   BIC: %s\n",
              format(x$loglik, digits = digits + 2L), x$df,
              format(x$aic, digits = digits + 2L),
              format(x$bic, digits = digits + 2L)))
  cat(sprintf("EM iterations: %d%s\n", x$iter,
              if (x$converged) "" else "  (NOT converged)"))
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         full = stats::setNames(c(cf$count, cf$zero),
                                c(paste0("count_", names(cf$count)),
                                  paste0("zero_", names(cf$zero)))),
         count = cf$count,
         zero = cf$zero)
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  V <- object$vcov
  p <- length(object$coefficients$count)
  idx <- switch(model, full = seq_len(nrow(V)), count = seq_len(p),
                zero = seq_len(nrow(V))[-seq_len(p)])
  V <- V[idx, idx, drop = FALSE]
  if (model != "full") dimnames(V) <- lapply(dimnames(V), sub,
                                             pattern = "^(count|zero)_",
                                             replacement = "")
  V
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  # df: count coefficients + zero coefficients + theta
  structure(object$loglik,
            df = length(object$coefficients$count) +
              length(object$coefficients$zero) + 1L,
            nobs = length(object$y), class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  mu <- object$mu
  pi_ <- object$pi
  r <- object$y - (1 - pi_) * mu
  if (type == "pearson") {
    # Var(Y) = (1 - pi) mu [1 + mu (pi + 1/theta)] for the ZINB mixture
    r <- r / sqrt((1 - pi_) * mu * (1 + mu * (pi_ + 1 / object$theta)))
  }
  r
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) length(object$y)

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero"),
                            offset = NULL, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    mu <- object$mu
    pi_ <- object$pi
  } else {
    # Count part: rebuild X with the fit's terms and factor levels
    ttX <- stats::delete.response(object$terms$count)
    mfX <- stats::model.frame(ttX, newdata, na.action = stats::na.pass,
                              xlev = object$xlevels$count)
    X <- stats::model.matrix(ttX, mfX, contrasts.arg = object$contrasts$count)
    off <- stats::model.offset(mfX)
    if (is.null(off)) off <- rep(0, nrow(X))
    if (!is.null(offset)) off <- off + offset
    mu <- exp(drop(X %*% object$coefficients$count) + off)

    # Zero part: rebuild Z the same way
    ttZ <- stats::delete.response(object$terms$zero)
    mfZ <- stats::model.frame(ttZ, newdata, na.action = stats::na.pass,
                              xlev = object$xlevels$zero)
    Z <- stats::model.matrix(ttZ, mfZ, contrasts.arg = object$contrasts$zero)
    pi_ <- stats::plogis(drop(Z %*% object$coefficients$zero))
  }
  switch(type, response = (1 - pi_) * mu, count = mu, zero = pi_)
}
