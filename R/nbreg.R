#' Control parameters for `nbreg()`
#'
#' @param tol Convergence tolerance for the outer loop: stop when the
#'   log-likelihood changes by less than `tol` between iterations.
#' @param maxit Maximum number of outer iterations (alternating the IRLS
#'   update for beta and the Newton-Raphson update for theta).
#' @param irls_tol,irls_maxit Tolerance and maximum iterations for each inner
#'   IRLS fit of beta (theta fixed).
#' @param theta_tol,theta_maxit Tolerance and maximum iterations for each inner
#'   Newton-Raphson fit of log(theta) (beta fixed).
#' @param theta_max Warn if the estimated theta exceeds this value (the model
#'   is then close to Poisson).
#' @param trace Logical; if `TRUE`, print the log-likelihood and theta at each
#'   outer iteration.
#'
#' @return A list of control parameters.
#' @examples
#' nbreg_control(tol = 1e-10, trace = TRUE)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100,
                          irls_tol = 1e-10, irls_maxit = 100,
                          theta_tol = 1e-10, theta_maxit = 100,
                          theta_max = 1e6, trace = FALSE) {
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0)
    stop("'tol' must be a single positive number.", call. = FALSE)
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1)
    stop("'maxit' must be a single number >= 1.", call. = FALSE)
  list(tol = tol, maxit = as.integer(maxit),
       irls_tol = irls_tol, irls_maxit = as.integer(irls_maxit),
       theta_tol = theta_tol, theta_maxit = as.integer(theta_maxit),
       theta_max = theta_max, trace = isTRUE(trace))
}

#' Negative binomial regression by IRLS and Newton-Raphson
#'
#' Fits an NB2 regression model with log link,
#' \eqn{\log \mu_i = x_i^\top \beta + o_i} and
#' \eqn{\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta}, by alternating
#' iteratively reweighted least squares (IRLS) for \eqn{\beta} (with
#' \eqn{\theta} fixed) and Newton-Raphson for \eqn{\log\theta} (with
#' \eqn{\beta} fixed), until the log-likelihood converges. This is the same
#' model and parameterization as `MASS::glm.nb()`. No general-purpose
#' optimizer is used.
#'
#' @param formula A model formula, e.g. `y ~ x1 + x2`. Offsets may be given
#'   in the formula with `offset()`.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset vector (e.g. log exposure), added to the
#'   linear predictor.
#' @param control A list from [nbreg_control()].
#'
#' @return An object of class `"nbreg"`, a list with components including
#'   `coefficients`, `theta`, `loglik`, `iter`, `converged`,
#'   `fitted.values`, `linear.predictors`, `weights` (final IRLS working
#'   weights), `y`, `x`, `offset`, `terms`, and `call`.
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- rnorm(n)
#' y <- rnbinom(n, size = 2, mu = exp(0.5 + 0.8 * x))
#' fit <- nbreg(y ~ x, data = data.frame(y = y, x = x))
#' fit$coefficients
#' fit$theta
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control()) {
  cl <- match.call()

  # Build the model frame the same way glm() does, so that both an `offset`
  # argument and offset() terms in the formula are picked up.
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mt, mf)
  off <- as.vector(stats::model.offset(mf))
  if (is.null(off)) off <- rep(0, length(y))

  if (is.null(y))
    stop("The formula must have a response, e.g. y ~ x.", call. = FALSE)
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8))
    stop("The response must be non-negative integer counts.", call. = FALSE)
  if (length(off) != length(y))
    stop("'offset' must have one value per observation.", call. = FALSE)

  fit <- nbreg_fit(X, y, offset = off, control = control)

  structure(c(fit, list(y = y, x = X, offset = off, terms = mt,
                        call = cl)),
            class = "nbreg")
}

# Outer loop for NB regression (Section 3.1(c)).
#
# 1. Start from a Poisson IRLS fit (theta = Inf in nb_irls_beta).
# 2. Initial theta from the method of moments: Var(Y) = mu + mu^2/theta gives
#      theta_0 = sum a_i mu_i^2 / sum a_i [(y_i - mu_i)^2 - mu_i]
#    (if the denominator is <= 0 there is no overdispersion; start large).
# 3. Alternate (a) IRLS for beta with theta fixed and (b) Newton-Raphson for
#    log(theta) with beta fixed, until the log-likelihood changes by < tol.
#
# beta_start / theta_start: warm starts. The ZINB EM algorithm calls this as
# its count-part M-step with prior_weights = 1 - tau, starting from the
# previous EM iteration's (beta, theta) instead of a fresh Poisson fit.
# warn = FALSE suppresses the convergence / large-theta warnings (the EM loop
# reports its own convergence).
#' @noRd
nbreg_fit <- function(X, y, offset = rep(0, length(y)),
                      prior_weights = rep(1, length(y)),
                      control = nbreg_control(),
                      beta_start = NULL, theta_start = NULL, warn = TRUE) {
  if (is.null(beta_start)) {
    # Step 1: Poisson starting values
    pois <- nb_irls_beta(X, y, theta = Inf, offset = offset,
                         prior_weights = prior_weights,
                         tol = control$irls_tol, maxit = control$irls_maxit)
    beta <- pois$beta
    mu <- pois$mu
  } else {
    beta <- beta_start
    mu <- exp(drop(X %*% beta) + offset)
  }

  if (is.null(theta_start)) {
    # Step 2: moment estimate of theta
    denom <- sum(prior_weights * ((y - mu)^2 - mu))
    theta <- if (denom > 0) sum(prior_weights * mu^2) / denom else 1e4
  } else {
    theta <- theta_start
  }

  ll_old <- nb_loglik(y, mu, theta, prior_weights)
  converged <- FALSE

  # Step 3: alternate (a) and (b)
  for (iter in seq_len(control$maxit)) {
    # (a) IRLS for beta, theta fixed (warm start from current beta)
    b <- nb_irls_beta(X, y, theta, offset = offset,
                      prior_weights = prior_weights, beta_start = beta,
                      tol = control$irls_tol, maxit = control$irls_maxit)
    beta <- b$beta
    mu <- b$mu

    # (b) Newton-Raphson for log(theta), beta fixed
    th <- nb_theta_nr(y, mu, theta_start = theta,
                      prior_weights = prior_weights,
                      tol = control$theta_tol, maxit = control$theta_maxit)
    theta <- th$theta

    ll <- nb_loglik(y, mu, theta, prior_weights)
    if (control$trace)
      cat(sprintf("Outer iter %3d: loglik = %.10f, theta = %.6f\n",
                  iter, ll, theta))

    if (abs(ll - ll_old) < control$tol) {
      converged <- TRUE
      break
    }
    ll_old <- ll
  }

  if (warn && !converged)
    warning(sprintf("nbreg() did not converge in %d iterations.",
                    control$maxit), call. = FALSE)
  if (warn && theta > control$theta_max)
    warning(sprintf(paste0("Estimated theta (%.3g) is very large: the data ",
                           "show little or no overdispersion, so the model is ",
                           "close to Poisson."), theta), call. = FALSE)

  # IRLS working weights at the final (beta, theta)
  w <- prior_weights * theta * mu / (theta + mu)

  # Standard errors (Section 3.1(d)):
  #   Var(beta_hat) = (X' W_hat X)^{-1}, treating theta as known (as glm.nb
  #   does); X'WX is inverted through its Cholesky factor.
  #   SE(theta_hat) = {-H(theta_hat)}^{-1/2}, H from the theta Newton step.
  vcov <- chol2inv(chol(crossprod(X, w * X)))
  dimnames(vcov) <- list(names(beta), names(beta))
  H <- nb_theta_score(y, mu, theta, prior_weights)$H
  se_theta <- if (H < 0) (-H)^(-1/2) else NA_real_

  list(coefficients = beta, vcov = vcov, theta = theta, se_theta = se_theta,
       loglik = ll, iter = iter, converged = converged,
       fitted.values = mu, linear.predictors = drop(X %*% beta) + offset,
       weights = w, prior.weights = prior_weights)
}
