#' Control parameters for `zinbreg()`
#'
#' @param tol EM convergence tolerance: stop when the observed-data
#'   log-likelihood changes by less than `tol` between EM iterations.
#' @param maxit Maximum number of EM iterations.
#' @param nb_control A list from [nbreg_control()] used for the count-part
#'   M-step (the weighted NB fit).
#' @param logit_tol,logit_maxit Tolerance and maximum iterations for the
#'   zero-part M-step (the IRLS logistic regression).
#' @param trace Logical; if `TRUE`, print the observed-data log-likelihood at
#'   each EM iteration.
#'
#' @return A list of control parameters.
#' @examples
#' zinbreg_control(tol = 1e-10, maxit = 1000)
#' @export
zinbreg_control <- function(tol = 1e-8, maxit = 5000,
                            nb_control = nbreg_control(),
                            logit_tol = 1e-10, logit_maxit = 100,
                            trace = FALSE) {
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0)
    stop("'tol' must be a single positive number.", call. = FALSE)
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1)
    stop("'maxit' must be a single number >= 1.", call. = FALSE)
  list(tol = tol, maxit = as.integer(maxit), nb_control = nb_control,
       logit_tol = logit_tol, logit_maxit = as.integer(logit_maxit),
       trace = isTRUE(trace))
}

#' Zero-inflated negative binomial regression by EM with IRLS M-steps
#'
#' Fits a zero-inflated negative binomial (ZINB) model. With probability
#' \eqn{\pi_i} an observation is a structural zero; otherwise it comes from
#' an NB2 distribution with mean \eqn{\mu_i} and dispersion \eqn{\theta}.
#' The count model is \eqn{\log \mu_i = x_i^\top\beta + o_i} and the zero model
#' is \eqn{\mathrm{logit}\,\pi_i = z_i^\top\gamma}. The observed-data
#' log-likelihood is maximized by an EM algorithm: the E-step computes the
#' posterior probability \eqn{\tau_i} that each zero is structural, and the
#' M-steps are a logistic IRLS fit for \eqn{\gamma} (fractional response
#' \eqn{\tau_i}) and a weighted NB fit for \eqn{(\beta, \theta)} with prior
#' weights \eqn{1 - \tau_i}. No general-purpose optimizer is used.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before `|`, zero-model covariates after it. Without `|`, the
#'   same covariates are used in both parts.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset for the count model (e.g. log exposure).
#'   An `offset()` term in the count part of the formula also works.
#' @param control A list from [zinbreg_control()].
#'
#' @return An object of class `"zinbreg"`, a list with components including
#'   `coefficients` (a list with `count` and `zero`), `theta`, `loglik`,
#'   `loglik_trace` (log-likelihood at every EM iteration), `iter`,
#'   `converged`, `fitted.values` (\eqn{(1-\pi_i)\mu_i}), `mu`, `pi`, `tau`,
#'   `y`, `x`, `z`, `offset`, `terms`, `formula`, and `call`.
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- rnorm(n)
#' structural_zero <- rbinom(n, 1, plogis(-1 + 0.5 * x))
#' y <- ifelse(structural_zero == 1, 0,
#'             rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)))
#' fit <- zinbreg(y ~ x | x, data = data.frame(y = y, x = x))
#' fit$coefficients
#' fit$theta
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control()) {
  cl <- match.call()

  # Two-part formula y ~ x | z. Without "|", use the count covariates for the
  # zero model as well (as pscl::zeroinfl does).
  ff <- Formula::as.Formula(formula)
  if (length(ff)[2L] < 2L) {
    ff <- Formula::as.Formula(stats::formula(ff),
                              stats::formula(ff, lhs = 0L, rhs = 1L))
  }

  # Model frame (picks up the `offset` argument as "(offset)")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$formula <- ff
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  y <- stats::model.response(mf, "numeric")
  if (is.null(y))
    stop("The formula must have a response, e.g. y ~ x | z.", call. = FALSE)
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8))
    stop("The response must be non-negative integer counts.", call. = FALSE)
  if (!any(y == 0))
    stop("The response has no zeros, so a zero-inflated model cannot be fitted.",
         call. = FALSE)

  mtX <- stats::terms(ff, rhs = 1L)
  mtZ <- stats::terms(ff, rhs = 2L)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)

  # Count-model offset: offset() terms in the count formula plus the
  # `offset` argument
  off <- rep(0, length(y))
  for (i in attr(mtX, "offset")) {
    off <- off + mf[[deparse(attr(mtX, "variables")[[i + 1L]])]]
  }
  if (!is.null(mf[["(offset)"]])) off <- off + mf[["(offset)"]]

  fit <- zinbreg_fit(X, Z, y, offset = off, control = control)

  structure(c(fit, list(y = y, x = X, z = Z, offset = off,
                        terms = list(count = mtX, zero = mtZ),
                        formula = ff, call = cl)),
            class = "zinbreg")
}

# EM algorithm for ZINB (Section 3.2).
#
# Start:  (beta, theta) from the NB fit; gamma from an IRLS logistic
#         regression of I(y_i = 0) on z_i.
# Repeat:
#   E-step:  tau_i = P(structural zero | y_i)                  (zinb_estep)
#   M-step (zero part):  logistic IRLS of tau on Z -> gamma    (logit_irls)
#   M-step (count part): NB fit with prior weights 1 - tau_i -> (beta, theta)
#                        (nbreg_fit, warm-started)
# until the observed-data log-likelihood changes by < tol. EM must never
# decrease the log-likelihood; a decrease larger than rounding error gives a
# warning (it would indicate a bug).
#' @noRd
zinbreg_fit <- function(X, Z, y, offset = rep(0, length(y)),
                        control = zinbreg_control()) {
  # Starting values
  nb <- nbreg_fit(X, y, offset = offset, control = control$nb_control,
                  warn = FALSE)
  beta <- nb$coefficients
  theta <- nb$theta
  mu <- nb$fitted.values
  lz <- logit_irls(Z, as.numeric(y == 0), tol = control$logit_tol,
                   maxit = control$logit_maxit)
  gamma <- lz$gamma
  eta_zero <- lz$eta

  ll_old <- zinb_loglik(y, mu, theta, eta_zero)
  ll_trace <- ll_old
  converged <- FALSE

  for (iter in seq_len(control$maxit)) {
    # E-step
    tau <- zinb_estep(y, mu, theta, eta_zero)

    # M-step, zero part: maximize sum tau log pi + (1 - tau) log(1 - pi)
    lz <- logit_irls(Z, tau, gamma_start = gamma, tol = control$logit_tol,
                     maxit = control$logit_maxit)
    gamma <- lz$gamma
    eta_zero <- lz$eta

    # M-step, count part: maximize sum (1 - tau) log f_NB(y; mu, theta)
    nb <- nbreg_fit(X, y, offset = offset, prior_weights = 1 - tau,
                    control = control$nb_control,
                    beta_start = beta, theta_start = theta, warn = FALSE)
    beta <- nb$coefficients
    theta <- nb$theta
    mu <- nb$fitted.values

    # Observed-data log-likelihood
    ll <- zinb_loglik(y, mu, theta, eta_zero)
    ll_trace <- c(ll_trace, ll)
    if (control$trace)
      cat(sprintf("EM iter %4d: loglik = %.10f, theta = %.6f\n",
                  iter, ll, theta))
    if (ll < ll_old - 1e-8)
      warning(sprintf("EM log-likelihood decreased at iteration %d (%.10f -> %.10f).",
                      iter, ll_old, ll), call. = FALSE)

    if (abs(ll - ll_old) < control$tol) {
      converged <- TRUE
      break
    }
    ll_old <- ll
  }

  if (!converged)
    warning(sprintf("zinbreg() EM did not converge in %d iterations.",
                    control$maxit), call. = FALSE)
  if (theta > control$nb_control$theta_max)
    warning(sprintf(paste0("Estimated theta (%.3g) is very large: the count ",
                           "part is close to Poisson."), theta), call. = FALSE)

  # Standard errors from the observed information in (beta, gamma, log theta)
  vcov_full <- zinb_vcov(beta, gamma, theta, X, Z, y, offset)
  k <- nrow(vcov_full)
  se_logtheta <- sqrt(vcov_full[k, k])

  pi_ <- stats::plogis(eta_zero)
  list(coefficients = list(count = beta, zero = gamma),
       vcov = vcov_full[-k, -k, drop = FALSE], vcov_full = vcov_full,
       theta = theta, se_logtheta = se_logtheta,
       # delta method: SE(theta) = theta * SE(log theta)
       se_theta = theta * se_logtheta,
       loglik = ll, loglik_trace = ll_trace, iter = iter,
       converged = converged, fitted.values = (1 - pi_) * mu, mu = mu,
       pi = pi_, tau = zinb_estep(y, mu, theta, eta_zero),
       linear.predictors = list(count = drop(X %*% beta) + offset,
                                zero = eta_zero))
}
