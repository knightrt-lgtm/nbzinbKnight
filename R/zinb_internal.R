# Internal building blocks for ZINB regression.
# Equation references are to project-docs/STA651_Final_Project.pdf.

# log f_NB(0; mu, theta) = theta * [log theta - log(theta + mu)]
# (Section 2.2: f_NB(0; mu, theta) = {theta / (theta + mu)}^theta)
#' @noRd
nb_log_f0 <- function(mu, theta) {
  theta * (log(theta) - log(theta + mu))
}

# log(exp(a) + exp(b)) computed without overflow/underflow
#' @noRd
log_sum_exp <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

# Observed-data ZINB log-likelihood (Section 2.2):
#   l = sum_{y_i = 0} log[ pi_i + (1 - pi_i) f_NB(0; mu_i, theta) ]
#     + sum_{y_i > 0} [ log(1 - pi_i) + log f_NB(y_i; mu_i, theta) ]
# with logit(pi_i) = eta_zero_i = z_i' gamma.
# log(pi) and log(1 - pi) come from plogis(..., log.p = TRUE) so they stay
# accurate when pi is very close to 0 or 1; the y = 0 term uses log_sum_exp.
#' @noRd
zinb_loglik <- function(y, mu, theta, eta_zero) {
  log_pi <- stats::plogis(eta_zero, log.p = TRUE)       # log pi_i
  log_1mpi <- stats::plogis(-eta_zero, log.p = TRUE)    # log(1 - pi_i)
  zero <- y == 0
  ll0 <- log_sum_exp(log_pi[zero], log_1mpi[zero] + nb_log_f0(mu[zero], theta))
  ll1 <- log_1mpi[!zero] +
    stats::dnbinom(y[!zero], size = theta, mu = mu[!zero], log = TRUE)
  sum(ll0) + sum(ll1)
}

# E-step (Section 3.2): posterior probability of being a structural zero
#   tau_i = pi_i / [pi_i + (1 - pi_i) f_NB(0; mu_i, theta)]   if y_i = 0
#   tau_i = 0                                                 if y_i > 0
# Computed on the log scale: tau_i = exp(log pi_i - log[denominator]).
#' @noRd
zinb_estep <- function(y, mu, theta, eta_zero) {
  tau <- numeric(length(y))
  zero <- y == 0
  log_pi <- stats::plogis(eta_zero[zero], log.p = TRUE)
  log_1mpi <- stats::plogis(-eta_zero[zero], log.p = TRUE)
  log_denom <- log_sum_exp(log_pi, log_1mpi + nb_log_f0(mu[zero], theta))
  tau[zero] <- exp(log_pi - log_denom)
  tau
}

# M-step for the zero part (Section 3.2): logistic regression of the
# fractional response tau_i on z_i, maximizing
#   sum_i [ tau_i log pi_i + (1 - tau_i) log(1 - pi_i) ]
# by IRLS with
#   working weights  w_i = pi_i (1 - pi_i)
#   working response z_i = z_i' gamma + (tau_i - pi_i) / [pi_i (1 - pi_i)]
#   update           gamma = (Z' W Z)^{-1} Z' W z   (QR of sqrt(W) Z)
#
# If gamma_start is NULL, start from pi_i = (tau_i + 0.5) / 2, which is
# strictly between 0 and 1 (the same idea glm() uses for binomial data).
# Returns gamma, eta (= Z gamma), pi, final weights, iterations, convergence.
#' @noRd
logit_irls <- function(Z, tau, gamma_start = NULL, tol = 1e-10, maxit = 100) {
  if (is.null(gamma_start)) {
    p <- (tau + 0.5) / 2
    eta <- stats::qlogis(p)
    gamma <- NULL
  } else {
    gamma <- gamma_start
    eta <- drop(Z %*% gamma)
    p <- stats::plogis(eta)
  }

  converged <- FALSE
  for (iter in seq_len(maxit)) {
    # Working weights and working response
    w <- p * (1 - p)
    z <- eta + (tau - p) / w

    # Weighted least squares via QR of sqrt(W) Z
    sw <- sqrt(w)
    gamma_new <- qr.coef(qr(sw * Z), sw * z)

    eta <- drop(Z %*% gamma_new)
    p <- stats::plogis(eta)

    if (!is.null(gamma) && max(abs(gamma_new - gamma)) < tol) {
      gamma <- gamma_new
      converged <- TRUE
      break
    }
    gamma <- gamma_new
  }

  names(gamma) <- colnames(Z)
  list(gamma = gamma, eta = eta, pi = p, weights = p * (1 - p),
       iter = iter, converged = converged)
}
