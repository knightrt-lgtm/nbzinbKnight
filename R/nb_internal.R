# Internal building blocks for NB2 regression.
# Equation references are to project-docs/STA651_Final_Project.pdf.

# NB2 log-likelihood (Section 2.1), with optional prior weights a_i:
#   l(beta, theta) = sum_i a_i [ log Gamma(y_i + theta) - log Gamma(theta)
#                                - log y_i! + theta log theta + y_i log mu_i
#                                - (theta + y_i) log(theta + mu_i) ]
# log y_i! is computed as lgamma(y_i + 1).
#' @noRd
nb_loglik <- function(y, mu, theta, prior_weights = rep(1, length(y))) {
  sum(prior_weights * (lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
                         theta * log(theta) + y * log(mu) -
                         (theta + y) * log(theta + mu)))
}

# IRLS (Fisher scoring) for beta with theta held fixed (Section 3.1(a)).
#
# Log link, variance V(mu) = mu + mu^2 / theta. Each iteration:
#   working weights  w_i = a_i * theta * mu_i / (theta + mu_i)
#   working response z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
#   update           beta = (X' W X)^{-1} X' W z
# The weighted least-squares step is solved with a QR decomposition of
# sqrt(W) X, which gives the same beta as (X'WX)^{-1} X'Wz but is more
# numerically stable than forming X'WX directly.
#
# Arguments:
#   X             n x p model matrix
#   y             response counts
#   theta         fixed dispersion parameter (> 0)
#   offset        o_i (default 0)
#   prior_weights a_i (default 1; ZINB will pass 1 - tau_i)
#   beta_start    starting beta; if NULL, start from mu_i = (y_i + ybar) / 2
#   tol, maxit    stop when max |beta_new - beta_old| < tol
#
# Returns a list with beta, eta, mu, the final working weights w,
# the number of iterations, and a convergence flag.
#' @noRd
nb_irls_beta <- function(X, y, theta,
                         offset = rep(0, length(y)),
                         prior_weights = rep(1, length(y)),
                         beta_start = NULL,
                         tol = 1e-10, maxit = 100) {
  if (is.null(beta_start)) {
    # No beta yet: start from fitted means pulled halfway toward the weighted
    # mean of y (positive even when y_i = 0, so log(mu) is finite).
    ybar <- sum(prior_weights * y) / sum(prior_weights)
    mu <- (y + ybar) / 2
    eta <- log(mu)
    beta <- NULL
  } else {
    beta <- beta_start
    eta <- drop(X %*% beta) + offset
    mu <- exp(eta)
  }

  converged <- FALSE
  for (iter in seq_len(maxit)) {
    # Working weights and working response (Section 3.1(a))
    w <- prior_weights * theta * mu / (theta + mu)
    z <- (eta - offset) + (y - mu) / mu

    # Weighted least squares: beta = (X'WX)^{-1} X'Wz, via QR of sqrt(W) X
    sw <- sqrt(w)
    beta_new <- qr.coef(qr(sw * X), sw * z)

    # Update linear predictor and mean
    eta <- drop(X %*% beta_new) + offset
    mu <- exp(eta)

    # Convergence: change in beta below tolerance
    if (!is.null(beta) && max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }

  # Weights at the final beta (used later for Var(beta) = (X'WX)^{-1})
  w <- prior_weights * theta * mu / (theta + mu)
  names(beta) <- colnames(X)

  list(beta = beta, eta = eta, mu = mu, weights = w,
       iter = iter, converged = converged)
}
