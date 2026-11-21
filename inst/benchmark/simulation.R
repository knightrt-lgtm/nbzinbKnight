# Simulation study: does zinbreg() recover known parameters?
# (project description, Section 5.1)
#
# Design: R = 200 data sets, n = 500 each, from a ZINB model with one
# covariate in each part:
#   count model: log mu_i  = beta0 + beta1 * x_i,      x_i ~ N(0, 1)
#   zero model:  logit pi_i = gamma0 + gamma1 * z_i,   z_i ~ N(0, 1)
#   Y_i = 0 with probability pi_i, otherwise Y_i ~ NB(mu_i, theta)
# True values: beta = (1, 0.5), gamma = (-0.5, 1), theta = 2.
#
# For each parameter (and log theta) we report the bias of the estimate, the
# mean of the estimated SEs vs the Monte Carlo SD of the estimates, and the
# empirical coverage of 95% Wald intervals estimate +/- 1.96 * SE.
#
# Run from the package root after installing the package:
#   Rscript inst/benchmark/simulation.R
# Results are saved to inst/benchmark/simulation_results.rds for the report.

library(nbzinbKnight)
set.seed(20261121)

R <- 200
n <- 500
truth <- c(`count_(Intercept)` = 1, count_x = 0.5,
           `zero_(Intercept)` = -0.5, zero_z = 1, `Log(theta)` = log(2))
theta_true <- exp(truth[["Log(theta)"]])

simulate_zinb <- function(n) {
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  mu <- exp(truth[["count_(Intercept)"]] + truth[["count_x"]] * x)
  pi_ <- stats::plogis(truth[["zero_(Intercept)"]] + truth[["zero_z"]] * z)
  structural_zero <- stats::rbinom(n, 1, pi_)
  y <- ifelse(structural_zero == 1, 0,
              stats::rnbinom(n, size = theta_true, mu = mu))
  data.frame(y = y, x = x, z = z)
}

est <- se <- matrix(NA_real_, R, length(truth),
                    dimnames = list(NULL, names(truth)))
converged <- logical(R)
iters <- integer(R)
zero_share <- numeric(R)

t0 <- proc.time()[["elapsed"]]
for (r in seq_len(R)) {
  d <- simulate_zinb(n)
  zero_share[r] <- mean(d$y == 0)
  fit <- tryCatch(suppressWarnings(zinbreg(y ~ x | z, data = d)),
                  error = function(e) NULL)
  if (is.null(fit)) next
  converged[r] <- fit$converged
  iters[r] <- fit$iter
  est[r, ] <- c(coef(fit), log(fit$theta))
  se[r, ] <- c(sqrt(diag(vcov(fit))), fit$se_logtheta)
}
elapsed <- proc.time()[["elapsed"]] - t0

# Use replications that converged with finite SEs
ok <- converged & stats::complete.cases(est) & stats::complete.cases(se)
E <- est[ok, , drop = FALSE]
S <- se[ok, , drop = FALSE]
lower <- E - 1.96 * S
upper <- E + 1.96 * S
covered <- sweep(lower, 2, truth, "<=") & sweep(upper, 2, truth, ">=")

summary_table <- data.frame(
  parameter = names(truth),
  true = unname(truth),
  mean_estimate = colMeans(E),
  bias = colMeans(E) - truth,
  mc_sd = apply(E, 2, stats::sd),
  mean_se = colMeans(S),
  coverage_95 = colMeans(covered),
  row.names = NULL
)

results <- list(
  design = list(R = R, n = n, truth = truth, seed = 20261121),
  summary = summary_table,
  n_used = sum(ok), n_converged = sum(converged),
  iterations = summary(iters[ok]),
  mean_zero_share = mean(zero_share),
  estimates = est, std_errors = se, converged = converged,
  elapsed_seconds = elapsed
)
saveRDS(results, file.path("inst", "benchmark", "simulation_results.rds"))

cat(sprintf("Replications used: %d of %d (converged: %d); mean share of zeros: %.3f\n",
            sum(ok), R, sum(converged), mean(zero_share)))
cat(sprintf("EM iterations: median %d, range %d-%d; total time %.0f s\n",
            as.integer(stats::median(iters[ok])), min(iters[ok]),
            max(iters[ok]), elapsed))
# Monte Carlo SE of a coverage proportion near 0.95 with this many replications
cat(sprintf("Monte Carlo SE of coverage near 0.95: %.3f\n",
            sqrt(0.95 * 0.05 / sum(ok))))
print(summary_table, digits = 3, row.names = FALSE)
