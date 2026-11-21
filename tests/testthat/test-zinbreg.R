# Tests for zinbreg(). Reference fits use pscl::zeroinfl(dist = "negbin")
# with a tight tolerance. Required agreement (Table 4, ZINB column):
#   coefficients <= 1e-3, log(theta) <= 1e-3, log-likelihood <= 1e-3.
# (Standard errors are added and tested in a later checkpoint.)

ref_zinb <- function(formula, data) {
  pscl::zeroinfl(formula, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12))
}

expect_matches_zeroinfl <- function(ours, ref) {
  expect_true(ours$converged)
  expect_equal(names(ours$coefficients$count), names(ref$coefficients$count))
  expect_equal(names(ours$coefficients$zero), names(ref$coefficients$zero))
  expect_lt(max(abs(ours$coefficients$count - ref$coefficients$count)), 1e-3)
  expect_lt(max(abs(ours$coefficients$zero - ref$coefficients$zero)), 1e-3)
  expect_lt(abs(log(ours$theta) - log(ref$theta)), 1e-3)
  expect_lt(abs(ours$loglik - as.numeric(stats::logLik(ref))), 1e-3)
}

test_that("zinbreg matches zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ours <- zinbreg(f, data = bioChemists)
  expect_matches_zeroinfl(ours, ref_zinb(f, bioChemists))
})

test_that("zinbreg matches zeroinfl on NMES1988", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  ours <- zinbreg(f, data = NMES1988)
  expect_matches_zeroinfl(ours, ref_zinb(f, NMES1988))
})

test_that("the EM log-likelihood never decreases", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  ours <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  expect_equal(length(ours$loglik_trace), ours$iter + 1)
  expect_true(all(diff(ours$loglik_trace) >= -1e-10))
})

test_that("a formula without | uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  a <- zinbreg(art ~ fem + kid5 + ment, data = bioChemists)
  b <- zinbreg(art ~ fem + kid5 + ment | fem + kid5 + ment, data = bioChemists)
  expect_equal(a$coefficients, b$coefficients)
  expect_equal(colnames(a$z), colnames(a$x))
})

test_that("different covariates per part, and offsets, agree with zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  set.seed(2)
  bioChemists$expo <- stats::runif(nrow(bioChemists), 0.5, 3)
  f <- art ~ fem + kid5 + ment | ment

  ours_arg <- zinbreg(f, data = bioChemists, offset = log(expo))
  ours_formula <- zinbreg(art ~ fem + kid5 + ment + offset(log(expo)) | ment,
                          data = bioChemists)
  # zeroinfl evaluates `offset` lazily, so call it directly (not via ref_zinb)
  ref <- pscl::zeroinfl(f, data = bioChemists, dist = "negbin",
                        offset = log(expo),
                        control = pscl::zeroinfl.control(reltol = 1e-12))

  expect_length(ours_arg$coefficients$zero, 2)
  expect_equal(ours_arg$coefficients, ours_formula$coefficients)
  expect_matches_zeroinfl(ours_arg, ref)
})

test_that("fitted values, pi, and tau are consistent", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
  expect_equal(fit$fitted.values, (1 - fit$pi) * fit$mu)
  expect_true(all(fit$pi > 0 & fit$pi < 1))
  expect_true(all(fit$tau[bioChemists$art > 0] == 0))
  expect_true(all(fit$tau >= 0 & fit$tau <= 1))
})

test_that("the log-likelihood function matches zeroinfl at its estimates", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + kid5 + ment | ment
  ref <- ref_zinb(f, bioChemists)
  X <- stats::model.matrix(~ fem + kid5 + ment, bioChemists)
  Z <- stats::model.matrix(~ ment, bioChemists)
  mu <- exp(drop(X %*% ref$coefficients$count))
  eta_zero <- drop(Z %*% ref$coefficients$zero)
  expect_equal(zinb_loglik(bioChemists$art, mu, ref$theta, eta_zero),
               as.numeric(stats::logLik(ref)), tolerance = 1e-10)
})

test_that("bad input gives informative errors", {
  d <- data.frame(y = c(0, 1, 2, 3, 0), x = 1:5)
  expect_error(zinbreg(y ~ x, data = transform(d, y = c(-1, 1, 2, 3, 0))),
               "non-negative integer")
  expect_error(zinbreg(y ~ x, data = transform(d, y = c(0.5, 1, 2, 3, 0))),
               "non-negative integer")
  expect_error(zinbreg(y ~ x, data = transform(d, y = c(1, 1, 2, 3, 4))),
               "no zeros")
  expect_error(zinbreg_control(tol = 0), "tol")
  expect_error(zinbreg_control(maxit = 0), "maxit")
})

test_that("zinbreg warns when EM hits maxit without converging", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  expect_warning(
    fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists,
                   control = zinbreg_control(maxit = 2)),
    "did not converge")
  expect_false(fit$converged)
})
