# Tests for nbreg(). Reference fits use MASS::glm.nb() with tight tolerances.
# Required agreement (project description, Table 4, NB column):
#   coefficients <= 1e-5, theta <= 1e-4, log-likelihood <= 1e-4.
#   standard errors <= 1e-5 (SE(theta) checked to 1e-4).

ref_control <- function() stats::glm.control(epsilon = 1e-10, maxit = 100)

expect_matches_glm_nb <- function(ours, ref) {
  expect_true(ours$converged)
  expect_equal(names(ours$coefficients), names(stats::coef(ref)))
  expect_lt(max(abs(ours$coefficients - stats::coef(ref))), 1e-5)
  expect_lt(abs(ours$theta - ref$theta), 1e-4)
  expect_lt(max(abs(sqrt(diag(ours$vcov)) - sqrt(diag(stats::vcov(ref))))), 1e-5)
  expect_lt(abs(ours$se_theta - ref$SE.theta), 1e-4)
  expect_lt(abs(ours$loglik - as.numeric(stats::logLik(ref))), 1e-4)
}

test_that("nbreg matches glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment
  ours <- nbreg(f, data = bioChemists)
  ref <- MASS::glm.nb(f, data = bioChemists, control = ref_control())
  expect_matches_glm_nb(ours, ref)
})

test_that("nbreg matches glm.nb on NMES1988", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  ours <- nbreg(f, data = NMES1988)
  ref <- MASS::glm.nb(f, data = NMES1988, control = ref_control())
  expect_matches_glm_nb(ours, ref)
})

test_that("nbreg works with a single covariate", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  ours <- nbreg(art ~ ment, data = bioChemists)
  ref <- MASS::glm.nb(art ~ ment, data = bioChemists, control = ref_control())
  expect_length(ours$coefficients, 2)
  expect_matches_glm_nb(ours, ref)
})

test_that("offset argument and offset() in the formula agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  set.seed(2)
  bioChemists$expo <- stats::runif(nrow(bioChemists), 0.5, 3)
  f <- art ~ fem + mar + kid5 + phd + ment
  f_off <- art ~ fem + mar + kid5 + phd + ment + offset(log(expo))

  ours_arg <- nbreg(f, data = bioChemists, offset = log(expo))
  ours_formula <- nbreg(f_off, data = bioChemists)
  ref <- MASS::glm.nb(f_off, data = bioChemists, control = ref_control())

  expect_equal(ours_arg$coefficients, ours_formula$coefficients)
  expect_equal(ours_arg$offset, log(bioChemists$expo))
  expect_matches_glm_nb(ours_arg, ref)
})

test_that("the log-likelihood never decreases across outer iterations", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  out <- utils::capture.output(
    fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
                 control = nbreg_control(trace = TRUE)))
  out <- grep("^Outer iter", out, value = TRUE)
  ll <- as.numeric(sub(".*loglik = ([-0-9.]+),.*", "\\1", out))
  expect_equal(length(ll), fit$iter)
  expect_true(all(diff(ll) >= -1e-10))
})

test_that("bad input gives informative errors", {
  d <- data.frame(y = c(0, 1, 2, 3, 4), x = 1:5)
  expect_error(nbreg(y ~ x, data = transform(d, y = c(-1, 1, 2, 3, 4))),
               "non-negative integer")
  expect_error(nbreg(y ~ x, data = transform(d, y = c(0.5, 1, 2, 3, 4))),
               "non-negative integer")
  expect_error(nbreg(~ x, data = d), "response")
  expect_error(nbreg_control(tol = -1), "tol")
  expect_error(nbreg_control(maxit = 0), "maxit")
})

test_that("nbreg warns when it hits maxit without converging", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  expect_warning(
    fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
                 control = nbreg_control(maxit = 1)),
    "did not converge")
  expect_false(fit$converged)
})

test_that("nbreg warns when theta is very large (Poisson data)", {
  set.seed(10)
  n <- 500
  x <- stats::rnorm(n)
  # Underdispersed counts: no overdispersion at all, so theta -> infinity
  y <- stats::rbinom(n, size = 10, prob = stats::plogis(-1 + 0.3 * x))
  expect_warning(nbreg(y ~ x, data = data.frame(y = y, x = x)),
                 "close to Poisson")
})
