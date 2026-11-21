# Tests for the S3 methods of "nbreg" objects, checked against the
# corresponding methods for MASS::glm.nb fits.

fit_both <- function() {
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment
  list(d = bioChemists,
       ours = nbreg(f, data = bioChemists),
       ref = MASS::glm.nb(f, data = bioChemists,
                          control = stats::glm.control(epsilon = 1e-10,
                                                       maxit = 100)))
}

test_that("coef, vcov, nobs, fitted, residuals match glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  b <- fit_both()
  expect_equal(coef(b$ours), coef(b$ref), tolerance = 1e-6)
  expect_equal(vcov(b$ours), vcov(b$ref), tolerance = 1e-6)
  expect_equal(nobs(b$ours), nobs(b$ref))
  expect_equal(fitted(b$ours), fitted(b$ref), tolerance = 1e-6)
  expect_equal(residuals(b$ours, type = "response"),
               residuals(b$ref, type = "response"), tolerance = 1e-6)
  expect_equal(residuals(b$ours, type = "pearson"),
               residuals(b$ref, type = "pearson"), tolerance = 1e-6)
  expect_equal(residuals(b$ours), residuals(b$ours, type = "pearson"))
})

test_that("logLik has df so AIC and BIC work", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  b <- fit_both()
  ll <- logLik(b$ours)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), attr(logLik(b$ref), "df"))
  expect_equal(AIC(b$ours), AIC(b$ref), tolerance = 1e-8)
  expect_equal(BIC(b$ours), BIC(b$ref), tolerance = 1e-8)
})

test_that("summary gives the same coefficient table as glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  b <- fit_both()
  s <- summary(b$ours)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(unname(s$coefficients), unname(summary(b$ref)$coefficients),
               tolerance = 1e-6)
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_output(print(s), "Theta")
  expect_output(print(b$ours), "Coefficients")
})

test_that("predict works with and without newdata, including factors", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  b <- fit_both()
  expect_equal(predict(b$ours), fitted(b$ours))
  nd <- b$d[c(1, 50, 400, 900), ]
  expect_equal(predict(b$ours, nd),
               predict(b$ref, nd, type = "response"), tolerance = 1e-6)
  expect_equal(predict(b$ours, nd, type = "link"),
               predict(b$ref, nd, type = "link"), tolerance = 1e-6)
  one <- data.frame(fem = "Women", mar = "Married", kid5 = 0, phd = 3,
                    ment = 10)
  expect_equal(predict(b$ours, one), predict(b$ref, one, type = "response"),
               tolerance = 1e-6)
})

test_that("predict uses offsets from the formula and the argument", {
  set.seed(3)
  d <- data.frame(x = stats::rnorm(200), expo = stats::runif(200, 1, 4))
  d$y <- stats::rnbinom(200, size = 2, mu = d$expo * exp(0.2 + 0.5 * d$x))
  fit_f <- nbreg(y ~ x + offset(log(expo)), data = d)
  nd <- data.frame(x = 0, expo = 2)
  expect_equal(unname(predict(fit_f, nd, type = "link")),
               unname(coef(fit_f)[1]) + log(2))
  fit_a <- nbreg(y ~ x, data = d, offset = log(expo))
  expect_equal(unname(predict(fit_a, data.frame(x = 0), offset = log(2))),
               exp(unname(coef(fit_a)[1]) + log(2)))
})
