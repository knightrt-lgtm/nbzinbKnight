# Tests for the S3 methods of "zinbreg" objects, checked against the
# corresponding methods for pscl::zeroinfl fits. EM stops about 1e-4 from
# zeroinfl's optimum, so comparisons use the Table 4 ZINB tolerance (1e-3).

zinb_fit_both <- function() {
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment | fem + kid5 + ment
  list(d = bioChemists,
       ours = zinbreg(f, data = bioChemists),
       ref = pscl::zeroinfl(f, data = bioChemists, dist = "negbin",
                            control = pscl::zeroinfl.control(reltol = 1e-12)))
}

max_diff <- function(a, b) max(abs(a - b))

test_that("coef and vcov match zeroinfl, including labels", {
  skip_if_not_installed("pscl")
  b <- zinb_fit_both()
  expect_equal(names(coef(b$ours)), names(coef(b$ref)))
  expect_lt(max_diff(coef(b$ours), coef(b$ref)), 1e-3)
  expect_lt(max_diff(coef(b$ours, "count"), coef(b$ref, "count")), 1e-3)
  expect_lt(max_diff(coef(b$ours, "zero"), coef(b$ref, "zero")), 1e-3)
  expect_equal(dimnames(vcov(b$ours)), dimnames(vcov(b$ref)))
  expect_lt(max_diff(vcov(b$ours), vcov(b$ref)), 1e-3)
  expect_equal(dimnames(vcov(b$ours, "zero"))[[1]], names(coef(b$ours, "zero")))
})

test_that("logLik has df so AIC and BIC work", {
  skip_if_not_installed("pscl")
  b <- zinb_fit_both()
  expect_equal(attr(logLik(b$ours), "df"), attr(logLik(b$ref), "df"))
  expect_lt(abs(AIC(b$ours) - AIC(b$ref)), 1e-3)
  expect_lt(abs(BIC(b$ours) - BIC(b$ref)), 1e-3)
  expect_equal(nobs(b$ours), b$ref$n)
})

test_that("fitted and residuals match zeroinfl", {
  skip_if_not_installed("pscl")
  b <- zinb_fit_both()
  expect_lt(max_diff(fitted(b$ours), fitted(b$ref)), 1e-3)
  expect_lt(max_diff(residuals(b$ours, "response"),
                     residuals(b$ref, "response")), 1e-3)
  expect_lt(max_diff(residuals(b$ours, "pearson"),
                     residuals(b$ref, "pearson")), 1e-3)
  expect_equal(residuals(b$ours), residuals(b$ours, "pearson"))
})

test_that("summary has count (with Log(theta)) and zero tables", {
  skip_if_not_installed("pscl")
  b <- zinb_fit_both()
  s <- summary(b$ours)
  sr <- summary(b$ref)$coefficients
  expect_s3_class(s, "summary.zinbreg")
  expect_equal(dimnames(s$coefficients$count), dimnames(sr$count))
  expect_equal(dimnames(s$coefficients$zero), dimnames(sr$zero))
  expect_lt(max_diff(s$coefficients$count[, 1:2], sr$count[, 1:2]), 1e-3)
  expect_lt(max_diff(s$coefficients$zero[, 1:2], sr$zero[, 1:2]), 1e-3)
  expect_output(print(s), "Zero-inflation model coefficients")
  expect_output(print(b$ours), "Count model coefficients")
})

test_that("predict supports response, count, and zero, with newdata", {
  skip_if_not_installed("pscl")
  b <- zinb_fit_both()
  nd <- b$d[c(1, 50, 400, 900), ]
  for (ty in c("response", "count", "zero")) {
    expect_lt(max_diff(predict(b$ours, nd, type = ty),
                       predict(b$ref, nd, type = ty)), 1e-3)
    expect_lt(max_diff(predict(b$ours, type = ty),
                       predict(b$ref, type = ty)), 1e-3)
  }
  expect_equal(predict(b$ours, nd, type = "response"),
               (1 - predict(b$ours, nd, type = "zero")) *
                 predict(b$ours, nd, type = "count"))
  one <- data.frame(fem = "Women", mar = "Married", kid5 = 0, phd = 3,
                    ment = 10)
  expect_lt(max_diff(predict(b$ours, one), predict(b$ref, one)), 1e-3)
})
