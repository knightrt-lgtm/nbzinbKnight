# nbzinbKnight

<!-- badges: start -->
[![R-CMD-check](https://github.com/knightrt-lgtm/nbzinbKnight/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/knightrt-lgtm/nbzinbKnight/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**nbzinbKnight** fits two models for overdispersed count data, implementing
the estimation from scratch rather than with a general-purpose optimizer:

- **Negative binomial (NB2) regression**, `nbreg()`: iteratively reweighted
  least squares (IRLS) for the coefficients, alternated with Newton–Raphson
  for log(θ). Same model as `MASS::glm.nb()`.
- **Zero-inflated negative binomial (ZINB) regression**, `zinbreg()`: an EM
  algorithm whose M-steps are IRLS fits (a logistic regression for the
  zero-inflation part and a weighted NB fit for the count part). Same model as
  `pscl::zeroinfl(..., dist = "negbin")`.

Both functions reproduce the reference implementations on the `bioChemists`
and `NMES1988` data (see the comparison tests in `tests/testthat/`). This package was written for the STA 651 final project.

**Website:** <https://knightrt-lgtm.github.io/nbzinbKnight/>, with the
function reference and the
[project report](https://knightrt-lgtm.github.io/nbzinbKnight/articles/report.html).

## Installation

```r
# install.packages("devtools")
devtools::install_github("knightrt-lgtm/nbzinbKnight")
```

## Usage

```r
library(nbzinbKnight)
data("bioChemists", package = "pscl")

# Negative binomial regression
fit_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_nb)

# Zero-inflated negative binomial: count covariates | zero-model covariates
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                    data = bioChemists)
summary(fit_zinb)

# Model comparison and prediction
AIC(fit_nb, fit_zinb)
predict(fit_zinb, newdata = bioChemists[1:5, ], type = "zero")  # P(structural zero)
```

Both fits support `print`, `summary`, `coef`, `vcov`, `logLik`, `AIC`/`BIC`,
`fitted`, `residuals` (Pearson and response), `nobs`, and `predict`.
`predict()` for ZINB fits supports `type = "response"` ((1 − π)μ),
`"count"` (μ), and `"zero"` (π).

## R CMD check

`devtools::check()` finishes with 0 errors, 0 warnings, and 0 notes.
