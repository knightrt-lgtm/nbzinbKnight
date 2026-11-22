# Zero-inflated negative binomial regression by EM with IRLS M-steps

Fits a zero-inflated negative binomial (ZINB) model. With probability
\\\pi_i\\ an observation is a structural zero; otherwise it comes from
an NB2 distribution with mean \\\mu_i\\ and dispersion \\\theta\\. The
count model is \\\log \mu_i = x_i^\top\beta + o_i\\ and the zero model
is \\\mathrm{logit}\\\pi_i = z_i^\top\gamma\\. The observed-data
log-likelihood is maximized by an EM algorithm: the E-step computes the
posterior probability \\\tau_i\\ that each zero is structural, and the
M-steps are a logistic IRLS fit for \\\gamma\\ (fractional response
\\\tau_i\\) and a weighted NB fit for \\(\beta, \theta)\\ with prior
weights \\1 - \tau_i\\. No general-purpose optimizer is used.

## Usage

``` r
zinbreg(formula, data, offset, control = zinbreg_control())
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before `|`, zero-model covariates after it. Without `|`, the same
  covariates are used in both parts.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset for the count model (e.g. log exposure). An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the count
  part of the formula also works.

- control:

  A list from
  [`zinbreg_control()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/zinbreg_control.md).

## Value

An object of class `"zinbreg"`, a list with components including
`coefficients` (a list with `count` and `zero`), `theta`, `loglik`,
`loglik_trace` (log-likelihood at every EM iteration), `iter`,
`converged`, `fitted.values` (\\(1-\pi_i)\mu_i\\), `mu`, `pi`, `tau`,
`y`, `x`, `z`, `offset`, `terms`, `formula`, and `call`.

## Examples

``` r
set.seed(1)
n <- 500
x <- rnorm(n)
structural_zero <- rbinom(n, 1, plogis(-1 + 0.5 * x))
y <- ifelse(structural_zero == 1, 0,
            rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)))
fit <- zinbreg(y ~ x | x, data = data.frame(y = y, x = x))
fit$coefficients
#> $count
#> (Intercept)           x 
#>   1.1014731   0.4487059 
#> 
#> $zero
#> (Intercept)           x 
#>  -0.8283902   0.5317871 
#> 
fit$theta
#> [1] 2.322346
```
