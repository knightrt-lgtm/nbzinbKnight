# Negative binomial regression by IRLS and Newton-Raphson

Fits an NB2 regression model with log link, \\\log \mu_i = x_i^\top
\beta + o_i\\ and \\\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta\\, by
alternating iteratively reweighted least squares (IRLS) for \\\beta\\
(with \\\theta\\ fixed) and Newton-Raphson for \\\log\theta\\ (with
\\\beta\\ fixed), until the log-likelihood converges. This is the same
model and parameterization as
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html). No
general-purpose optimizer is used.

## Usage

``` r
nbreg(formula, data, offset, control = nbreg_control())
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`. Offsets may be given in the
  formula with [`offset()`](https://rdrr.io/r/stats/offset.html).

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset vector (e.g. log exposure), added to the linear
  predictor.

- control:

  A list from
  [`nbreg_control()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/nbreg_control.md).

## Value

An object of class `"nbreg"`, a list with components including
`coefficients`, `theta`, `loglik`, `iter`, `converged`, `fitted.values`,
`linear.predictors`, `weights` (final IRLS working weights), `y`, `x`,
`offset`, `terms`, and `call`.

## Examples

``` r
set.seed(1)
n <- 300
x <- rnorm(n)
y <- rnbinom(n, size = 2, mu = exp(0.5 + 0.8 * x))
fit <- nbreg(y ~ x, data = data.frame(y = y, x = x))
fit$coefficients
#> (Intercept)           x 
#>   0.4508058   0.7668655 
fit$theta
#> [1] 1.950002
```
