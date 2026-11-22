# Methods for `zinbreg` fits

Standard extractor and summary methods for objects returned by
[`zinbreg()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/zinbreg.md).
Coefficient vectors and covariance matrices use the labels `count_*` and
`zero_*`, as
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) does.

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
predict(
  object,
  newdata,
  type = c("response", "count", "zero"),
  offset = NULL,
  ...
)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- object, x:

  An object of class `"zinbreg"` (or `"summary.zinbreg"` for the print
  method of a summary).

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (default,
  count and zero parts together), `"count"`, or `"zero"`.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"pearson"` (default), \\(y_i - \hat E
  Y_i)/\sqrt{\widehat{\mathrm{Var}}(Y_i)}\\, or `"response"`, \\y_i -
  \hat E Y_i\\, where \\E Y_i = (1-\pi_i)\mu_i\\ and \\\mathrm{Var}(Y_i)
  = (1-\pi_i)\mu_i\\1 + \mu_i(\pi_i + 1/\theta)\\\\. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (default) for the mean \\(1-\hat\pi)\hat\mu\\, `"count"` for the
  count-model mean \\\hat\mu\\, or `"zero"` for the zero-inflation
  probability \\\hat\pi\\.

- newdata:

  Optional data frame for
  [`predict()`](https://rdrr.io/r/stats/predict.html). If omitted,
  predictions are for the data used in the fit.

- offset:

  For [`predict()`](https://rdrr.io/r/stats/predict.html) with
  `newdata`: an optional count-model offset vector (in addition to any
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the
  formula).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named coefficient
vector; [`vcov()`](https://rdrr.io/r/stats/vcov.html) its covariance
matrix; [`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of
class `"logLik"` with `df` = number of count and zero coefficients + 1
(for \\\theta\\) and `nobs`;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) the fitted
means \\(1-\hat\pi)\hat\mu\\;
[`residuals()`](https://rdrr.io/r/stats/residuals.html) a residual
vector; [`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
observations; [`predict()`](https://rdrr.io/r/stats/predict.html) a
vector of predictions;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`"summary.zinbreg"` with separate coefficient tables for the count and
zero parts.

## Examples

``` r
set.seed(1)
n <- 500
d <- data.frame(x = rnorm(n))
d$y <- ifelse(rbinom(n, 1, plogis(-1 + 0.5 * d$x)) == 1, 0,
              rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
fit <- zinbreg(y ~ x | x, data = d)
fit
#> 
#> Call:
#> zinbreg(formula = y ~ x | x, data = d)
#> 
#> Zero-inflated negative binomial regression, fitted by EM with IRLS M-steps
#> 
#> Count model coefficients (negbin with log link):
#> (Intercept)            x  
#>      1.1015       0.4487  
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#> (Intercept)            x  
#>     -0.8284       0.5318  
#> 
#> Theta: 2.322   Log-likelihood: -927.076 on 5 df
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | x, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -0.8304 -0.7848 -0.3948  0.5027  4.9886 
#> 
#> Count model coefficients (negbin with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.10147    0.06142  17.934  < 2e-16 ***
#> x            0.44871    0.05427   8.268  < 2e-16 ***
#> Log(theta)   0.84258    0.22107   3.811 0.000138 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.8284     0.1807  -4.584 4.55e-06 ***
#> x             0.5318     0.1529   3.477 0.000507 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.322  (SE 0.5134)
#> Log-likelihood: -927.076 on 5 df   AIC: 1864.15   BIC: 1885.23
#> EM iterations: 43
coef(fit, model = "zero")
#> (Intercept)           x 
#>  -0.8283902   0.5317871 
AIC(fit)
#> [1] 1864.153
predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "zero")
#>         1         2         3 
#> 0.2042115 0.3039856 0.4263881 
```
