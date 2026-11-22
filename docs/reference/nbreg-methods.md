# Methods for `nbreg` fits

Standard extractor and summary methods for objects returned by
[`nbreg()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("response", "link"), offset = NULL, ...)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- object, x:

  An object of class `"nbreg"` (or `"summary.nbreg"` for the print
  method of a summary).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"pearson"` (default), \\(y_i - \hat\mu_i)/\sqrt{\hat\mu_i +
  \hat\mu_i^2/\hat\theta}\\, or `"response"`, \\y_i - \hat\mu_i\\. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (default) for \\\hat\mu\\, or `"link"` for the linear predictor
  \\\hat\eta = \log\hat\mu\\.

- newdata:

  Optional data frame for
  [`predict()`](https://rdrr.io/r/stats/predict.html). If omitted,
  predictions are for the data used in the fit.

- offset:

  For [`predict()`](https://rdrr.io/r/stats/predict.html) with
  `newdata`: an optional offset vector (in addition to any
  [`offset()`](https://rdrr.io/r/stats/offset.html) term in the
  formula).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) the coefficient vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) its covariance matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`"logLik"` with `df` = number of coefficients + 1 (for \\\theta\\) and
`nobs`; [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) the
fitted means; [`residuals()`](https://rdrr.io/r/stats/residuals.html) a
residual vector; [`nobs()`](https://rdrr.io/r/stats/nobs.html) the
number of observations;
[`predict()`](https://rdrr.io/r/stats/predict.html) a vector of
predictions; [`summary()`](https://rdrr.io/r/base/summary.html) an
object of class `"summary.nbreg"` containing the coefficient table.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(300))
d$y <- rnbinom(300, size = 2, mu = exp(0.5 + 0.8 * d$x))
fit <- nbreg(y ~ x, data = d)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Negative binomial (NB2) regression, log link, fitted by IRLS + Newton-Raphson
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.4508       0.7669  
#> 
#> Theta: 1.95   Log-likelihood: -526.892 on 3 df
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1939 -0.7791 -0.2593  0.4389  5.2098 
#> 
#> Coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.45081    0.06613   6.817 9.32e-12 ***
#> x            0.76687    0.06699  11.448  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 1.95  (SE 0.3699)
#> Log-likelihood: -526.892 on 3 df   AIC: 1059.78   BIC: 1070.89
#> Outer iterations (IRLS + Newton-Raphson): 3
coef(fit)
#> (Intercept)           x 
#>   0.4508058   0.7668655 
AIC(fit)
#> [1] 1059.783
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -0.8050751  2.2472901 -0.7620322  1.2691926 -0.5032198  0.1494196 
predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>        1        2        3 
#> 0.729016 1.569576 3.379309 
```
