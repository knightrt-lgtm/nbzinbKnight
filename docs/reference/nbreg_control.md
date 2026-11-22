# Control parameters for `nbreg()`

Control parameters for
[`nbreg()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/nbreg.md)

## Usage

``` r
nbreg_control(
  tol = 1e-08,
  maxit = 100,
  irls_tol = 1e-10,
  irls_maxit = 100,
  theta_tol = 1e-10,
  theta_maxit = 100,
  theta_max = 1e+06,
  trace = FALSE
)
```

## Arguments

- tol:

  Convergence tolerance for the outer loop: stop when the log-likelihood
  changes by less than `tol` between iterations.

- maxit:

  Maximum number of outer iterations (alternating the IRLS update for
  beta and the Newton-Raphson update for theta).

- irls_tol, irls_maxit:

  Tolerance and maximum iterations for each inner IRLS fit of beta
  (theta fixed).

- theta_tol, theta_maxit:

  Tolerance and maximum iterations for each inner Newton-Raphson fit of
  log(theta) (beta fixed).

- theta_max:

  Warn if the estimated theta exceeds this value (the model is then
  close to Poisson).

- trace:

  Logical; if `TRUE`, print the log-likelihood and theta at each outer
  iteration.

## Value

A list of control parameters.

## Examples

``` r
nbreg_control(tol = 1e-10, trace = TRUE)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 100
#> 
#> $irls_tol
#> [1] 1e-10
#> 
#> $irls_maxit
#> [1] 100
#> 
#> $theta_tol
#> [1] 1e-10
#> 
#> $theta_maxit
#> [1] 100
#> 
#> $theta_max
#> [1] 1e+06
#> 
#> $trace
#> [1] TRUE
#> 
```
