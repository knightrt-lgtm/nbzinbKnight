# Control parameters for `zinbreg()`

Control parameters for
[`zinbreg()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/zinbreg.md)

## Usage

``` r
zinbreg_control(
  tol = 1e-08,
  maxit = 5000,
  nb_control = nbreg_control(),
  logit_tol = 1e-10,
  logit_maxit = 100,
  trace = FALSE
)
```

## Arguments

- tol:

  EM convergence tolerance: stop when the observed-data log-likelihood
  changes by less than `tol` between EM iterations.

- maxit:

  Maximum number of EM iterations.

- nb_control:

  A list from
  [`nbreg_control()`](https://knightrt-lgtm.github.io/nbzinbKnight/reference/nbreg_control.md)
  used for the count-part M-step (the weighted NB fit).

- logit_tol, logit_maxit:

  Tolerance and maximum iterations for the zero-part M-step (the IRLS
  logistic regression).

- trace:

  Logical; if `TRUE`, print the observed-data log-likelihood at each EM
  iteration.

## Value

A list of control parameters.

## Examples

``` r
zinbreg_control(tol = 1e-10, maxit = 1000)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 1000
#> 
#> $nb_control
#> $nb_control$tol
#> [1] 1e-08
#> 
#> $nb_control$maxit
#> [1] 100
#> 
#> $nb_control$irls_tol
#> [1] 1e-10
#> 
#> $nb_control$irls_maxit
#> [1] 100
#> 
#> $nb_control$theta_tol
#> [1] 1e-10
#> 
#> $nb_control$theta_maxit
#> [1] 100
#> 
#> $nb_control$theta_max
#> [1] 1e+06
#> 
#> $nb_control$trace
#> [1] FALSE
#> 
#> 
#> $logit_tol
#> [1] 1e-10
#> 
#> $logit_maxit
#> [1] 100
#> 
#> $trace
#> [1] FALSE
#> 
```
