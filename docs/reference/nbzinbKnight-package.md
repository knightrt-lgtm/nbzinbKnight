# nbzinbKnight: Negative Binomial and Zero-Inflated Negative Binomial Regression

Fits negative binomial (NB2) regression by iteratively reweighted least
squares (IRLS) for the coefficients, alternated with Newton-Raphson
updates for the dispersion parameter theta. Fits zero-inflated negative
binomial regression by an EM algorithm whose M-steps are IRLS fits.
Results are checked against MASS::glm.nb() and pscl::zeroinfl().
Developed as the STA 651 final project.

## Author

**Maintainer**: Robert Knight <knightrt@miamioh.edu>
