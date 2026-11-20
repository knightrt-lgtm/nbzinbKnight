# STA 651 final project: NB and ZINB regression R package (nbzinbKnight)

## Rules (must follow)
- Fitting code in R/ must NOT call optim, optimize/optimise, nlm, nlminb,
  uniroot, glm, glm.fit, lm, lm.fit, lm.wfit, MASS::glm.nb, MASS::theta.ml
  (or any MASS fitting code), pscl::zeroinfl, countreg, glmmTMB, gamlss, VGAM,
  optimx, maxLik, stats4::mle, bbmle, or any other optimizer.
- Allowed in R/: solve, qr, chol, crossprod, %*%; lgamma, digamma, trigamma,
  dnbinom, plogis, qlogis; model.frame, model.matrix, model.response,
  model.offset; the Formula package for two-part formulas;
  numDeriv::hessian for standard errors only.
- NB: IRLS for beta alternated with Newton-Raphson on log(theta)
  (project-docs/STA651_Final_Project.pdf, Section 3.1).
- ZINB: EM algorithm; M-steps are IRLS (weighted logistic for gamma +
  weighted NB with prior weights 1 - tau for beta, theta) (Section 3.2).
  The observed-data log-likelihood must never decrease.
- ZINB standard errors: numerical Hessian via numDeriv::hessian of our own
  log-likelihood in (beta, gamma, log theta). The analytic-Hessian bonus is
  deferred until the base project is finished.
- MASS::glm.nb, pscl::zeroinfl, and glm may be used ONLY in tests/,
  inst/benchmark/, and the report, for comparison.
- Follow the equations in project-docs/STA651_Final_Project.pdf, Section 3.
- Comment estimation code with the equation it implements.
- Explain your reasoning to me; I must be able to explain every line in the
  oral walkthrough.

## How we work
- Work in small steps by the Section 11 schedule checkpoints; check in with me
  at the end of every checkpoint (more often is fine). Do not build large
  chunks in one go.
- Make small git commits with clear messages after each working step
  (at least ~10 over the project).
- Before each commit, search R/ for prohibited functions.

## Workflow commands
- After editing roxygen comments: Rscript -e 'devtools::document()'
- Run tests: Rscript -e 'devtools::test()'
- Full check (must be 0 errors, 0 warnings): Rscript -e 'devtools::check()'
- Report PDF: quarto render vignettes/articles/report.qmd --to pdf
- Website: Rscript -e 'pkgdown::build_site()'
- Close RStudio before updating installed packages (it locks package DLLs).
