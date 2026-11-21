# Benchmark: nbzinbKnight vs the reference implementations
# (project description, Section 5 and Table 4).
#
# For each data set (pscl::bioChemists, AER::NMES1988) and each model
# (NB vs MASS::glm.nb, ZINB vs pscl::zeroinfl) this script builds a table of
# our estimates and the reference estimates side by side (coefficients, SEs,
# theta or log(theta), log-likelihood), the maximum absolute difference per
# quantity against the Table 4 tolerance, and computation times.
#
# Run from the package root after installing the package:
#   Rscript inst/benchmark/benchmark.R
# Results are saved to inst/benchmark/benchmark_results.rds for the report.
# The reference functions are used here ONLY for comparison.

library(nbzinbKnight)
set.seed(651)

data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

models <- list(
  bioChemists = list(
    data = bioChemists,
    nb = art ~ fem + mar + kid5 + phd + ment,
    zinb = art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ),
  NMES1988 = list(
    data = NMES1988,
    nb = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance
  )
)

# Table 4 tolerances (maximum absolute difference)
tol <- list(
  NB = c(Coefficients = 1e-5, `Standard errors` = 1e-5,
         Theta = 1e-4, `Log-likelihood` = 1e-4),
  ZINB = c(Coefficients = 1e-3, `Standard errors` = 1e-3,
           `Log(theta)` = 1e-3, `Log-likelihood` = 1e-3)
)

ref_nb_control <- stats::glm.control(epsilon = 1e-10, maxit = 100)
ref_zinb_control <- pscl::zeroinfl.control(reltol = 1e-12)

# Side-by-side table: one row per parameter, columns ours / reference / diff
side_by_side <- function(quantity, term, ours, ref) {
  data.frame(quantity = quantity, term = term, ours = unname(ours),
             reference = unname(ref), abs_diff = abs(unname(ours - ref)),
             row.names = NULL)
}

# Summary: maximum absolute difference per quantity vs its tolerance
max_diffs <- function(tab, tolerance) {
  md <- tapply(tab$abs_diff, tab$quantity, max)[names(tolerance)]
  data.frame(quantity = names(tolerance), max_abs_diff = unname(md),
             tolerance = unname(tolerance), pass = unname(md <= tolerance),
             row.names = NULL)
}

# Median elapsed time (seconds) over `reps` runs
time_it <- function(expr_fun, reps) {
  median(vapply(seq_len(reps),
                function(i) system.time(expr_fun())[["elapsed"]], numeric(1)))
}

results <- list()
for (ds in names(models)) {
  m <- models[[ds]]
  d <- m$data

  # ---- NB ----
  ours <- nbreg(m$nb, data = d)
  ref <- MASS::glm.nb(m$nb, data = d, control = ref_nb_control)
  nb_tab <- rbind(
    side_by_side("Coefficients", names(coef(ours)), coef(ours), coef(ref)),
    side_by_side("Standard errors", names(coef(ours)),
                 sqrt(diag(vcov(ours))), sqrt(diag(vcov(ref)))),
    side_by_side("Theta", "theta", ours$theta, ref$theta),
    side_by_side("Theta SE", "SE(theta)", ours$se_theta, ref$SE.theta),
    side_by_side("Log-likelihood", "logLik", as.numeric(logLik(ours)),
                 as.numeric(logLik(ref)))
  )
  nb_time <- c(
    nbzinbKnight = time_it(function() nbreg(m$nb, data = d), reps = 10),
    reference = time_it(function() MASS::glm.nb(m$nb, data = d,
                                                 control = ref_nb_control),
                        reps = 10))

  # ---- ZINB ----
  ours_z <- zinbreg(m$zinb, data = d)
  ref_z <- pscl::zeroinfl(m$zinb, data = d, dist = "negbin",
                          control = ref_zinb_control)
  zinb_tab <- rbind(
    side_by_side("Coefficients", names(coef(ours_z)), coef(ours_z),
                 coef(ref_z)),
    side_by_side("Standard errors", names(coef(ours_z)),
                 sqrt(diag(vcov(ours_z))), sqrt(diag(vcov(ref_z)))),
    side_by_side("Log(theta)", "Log(theta)", log(ours_z$theta),
                 log(ref_z$theta)),
    side_by_side("Log(theta) SE", "SE(Log(theta))", ours_z$se_logtheta,
                 ref_z$SE.logtheta),
    side_by_side("Log-likelihood", "logLik", as.numeric(logLik(ours_z)),
                 as.numeric(logLik(ref_z)))
  )
  zinb_time <- c(
    nbzinbKnight = time_it(function() zinbreg(m$zinb, data = d), reps = 3),
    reference = time_it(function() pscl::zeroinfl(m$zinb, data = d,
                                                  dist = "negbin",
                                                  control = ref_zinb_control),
                        reps = 3))

  results[[ds]] <- list(
    n = nrow(d),
    nb = list(table = nb_tab, max_diff = max_diffs(nb_tab, tol$NB),
              time = nb_time, iter = ours$iter),
    zinb = list(table = zinb_tab, max_diff = max_diffs(zinb_tab, tol$ZINB),
                time = zinb_time, iter = ours_z$iter,
                monotone = all(diff(ours_z$loglik_trace) >= -1e-10))
  )
}

results$session <- list(date = Sys.Date(), R = R.version.string,
                        nbzinbKnight = as.character(packageVersion("nbzinbKnight")),
                        MASS = as.character(packageVersion("MASS")),
                        pscl = as.character(packageVersion("pscl")))

out <- file.path("inst", "benchmark", "benchmark_results.rds")
saveRDS(results, out)

# Console summary
for (ds in setdiff(names(results), "session")) {
  for (mod in c("nb", "zinb")) {
    r <- results[[ds]][[mod]]
    cat(sprintf("\n== %s, %s (iterations: %d) ==\n", ds, toupper(mod), r$iter))
    print(r$max_diff, digits = 3, row.names = FALSE)
    cat(sprintf("Median time (s): nbzinbKnight %.3f, reference %.3f\n",
                r$time[["nbzinbKnight"]], r$time[["reference"]]))
  }
}
cat("\nSaved", out, "\n")
