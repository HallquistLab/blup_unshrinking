#!/usr/bin/env Rscript
# End-to-end timings, excluding compilation. Run from the repository root.
source('R/fuller_mlm.R')
source('R/stage2_estimators.R')
Rcpp::sourceCpp('src/fuller_mlm_kernels.cpp', cacheDir = file.path(tempdir(), 'fuller_mlm_rcpp_cache'))
args <- commandArgs(trailingOnly = TRUE)
with_fsr <- '--with-fsr' %in% args
args <- args[args != '--with-fsr']
if (with_fsr) {
  source('vh_fsr/R/vh_fuller_q.R')
  Rcpp::sourceCpp('vh_fsr/src/vh_fuller_kernels.cpp',
                  cacheDir = file.path(tempdir(), 'vh_fsr_rcpp_cache'))
}
n <- if (length(args)) as.integer(args[1L]) else 1600L
reps <- if (length(args) > 1L) as.integer(args[2L]) else 10L
stopifnot(is.finite(n), n >= 8L, is.finite(reps), reps >= 1L)
set.seed(241002)
x <- matrix(rnorm(n * 2L), n)
h <- runif(n, 0.5, 1.5)
omega_x <- array(rep(c(0.12, 0.02, 0.02, 0.1), n) * rep(h, each = 4L), c(2L, 2L, n))
y <- 0.5 * x[, 1L] - 0.25 * x[, 2L] + rnorm(n)
dat <- data.frame(y = y, x0 = x[, 1L], x1 = x[, 2L],
                  s11 = omega_x[1L, 1L, ], s12 = omega_x[1L, 2L, ],
                  s22 = omega_x[2L, 2L, ], syy = 0.05 * h)
settings <- list(c('modified', 'modified'), c('fuller', 'modified'), c('fuller', 'fuller'))
runners <- list(
  existing_dual_variants = function() fit_fuller_dual_variants(
    dat, 'y', 'x0', 'x1', 's11', 's12', 's22', 'syy'),
  mlm_r_variants = function() fit_fuller_mlm_variants(y, x, omega_x, dat$syy, backend = 'r'),
  mlm_rcpp_variants = function() fit_fuller_mlm_variants(y, x, omega_x, dat$syy, backend = 'rcpp')
)
# Optional comparison with the separate, untracked factor-score prototype.
if (with_fsr) {
  runners$fsr_rcpp_variants <- function() lapply(settings, function(s) fit_vh_fuller_q(
    y, x, omega_x, dat$syy, preliminary_moment = s[1L], variance_bread = s[2L],
    backend = 'rcpp'))
}
results <- lapply(runners, function(run) {
  run() # warm-up
  median(replicate(reps, system.time(run())[['elapsed']]))
})
seconds <- unlist(results)
print(data.frame(method = names(seconds), n = n, repetitions = reps,
                 median_seconds = unname(seconds),
                 speedup_vs_existing = unname(seconds[1L] / seconds)), row.names = FALSE)
