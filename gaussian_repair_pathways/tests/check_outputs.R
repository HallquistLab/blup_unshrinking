#!/usr/bin/env Rscript
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
folder <- dirname(dirname(normalizePath(script)))
rd <- file.path(folder,"results")
read <- function(n) read.csv(file.path(rd,paste0(n,".csv")))
checks <- read("equivalence_checks")
stopifnot(nrow(checks)==4L,max(as.matrix(checks[-1]))<1e-9)
res <- read("structural_results")
stopifnot(all(res$status_code==0),all(is.finite(res$estimate)),all(res$se>0))
ex <- read("mplus_export_check")
stopifnot(ex$max_likelihood_score_difference<.02,ex$max_Fuller_coefficient_difference<.005)
mc <- read("monte_carlo_summary")
stopifnot(nrow(mc)==10L,all(mc$attempts==mc$successes),all(mc$coverage>=0 & mc$coverage<=1))
# Conditional statistical check with a broad, prespecified MC allowance.
repair <- mc[mc$method=="REPAIR / stabilized",]
stopifnot(all(abs(repair$bias)<.04+4*repair$bias_MCSE))
zero <- read("zero_ME_check"); stopifnot(max(abs(zero$OLS-zero$Fuller))<1e-10)
source(file.path(folder,"R/helpers.R"))
bad <- try(gaussian_person(c(1,2),matrix(1,2,2),diag(2),c(0,0),diag(2)),silent=TRUE)
stopifnot(inherits(bad,"try-error"))
for (software in c("brms","rstanarm")) {
  d <- read(paste0(software,"_diagnostics"))
  s <- read(paste0(software,"_summary"))
  stopifnot(sum(d$num_divergent)==0,sum(d$num_max_treedepth)==0,max(s$rhat,na.rm=TRUE)<1.01)
}
for (m in c("stress","factors")) {
  x <- readLines(file.path(folder,"mplus",paste0(m,".out")),warn=FALSE)
  stopifnot(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY",x,fixed=TRUE)))
}
cat("PASS: four pathway identities, exported Mplus scores, Fuller limit, GLS/MC outputs, identification guard, and software diagnostics.\n")
