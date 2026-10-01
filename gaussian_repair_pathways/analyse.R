#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
folder <- dirname(normalizePath(script)); root <- dirname(folder)
source(file.path(folder,"R/helpers.R")); load_repair(root)
out <- file.path(folder,"results")
read <- function(n) readRDS(file.path(out,paste0(n,".rds")))
csv <- function(x,n) write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE,na="")
stress <- read("stress_data"); factors <- read("factor_data")
scores <- list(); results <- list(); checks <- list(); plugins <- list()
for (backend in c("brms","rstanarm","mplus")) {
  plugin <- read(paste0(backend,"_plugin")); plugins[[backend]] <- plugin
  d <- stress_scores(stress,plugin); scores[[backend]] <- d
  z <- analyse_scores(d); z$model <- backend; results[[backend]] <- z
  checks[[backend]] <- equivalence_check(d,backend)
}
plugin <- read("mplus_factor_plugin")
d <- factor_scores(factors,plugin); scores[["Mplus factors"]] <- d
z <- analyse_scores(d); z$model <- "Mplus factors"; results[["Mplus factors"]] <- z
checks[["Mplus factors"]] <- equivalence_check(d,"Mplus factors")
results <- dplyr::bind_rows(results); checks <- dplyr::bind_rows(checks)
stopifnot(all(results$status_code==0),max(as.matrix(checks[-1]))<1e-9)
saveRDS(scores,file.path(out,"scores.rds")); csv(results,"structural_results"); csv(checks,"equivalence_checks")
csv(do.call(rbind,lapply(names(plugins),function(k) {
  p <- plugins[[k]]; data.frame(model=k,mu0=p$mu[1],mu1=p$mu[2],beta_bw=p$beta_bw,
    sigma2=p$sigma2,G00=p$G[1,1],G01=p$G[1,2],G11=p$G[2,2])
})),"plugin_comparison")

# A genuine saved-score pathway: use Mplus's exported ML factor scores.
# Mplus's printed parameter estimates have three decimals. Save that precision
# limitation separately rather than asserting exact agreement with full-precision scores.
mf <- read("mplus_factor_parsed")
saved <- mf$savedata
if (!all(c("ID","F1","F2") %in% toupper(names(saved)))) stop("Cannot locate Mplus saved factor scores")
names(saved) <- toupper(names(saved))
d <- scores[["Mplus factors"]]
saved <- saved[match(d$id,saved$ID),]
stopifnot(identical(as.numeric(saved$ID),as.numeric(d$id)))
recovered <- d
for(j in seq_len(nrow(d))) {
  C <- matrix(c(d$c00[j],d$c01[j],d$c01[j],d$c11[j]),2)
  P <- solve(plugin$G); J <- solve(C)-P
  x <- solve(J,solve(C,as.numeric(saved[j,c("F1","F2")]))-P %*% plugin$mu)
  recovered$x0[j] <- x[1]; recovered$x1[j] <- x[2]
}
fr <- fuller_result(recovered); fr$model <- "Mplus factors: saved scores"; fr$method <- paste0("REPAIR / ",fr$fuller_variant)
csv(fr,"mplus_saved_score_results")
saved_check <- data.frame(
  max_export_vs_reconstructed_posterior_mean=max(abs(as.matrix(saved[c("F1","F2")])-as.matrix(d[c("m0","m1")]))),
  max_likelihood_score_difference=max(abs(as.matrix(recovered[c("x0","x1")])-as.matrix(d[c("x0","x1")]))),
  max_Fuller_coefficient_difference=max(abs(fr$estimate-fuller_result(d)$estimate)),
  max_Fuller_SE_difference=max(abs(fr$se-fuller_result(d)$se)))
csv(saved_check,"mplus_export_check")
stopifnot(saved_check$max_likelihood_score_difference<.02,saved_check$max_Fuller_coefficient_difference<.005)

# Deliberately incorrect use of the marginal MCMC covariance, for illustration.
plugin <- plugins$brms; d <- scores$brms
bad <- do.call(rbind,lapply(plugin$marginal,function(z) {
  J <- solve(z$C)-solve(plugin$G); ev <- eigen(J,symmetric=TRUE)$values
  i <- match(z$id,d$id)
  if(min(ev)<=0) return(data.frame(id=z$id,positive=FALSE,slope_difference=NA_real_,ME_difference=NA_real_))
  V <- solve(J)
  # z$m is a random deviation; add the same population mean before comparison.
  x <- drop(V %*% solve(z$C,z$m)) + plugin$mu
  data.frame(id=z$id,positive=TRUE,slope_difference=x[2]-d$x1[i],ME_difference=V[2,2]-d$v11[i])
}))
csv(bad,"marginal_covariance_misuse")

# Independent numerical validation with a dense, nondiagonal residual covariance.
set.seed(731)
H <- cbind(1,seq(-1,1,length.out=9)); R <- 1.7*.45^abs(outer(1:9,1:9,"-"))
G <- matrix(c(.8,.13,.13,.4),2); mu <- c(.4,-.2); theta <- c(.7,.6)
K <- solve(crossprod(H,solve(R,H)),t(H)%*%solve(R))
x <- gaussian_person(drop(H%*%theta),H,R,mu,G)
stopifnot(max(abs(K%*%H-diag(2)))<1e-12,max(abs(K%*%R%*%t(K)-x$Psi))<1e-12)
set.seed(732); B <- 20000
err <- matrix(rnorm(B*9),B,9)%*%chol(R)
actual <- err %*% t(K)
covcheck <- data.frame(quantity=c("variance 1","covariance 12","variance 2"),
  analytical=c(x$Psi[1,1],x$Psi[1,2],x$Psi[2,2]),
  empirical=c(var(actual[,1]),cov(actual[,1],actual[,2]),var(actual[,2])))
csv(covcheck,"GLS_covariance_check")
stopifnot(max(abs(covcheck$empirical-covcheck$analytical))<.04)

# Boundary case: no predictor measurement error reduces to OLS coefficients.
zero <- scores$brms; zero$v00 <- zero$v01 <- zero$v11 <- 0
zero_test <- data.frame(OLS=ols_result(zero)$estimate,Fuller=fuller_result(zero)$estimate)
csv(zero_test,"zero_ME_check"); stopifnot(max(abs(zero_test$OLS-zero_test$Fuller))<1e-10)

# Repeated cohorts validate the measurement-error/structural pipeline, not
# repeated Bayesian model fitting. Shared nuisance means/loadings are fixed at
# truth; the stress residual variance is re-estimated by pooling OLS residuals.
reps <- if(length(args)) as.integer(args[1]) else 250L
mc <- vector("list",2*reps)
for (b in seq_len(reps)) {
  s <- simulate_stress(300,20262000+b)
  pieces <- split(s$long,s$long$id)
  sse <- sum(vapply(pieces,function(d)sum(lm.fit(cbind(1,d$stress_w),d$na)$residuals^2),numeric(1)))
  p <- s$truth; p$sigma2 <- sse/(nrow(s$long)-2*nrow(s$person))
  ds <- stress_scores(s,p)
  fs <- simulate_factors(300,20265000+b)
  df <- factor_scores(fs,fs$truth)
  for (k in 1:2) {
    z <- analyse_scores(if(k==1) ds else df)
    z$replication <- b; z$design <- c("Stress: pooled residual variance","Factors: fixed measurement parameters")[k]
    z$sigma2 <- if(k==1) p$sigma2 else NA_real_
    mc[[2*(b-1)+k]] <- z
  }
  if(b%%25L==0L) message("Monte Carlo: ",b," / ",reps)
}
mc <- dplyr::bind_rows(mc); csv(mc,"monte_carlo_results")
summary <- mc |>
  dplyr::mutate(success=status_code==0 & is.finite(estimate) & is.finite(se)) |>
  dplyr::group_by(design,method) |>
  dplyr::summarise(attempts=dplyr::n(),successes=sum(success),
    mean_estimate=mean(estimate[success]),bias=mean(estimate[success]-.6),
    bias_MCSE=sd(estimate[success])/sqrt(sum(success)),
    RMSE=sqrt(mean((estimate[success]-.6)^2)),empirical_SD=sd(estimate[success]),mean_SE=mean(se[success]),
    coverage=mean(covered[success]),coverage_MCSE=sqrt(coverage*(1-coverage)/sum(success)),.groups="drop")
csv(summary,"monte_carlo_summary")
stopifnot(all(mc$status_code==0))
capture.output(sessionInfo(),file=file.path(out,"session_analysis.txt"))
paths <- c(file.path(folder,c("R/helpers.R","run_models.R","analyse.R")),
           file.path(root,c("R/core_utils.R","R/stage2_estimators.R")))
csv(data.frame(file=sub(paste0(root,"/"),"",paths,fixed=TRUE),md5=unname(tools::md5sum(paths))),"source_hashes")
saveRDS(list(seed=20260918,reps=reps,stress_n=300,factor_n=500,mc_n=300),file.path(out,"config.rds"))
message("All pathway and covariance checks passed.")
