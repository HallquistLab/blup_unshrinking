#!/usr/bin/env Rscript
# Run from any directory. Optional argument: data, brms, rstanarm, mplus, or all.
args <- commandArgs(trailingOnly=TRUE)
phase <- if(length(args)) args[1] else "all"
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
folder <- dirname(normalizePath(script)); root <- dirname(folder)
source(file.path(folder,"R/helpers.R"))
out <- file.path(folder,"results"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(folder,"cache"); dir.create(cache,showWarnings=FALSE)
if (phase %in% c("data","all")) {
  saveRDS(simulate_stress(300,20260918),file.path(out,"stress_data.rds"))
  saveRDS(simulate_factors(500,20260919),file.path(out,"factor_data.rds"))
}
cohort <- readRDS(file.path(out,"stress_data.rds")); dat <- cohort$long
dat$id <- factor(dat$id)
write_csv <- function(x,n) write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE)

if (phase %in% c("brms","all")) {
  options(mc.cores=4)
  fit <- brms::brm(na ~ stress_w + stress_b + (1 + stress_w | id),data=dat,
    family=gaussian(),backend="cmdstanr",seed=20260918,chains=4,cores=4,
    iter=2000,warmup=1000,refresh=200,control=list(adapt_delta=.95,max_treedepth=12),
    prior=c(brms::prior(normal(0,2),class="b"),brms::prior(normal(0,3),class="Intercept"),
            brms::prior(exponential(1),class="sd"),brms::prior(exponential(1),class="sigma"),
            brms::prior(lkj(2),class="cor")),
    stan_model_args=list(cpp_options=list(PRECOMPILED_HEADERS="false")),
    file=file.path(cache,"brms_measurement"),file_refit="on_change")
  dr <- posterior::as_draws_df(fit)
  gv <- brms::VarCorr(fit,summary=FALSE)$id$cov
  G <- apply(gv,c(2,3),mean)
  re <- brms::ranef(fit,summary=FALSE)$id
  # Random-effect draws have person names; retain explicit IDs rather than assuming order.
  marginal <- lapply(seq_len(dim(re)[2]),function(j) {
    x <- re[,j,]; list(id=as.integer(dimnames(re)[[2]][j]),m=colMeans(x),C=cov(x))
  })
  plugin <- list(mu=c(mean(dr$b_Intercept),mean(dr$b_stress_w)),beta_bw=mean(dr$b_stress_b),
    sigma2=mean(dr$sigma^2),G=unname(G),marginal=marginal)
  saveRDS(plugin,file.path(out,"brms_plugin.rds"))
  write_csv(as.data.frame(posterior::summarise_draws(posterior::as_draws_array(fit))),"brms_summary")
  sp <- rstan::get_sampler_params(fit$fit,inc_warmup=FALSE)
  diag <- do.call(rbind,lapply(seq_along(sp),function(j) data.frame(chain=j,
    num_divergent=sum(sp[[j]][,"divergent__"]),num_max_treedepth=sum(sp[[j]][,"treedepth__"]>=12))))
  write_csv(diag,"brms_diagnostics")
  writeLines(brms::stancode(fit),file.path(out,"brms_measurement.stan"))
}

if (phase %in% c("rstanarm","all")) {
  options(mc.cores=4)
  fit <- rstanarm::stan_glmer(na ~ stress_w + stress_b + (1 + stress_w | id),data=dat,
    family=gaussian(),seed=20260918,chains=4,cores=4,iter=2000,warmup=1000,refresh=200,
    prior=rstanarm::normal(0,2,autoscale=FALSE),
    prior_intercept=rstanarm::normal(0,3,autoscale=FALSE),
    prior_aux=rstanarm::exponential(1,autoscale=FALSE),
    adapt_delta=.95,control=list(max_treedepth=12))
  saveRDS(fit,file.path(cache,"rstanarm_measurement.rds"))
  dr <- as.matrix(fit); saveRDS(colnames(dr),file.path(out,"rstanarm_parameter_names.rds"))
  cn <- colnames(dr)
  get <- function(nm) { if(!nm %in% cn) stop("Missing parameter ",nm); dr[,nm] }
  G <- matrix(c(mean(get("Sigma[id:(Intercept),(Intercept)]")),
                mean(get("Sigma[id:stress_w,(Intercept)]")),
                mean(get("Sigma[id:stress_w,(Intercept)]")),
                mean(get("Sigma[id:stress_w,stress_w]"))),2)
  plugin <- list(mu=c(mean(get("(Intercept)")),mean(get("stress_w"))),beta_bw=mean(get("stress_b")),
                 sigma2=mean(get("sigma")^2),G=G)
  saveRDS(plugin,file.path(out,"rstanarm_plugin.rds"))
  write_csv(as.data.frame(posterior::summarise_draws(posterior::as_draws_array(fit))),"rstanarm_summary")
  sp <- rstan::get_sampler_params(fit$stanfit,inc_warmup=FALSE)
  diag <- do.call(rbind,lapply(seq_along(sp),function(j) data.frame(chain=j,
    num_divergent=sum(sp[[j]][,"divergent__"]),num_max_treedepth=sum(sp[[j]][,"treedepth__"]>=12))))
  write_csv(diag,"rstanarm_diagnostics")
}

if (phase %in% c("mplus","all")) {
  md <- file.path(folder,"mplus"); dir.create(md,showWarnings=FALSE)
  write.table(cohort$long,file.path(md,"stress.dat"),row.names=FALSE,col.names=FALSE,quote=FALSE,na="-999")
  input <- c("TITLE: REPAIR Gaussian stress measurement model;", "DATA: FILE = stress.dat;",
    "VARIABLE: NAMES = id affect stress_w stress_b;", " USEVARIABLES = affect stress_w stress_b;",
    " CLUSTER = id; WITHIN = stress_w; BETWEEN = stress_b;", "ANALYSIS: TYPE = TWOLEVEL RANDOM;",
    " ESTIMATOR = BAYES; CHAINS = 4; PROCESSORS = 4; BITERATIONS = (10000); BSEED = 20260918;",
    "MODEL:"," %WITHIN%", " s | affect ON stress_w;", " affect (ve);", " %BETWEEN%",
    " affect ON stress_b;", " [affect s];", " affect s;", " affect WITH s;", "OUTPUT: TECH1 TECH8;",
    "SAVEDATA: FILE = stress_scores.dat; SAVE = FSCORES (1000);", " BPARAMETERS = stress_parameters.dat;",
    " FORMAT = F20.10;")
  writeLines(input,file.path(md,"stress.inp"))
  fit <- run_mplus(file.path(md,"stress.inp"))
  saveRDS(fit,file.path(out,"mplus_stress_parsed.rds"))
  write_csv(fit$parameters$unstandardized,"mplus_stress_parameters")
  # Mplus Bayes reports posterior medians. Use one coherent set of plug-in summaries;
  # these are explicitly distinct from the posterior means used for brms/rstanarm.
  p <- function(h,n,l) mplus_value(fit,h,n,l)
  plugin <- list(mu=c(p("Intercepts","AFFECT","Between"),p("Means","S","Between")),
     beta_bw=p("AFFECT.ON","STRESS_B","Between"),sigma2=p("Residual.Variances","AFFECT","Within"),
     G=matrix(c(p("Residual.Variances","AFFECT","Between"),p("AFFECT.WITH","S","Between"),
                p("AFFECT.WITH","S","Between"),p("Variances","S","Between")),2))
  saveRDS(plugin,file.path(out,"mplus_plugin.rds"))

  fc <- readRDS(file.path(out,"factor_data.rds"))
  write.table(fc$data,file.path(md,"factors.dat"),row.names=FALSE,col.names=FALSE,quote=FALSE,na="-999")
  input <- c("TITLE: REPAIR two-factor Gaussian measurement model;", "DATA: FILE = factors.dat;",
    "VARIABLE: NAMES = id y1-y6; USEVARIABLES = y1-y6; IDVARIABLE = id; MISSING = ALL(-999);",
    "ANALYSIS: ESTIMATOR = ML;", "MODEL:", " f1 BY y1@1 y2 y3;", " f2 BY y4@1 y5 y6;",
    " f1 WITH f2;", " y3 WITH y4;", "OUTPUT: TECH1 TECH4;",
    "SAVEDATA: FILE = factor_scores.dat; SAVE = FSCORES; FORMAT = F20.10;")
  writeLines(input,file.path(md,"factors.inp"))
  ff <- run_mplus(file.path(md,"factors.inp"))
  saveRDS(ff,file.path(out,"mplus_factor_parsed.rds"))
  write_csv(ff$parameters$unstandardized,"mplus_factor_parameters")
  p <- function(h,n) mplus_value(ff,h,n)
  L <- matrix(0,6,2); L[1:3,1] <- vapply(paste0("Y",1:3),function(y)p("F1.BY",y),numeric(1))
  L[4:6,2] <- vapply(paste0("Y",4:6),function(y)p("F2.BY",y),numeric(1))
  R <- diag(vapply(paste0("Y",1:6),function(y)p("Residual.Variances",y),numeric(1)))
  R[3,4] <- R[4,3] <- p("Y3.WITH","Y4")
  plugin <- list(L=L,R=R,nu=vapply(paste0("Y",1:6),function(y)p("Intercepts",y),numeric(1)),
    mu=c(0,0),G=matrix(c(p("Variances","F1"),p("F1.WITH","F2"),p("F1.WITH","F2"),p("Variances","F2")),2))
  saveRDS(plugin,file.path(out,"mplus_factor_plugin.rds"))
}
capture.output(sessionInfo(),file=file.path(out,paste0("session_",phase,".txt")))
message("Completed: ",phase)
