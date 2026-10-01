# Gaussian measurement-model adapters. Fuller itself is sourced unchanged.
load_repair <- function(root) {
  source(file.path(root, "R/core_utils.R"), local = .GlobalEnv)
  source(file.path(root, "R/stage2_estimators.R"), local = .GlobalEnv)
}

simulate_stress <- function(n = 300L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- matrix(c(.64, .168, .168, .36), 2)
  theta <- sweep(matrix(rnorm(2*n), n, 2) %*% chol(G), 2, c(1, .5), "+")
  bw <- rnorm(n)
  person <- data.frame(id = seq_len(n), a0 = theta[,1], slope = theta[,2], stress_b = bw)
  person$depression <- 1 + .25*(person$a0-1) + .6*(person$slope-.5) + rnorm(n, 0, .8)
  long <- do.call(rbind, lapply(seq_len(n), function(j) {
    nj <- c(8L, 24L)[1L + (j %% 2L)]
    x <- as.numeric(scale(rnorm(nj)))
    data.frame(id = j, na = theta[j,1] + .4*bw[j] + theta[j,2]*x + rnorm(nj,0,1.5),
               stress_w = x, stress_b = bw[j])
  }))
  list(long = long, person = person, truth = list(mu = c(1,.5), G = G, sigma2 = 2.25,
       beta_bw = .4, gamma = c(.25,.6)))
}

simulate_factors <- function(n = 500L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  L <- matrix(c(1,.8,1.2,0,0,0, 0,0,0,1,.9,1.1), 6,2)
  G <- matrix(c(.64,.28,.28,1),2)
  R <- diag(c(.8,1,.9,.8,1.1,.9)); R[3,4] <- R[4,3] <- .25
  nu <- c(1,1.5,.5,2,1,1.8)
  eta <- matrix(rnorm(n*2),n,2) %*% chol(G)
  Y <- sweep(eta %*% t(L) + matrix(rnorm(n*6),n,6) %*% chol(R),2,nu,"+")
  # Planned MCAR: every fourth participant lacks indicator 2; all remain identified.
  Y[seq(4L,n,4L),2] <- NA_real_
  dat <- data.frame(id = seq_len(n), Y)
  names(dat)[-1] <- paste0("y",1:6)
  person <- data.frame(id = seq_len(n), a0 = eta[,1], slope = eta[,2],
                      depression = 1 + .25*eta[,1] + .6*eta[,2] + rnorm(n,0,.8))
  list(data = dat, person = person, truth = list(L=L,G=G,R=R,nu=nu,mu=c(0,0),gamma=c(.25,.6)))
}

spd_inverse <- function(x) chol2inv(chol((x+t(x))/2))

gaussian_person <- function(y,H,R,mu,G) {
  stopifnot(length(y)==nrow(H), length(mu)==ncol(H), qr(H)$rank==ncol(H))
  RiH <- solve(R,H)
  J <- crossprod(H,RiH); h <- crossprod(H,solve(R,y))
  Psi <- spd_inverse(J)
  direct <- drop(Psi %*% h)
  P <- spd_inverse(G); C <- spd_inverse(J+P)
  m <- drop(C %*% (h + P %*% mu))
  recovered_Psi <- spd_inverse(solve(C)-P)
  recovered <- drop(recovered_Psi %*% (solve(C,m)-P %*% mu))
  list(direct=direct, Psi=Psi, m=m, C=C, recovered=recovered,
       recovered_Psi=recovered_Psi, J=J, h=h)
}

score_row <- function(id,x) {
  data.frame(id=id,x0=x$direct[1],x1=x$direct[2],
    v00=x$Psi[1,1],v01=x$Psi[1,2],v11=x$Psi[2,2],
    m0=x$m[1],m1=x$m[2],c00=x$C[1,1],c01=x$C[1,2],c11=x$C[2,2],
    recovered0=x$recovered[1],recovered1=x$recovered[2],
    rv00=x$recovered_Psi[1,1],rv01=x$recovered_Psi[1,2],rv11=x$recovered_Psi[2,2])
}

stress_scores <- function(cohort, plugin) {
  d <- cohort$long
  rows <- lapply(split(d,d$id),function(z) {
    H <- cbind(1,z$stress_w)
    y <- z$na-plugin$beta_bw*z$stress_b
    score_row(z$id[1],gaussian_person(y,H,diag(plugin$sigma2,nrow(z)),plugin$mu,plugin$G))
  })
  s <- do.call(rbind,rows); s <- s[order(s$id),]
  merge(s,cohort$person[c("id","a0","slope","depression")],by="id",sort=TRUE)
}

factor_scores <- function(cohort,plugin) {
  s <- do.call(rbind,lapply(seq_len(nrow(cohort$data)),function(j) {
    y <- as.numeric(cohort$data[j,paste0("y",1:6)])
    keep <- is.finite(y)
    score_row(cohort$data$id[j],gaussian_person(y[keep]-plugin$nu[keep],
      plugin$L[keep,,drop=FALSE],plugin$R[keep,keep,drop=FALSE],plugin$mu,plugin$G))
  }))
  merge(s,cohort$person,by="id",sort=TRUE)
}

ols_result <- function(d, x0="x0", x1="x1") {
  f <- lm(reformulate(c(x0,x1),response="depression"),data=d)
  est <- unname(coef(f)[x1]); se <- sqrt(vcov(f)[x1,x1])
  data.frame(estimate=est,se=se,ci_low=est-qnorm(.975)*se,ci_high=est+qnorm(.975)*se,status_code=0L)
}

fuller_result <- function(d,variants=c("stabilized","fuller_equations")) {
  d$zero <- 0
  fit_fuller_dual_variants(d,outcome="depression",predictor_u0="x0",predictor_u1="x1",
    meas11="v00",meas12="v01",meas22="v11",outcome_meas_var="zero",
    predictor_outcome_meas_cov_u0="zero",predictor_outcome_meas_cov_u1="zero",
    variants=variants,skip_internal_scaling=TRUE)
}

analyse_scores <- function(d) {
  ff <- fuller_result(d); ff$method <- paste0("REPAIR / ",ff$fuller_variant)
  out <- dplyr::bind_rows(
    dplyr::mutate(ols_result(d,"a0","slope"),method="Oracle: true predictors"),
    dplyr::mutate(ols_result(d),method="Likelihood-only / naive OLS"),
    dplyr::mutate(ols_result(d,"m0","m1"),method="Conditional posterior mean / naive OLS"),ff)
  out$covered <- out$ci_low <= .6 & out$ci_high >= .6
  out
}

equivalence_check <- function(d,label) {
  r <- d; r$x0 <- d$recovered0; r$x1 <- d$recovered1
  r$v00 <- d$rv00; r$v01 <- d$rv01; r$v11 <- d$rv11
  f1 <- fuller_result(d); f2 <- fuller_result(r)
  data.frame(model=label,
    max_score_difference=max(abs(as.matrix(d[c("x0","x1")])-as.matrix(r[c("x0","x1")]))),
    max_ME_difference=max(abs(as.matrix(d[c("v00","v01","v11")])-as.matrix(r[c("v00","v01","v11")]))),
    max_Fuller_coefficient_difference=max(abs(f1$estimate-f2$estimate)),
    max_Fuller_SE_difference=max(abs(f1$se-f2$se)))
}

run_mplus <- function(path) {
  MplusAutomation::runModels(path,showOutput=FALSE,replaceOutfile="always",local_tmpdir=TRUE,
                            logFile=paste0(path,".run.log"))
  output <- sub("[.]inp$",".out",path)
  lines <- readLines(output,warn=FALSE)
  if (!any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY",lines,fixed=TRUE))) {
    stop("Mplus did not terminate normally: ",output)
  }
  MplusAutomation::readModels(output,quiet=TRUE)
}

mplus_value <- function(model,header,param,level=NULL) {
  p <- model$parameters$unstandardized
  take <- p$paramHeader==header & p$param==param
  if (!is.null(level)) take <- take & p$BetweenWithin==level
  ans <- p$est[take]
  if (length(ans)!=1L || !is.finite(ans)) stop("Cannot identify Mplus parameter: ",header," / ",param)
  ans
}
