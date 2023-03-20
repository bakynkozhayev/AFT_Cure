library(rstpm2)
library(pbmcapply)
library(parallel)
# functions ====================================================================
logit <- \(x) ln(p) - ln(1 - p)
expit <- \(x) 1 / (1 + exp(-x))
# transformation of scale from Crowther Lambert to R Weibul
scale_tr = \(lambda, gamma) c(lambda[1] ^ (-1 / gamma[1]), lambda[2] ^ (-1 /
                                                                          gamma[2]))
# data generation ==============================================================
sim_cure_Weibull <- \(
  n = 1e4,
  shape = c(1, 1),
  scale = c(1, 1),
  # proportion drawn from 1st Weibull
  mixing_par = 0.5,
  # proportion of cured
  cure_frac = 0.9,
  # -log acceleration factor
  beta = 1,
  # cure_frac = 1/(1+exp(-theta*X))
  theta = NULL,
  # cure fraction's dependence on covariates
  cf_dep = FALSE
) {
  X <- rbinom(n, 1, 0.5)
  t_uncured <- vuniroot(
    function(t)
      mixing_par * pweibull(t * exp(-beta * X), shape[1], scale[1],
                            lower.tail = FALSE) +
      (1 - mixing_par) *
      pweibull(t * exp(-beta * X), shape[2], scale[2],
               lower.tail = FALSE) - runif(n),
    lower = rep(-1e6, n),
    upper = rep(1e6, n)
  )$root
  t <- ifelse(rbinom(n, 1, cure_frac), Inf, t_uncured)
  
  
  censor_time = pmin(rexp(n, 0.01), 10)
  observed_time = pmin(t, censor_time)
  delta = observed_time < censor_time
  
  
  data.frame(observed_time = observed_time,
             delta = delta,
             X = X)
}

# functions to simulate datasets and fit various models=========================
# aft mixture cure with constrained splines on right boundary
sim_sets <-
  \(
    n_datasets = 1e3,
    n_obs = 1e4,
    shape = c(1, 1),
    scale = scale_tr(c(1, 1), c(1, 1)),
    mixing_par = 0.5,
    cure_frac = 0.1,
    beta = 1,
    theta = NULL,
    cf_dep = FALSE,
    mixture = TRUE, # mixture cure flag
    asympt = FALSE # constrained splines on the right boundary flag
  ) {
    pbmcmapply(FUN)
    dat <- sim_cure_Weibull(n = n_obs,
                            shape = shape,
                            scale = scale,
                            mixing_par = mixing_par,
                            cure_frac = cure_frac,
                            beta = beta,
                            theta = theta,
                            cf_dep = cf_dep
                            )
    fit <- ifelse()
  } 