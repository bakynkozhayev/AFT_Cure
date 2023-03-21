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

# function to simulate datasets and fit the models =============================
sim_fits <-
  \(
    sce = 1,
    df = 2,
    n_datasets = 1e3,
    n_obs = 1e4,
    cure_frac = 0.1,
    # mixture cure flag
    mixture = TRUE,
    # flag for spline constraint on the right boundary
    asymp = FALSE,
    beta = 1,
    theta = NULL,
    cf_dep = FALSE
  ) {
    # scenarios 1-4 ====
    dist <-
      switch(
        sce,
        
        list(
          mixing_par = 0.8,
          shape = c(3, 0.7),
          scale = scale_tr(c(0.1, 0.1), c(3, 1.6))
        ),
        list(
          mixing_par = 0.5,
          shape = c(3, 0.7),
          scale = scale_tr(c(1, 1), c(1.5, 0.5))
        ),
        list(
          mixing_par = 0.26,
          shape = c(3, 0.7),
          scale = scale_tr(c(0.02, 0.5), c(3, 0.7))
        ),
        list(
          mixing_par = 0.5,
          shape = c(1.2, 1.2),
          scale = scale_tr(c(0.1, 0.1), c(1.2, 1.2))
        )
      )
    # ====
    datasets <- list()
    for (dataset in 1:n_datasets) {
      dat <- sim_cure_Weibull(
        n = n_obs,
        shape = dist$shape,
        scale = dist$scale,
        mixing_par = dist$mixing_par,
        cure_frac = cure_frac,
        beta = beta,
        theta = theta,
        cf_dep = cf_dep
      )
      datasets[dataset] <- dat
    }
    for (dataset in 1:n_datasets) {
      fit <- aft_mixture(
        Surv(observed_time, delta) ~ X,
        data = datasets[[dataset]],
        df = df,
        mixture = mixture,
        cure = asymp,
        use.gr = TRUE
      )
    }
  }
test <-
  t(as.data.frame(
    sim_fits(1, 2, 10, 1e3, 0.1, TRUE, FALSE),
    row.names = c(
      "sim_number",
      "bias_beta",
      "bias_relative beta",
      "coverage_beta",
      "bias_cure frac",
      "bias_relative_cure_frac",
      "coverage_cure_frac",
      "AIC",
      "BIC"
    ),
    fix.empty.names = FALSE
  ))
