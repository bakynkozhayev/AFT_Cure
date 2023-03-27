RNGkind("L'Ecuyer-CMRG")
set.seed(1)
# number of observations = 1e4, number of datasets = 100, df = 2:9 
results_cf_0.9 <- sim_fits(cure_frac = 0.9, n_datasets = 100, cores = 30)
results_cf_0.1 <- sim_fits(cure_frac = 0.1, n_datasets = 100, cores = 30)
results_cf_0.5 <- sim_fits(cure_frac = 0.5, n_datasets = 100, cores = 30)
results_cf_0 <- sim_fits(cure_frac = 0, n_datasets = 100, cores = 30)