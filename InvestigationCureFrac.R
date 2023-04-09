## Cure fraction investigation============================================
beta = 1
df = 5
cure_frac = 0.95
a <-
  sim_cure_Weibull(
    cure_frac = cure_frac,
    beta = 1,
    mixing_par = 0.8,
    shape = c(3, 0.7),
    scale = scale_tr(c(0.1, 0.1), c(3, 1.6))
  )
sum(a$delta) / 1e4
dtasets <- rbind(dtasets, a)
# events
plot(density(
  subset(a, subset = delta == 1)$observed_time,
  from = 0,
  to = 10
), main = "pdf event")
# censor dist
plot(density(
  subset(a, subset = delta == 0)$observed_time,
  from = 0,
  to = 10
), main = "pdf censor")



fit <- aft_mixture(
  Surv(observed_time, delta) ~ X,
  data = a,
  mixture = TRUE,
  cure.formula = Surv(observed_time, delta) ~ 1,
  df = df
)
summary(fit)
expit(coef(fit)[2])


bias_beta <- as.numeric(beta - coef(fit)[1])
bias_relative_beta <-
  as.numeric((coef(fit)[1] - beta) / beta * 100)

bias_cure_frac <-
  as.numeric(cure_frac - expit(coef(fit)[2]))
bias_relative_cure_frac <-
  as.numeric((expit(coef(fit)[2]) - cure_frac) / cure_frac * 100)

coverage_beta <- beta  <=
  as.numeric(coef(fit)[1]) + qnorm(0.975) * sqrt(fit@vcov[1, 1]) &
  beta >= as.numeric(coef(fit)[1]) + qnorm(0.025) * sqrt(fit@vcov[1, 1])

coverage_cure_frac <- cure_frac <=
  expit(as.numeric(coef(fit)[2]) + qnorm(0.975) * sqrt(fit@vcov[2, 2])) &
  cure_frac >=
  expit(as.numeric(coef(fit)[2]) + qnorm(0.025) * sqrt(fit@vcov[2, 2]))


AIC <- AIC(fit)
BIC <- BIC(fit)

res <- rbind(
  res,
  data.frame(
    cure_frac = cure_frac,
    df = df,
    bias_beta = bias_beta,
    bias_relative_beta = bias_relative_beta,
    coverage_beta = coverage_beta,
    bias_cure_frac = bias_cure_frac,
    bias_relative_cure_frac = bias_relative_cure_frac,
    coverage_cure_frac = coverage_cure_frac,
    AIC = AIC,
    BIC = BIC
  )
)
res



plot(
  x = seq(1e-4, 10, length.out = 1e4),
  predict(
    fit,
    newdata =
      data.frame(
        observed_time = seq(1e-4, 10, length.out = 1e4),
        X = rep(0, 1e4)
      ),
    type = "surv"
  ),
  type = "l",
  col = "blue",
  ylim = c(0, 1),
  main = "survival"
)

par(new = TRUE)

plot(
  x = seq(1e-4, 10, length.out = 1e4),
  predict(
    fit,
    newdata =
      data.frame(
        observed_time = seq(1e-4, 10, length.out = 1e4),
        X = rep(1, 1e4)
      ),
    type = "surv"
  ),
  type = "l",
  col = "red",
  ylim = c(0, 1)
)

par(new = TRUE)

plot(survfit(Surv(observed_time, delta) ~ X, data = a),
     col = c("purple", "pink"),
     ylim = c(0, 1))
abline(h = cure_frac, col = "darkgreen")
abline(h = expit(coef(fit)[2]), col = "black")


##==============================================================================

suspect <- dtasets[(20e4 + 1):21e4,]
suspect_fit <- aft_mixture(
  Surv(observed_time, delta) ~ X,
  data = suspect,
  mixture = TRUE,
  cure.formula = Surv(observed_time, delta) ~ 1,
  df = df
)
summary(suspect_fit)
1 - coef(suspect_fit)[1]

expit(-coef(suspect_fit)[2])

## ======================================================
# Coverage cure fraction, does not count NA rows
sum(subset(results_cf_0.9, subset = (model == 1) & (sce == 1) &
             (df == 4))$coverage_cure_frac,
    na.rm = TRUE) / sum(complete.cases(subset(
      results_cf_0.9, subset = (model == 1) &
        (sce == 1) &
        (df == 4)
    )))


x <- datasets_sce_1_cf_0.9[[13]]

results_cf_0.9[results_cf_0.9$df == 4 &
                 results_cf_0.9$sim_num == 13 &
                 results_cf_0.9$sce == 1 & results_cf_0.9$model == 2, ]
sum(complete.cases(subset(results_cf_0.9, subset = model == 1 & sce == 1 & df == 10)))

newfit <-
  aft_mixture(
    Surv(observed_time, delta) ~ X,
    cure.formula = Surv(observed_time, delta) ~ 1,
    data = x,
    df = 4,
    cure = TRUE,
    mixture = TRUE,
    use.gr = FALSE
  )
coef(newfit)[[1]]-1
summary(newfit)
expit(coef(newfit))
sum(x$delta)

plot(
  x = seq(1e-4, 10, length.out = 1e4),
  predict(
    newfit,
    newdata =
      data.frame(
        observed_time = seq(1e-4, 10, length.out = 1e4),
        X = rep(0, 1e4)
      ),
    type = "surv"
  ),
  type = "l",
  col = "blue",
  ylim = c(0.87, 1),
  main = "survival"
)

par(new = TRUE)

plot(
  x = seq(1e-4, 10, length.out = 1e4),
  predict(
    newfit,
    newdata =
      data.frame(
        observed_time = seq(1e-4, 10, length.out = 1e4),
        X = rep(1, 1e4)
      ),
    type = "surv"
  ),
  type = "l",
  col = "red",
  ylim = c(0.87, 1)
)

par(new = TRUE)

plot(
  survfit(Surv(observed_time, delta) ~ X, data = x),
  col = c("purple", "pink"),
  ylim = c(0.87, 1)
)
abline(h = cure_frac, col = "darkgreen")
abline(h = expit(coef(newfit)[2]), col = "pink")


y <- datasets_sce_1_cf_0.9[[14]]

newfit2 <-
  aft_mixture(
    Surv(observed_time, delta) ~ X,
    cure.formula = Surv(observed_time, delta) ~ 1,
    data = y,
    df = 2,
    cure = TRUE,
    mixture = TRUE
  )
summary(newfit2)
expit(coef(newfit2)[2])
sum(y$delta)




cure_frac = 0.9
plot(
  x = seq(1e-4, 10, length.out = 1e4),
  predict(
    newfit2,
    newdata =
      data.frame(
        observed_time = seq(1e-4, 10, length.out = 1e4),
        X = rep(0, 1e4)
      ),
    type = "surv"
  ),
  type = "l",
  col = "blue",
  ylim = c(0.87, 1),
  main = "survival"
)

par(new = TRUE)

plot(
  x = seq(1e-4, 10, length.out = 1e4),
  predict(
    newfit2,
    newdata =
      data.frame(
        observed_time = seq(1e-4, 10, length.out = 1e4),
        X = rep(1, 1e4)
      ),
    type = "surv"
  ),
  type = "l",
  col = "red",
  ylim = c(0.87, 1)
)

par(new = TRUE)

plot(
  survfit(Surv(observed_time, delta) ~ X, data = y),
  col = c("purple", "pink"),
  ylim = c(0.87, 1)
)
abline(h = cure_frac, col = "darkgreen")
abline(h = expit(coef(newfit2)[2]), col = "black")
