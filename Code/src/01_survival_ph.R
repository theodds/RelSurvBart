## Load ------------------------------------------------------------------------

source("load.R")
data("LeukSurv", package = "spBayesSurv")
data("leuk_data", package = "BART4RS")

SEED <- digest::digest2int("01_survival_ph.R")

## Fit the proportional hazards BART model -------------------------------------

filename <- "cache/01_fitted_coxpe.rds"
if(file.exists(filename)) {
  fitted_coxpe <- readRDS(filename)
} else {
  set.seed(SEED + 1)
  fitted_coxpe <- BART4RS::coxpe_bart(
    formula = Surv(event_time, status) ~ age + sex + wbc + tpi,
    data = leuk_data,
    pop_haz = leuk_data$haz_rate,
    num_burn = 2500,
    num_save = 2500,
    num_thin = 10
  )
  saveRDS(fitted_coxpe, filename)
}

## Compute the leave-one-out expected log predictive density (ELPD) ------------

loo_ph <- loo(fitted_coxpe$loglik_obs)

## Compute posterior projections using the possum package ----------------------

r <- fitted_coxpe$r_train
r_hat <- fitted_coxpe$r_train_hat

possum_summary <- possum::additive_summary(
  r_hat ~ s(age) + factor(sex) + s(logwbc) + s(tpi),
  fhatSamples = t(r),
  fhat = r_hat,
  df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
)

possum_summary_nosex <- possum::additive_summary(
  r_hat ~ s(age) + s(wbc) + s(tpi),
  fhatSamples = t(r),
  fhat = r_hat,
  df = LeukSurv
)

possum_summary_noage <- possum::additive_summary(
  r_hat ~ s(wbc) + factor(sex) + s(tpi),
  fhatSamples = t(r),
  fhat = r_hat,
  df = LeukSurv
)

possum_summary_notpi <- possum::additive_summary(
  r_hat ~ s(wbc) + factor(sex) + s(age),
  fhatSamples = t(r),
  fhat = r_hat,
  df = LeukSurv
)

possum_summary_nowbc <- possum::additive_summary(
  r_hat ~ s(tpi) + factor(sex) + s(age),
  fhatSamples = t(r),
  fhat = r_hat,
  df = LeukSurv
)

# Plot the additive summary ----------------------------------------------------

add_sum <- additive_summary_plot_2(possum_summary) + xlab("") + theme_bw()
plot(add_sum)

# Compute and plot the summary R-squareds --------------------------------------

rsq_df <- data.frame(Rsq = c(
  possum_summary$summaryRsq,
  possum_summary_nosex$summaryRsq,
  possum_summary_noage$summaryRsq,
  possum_summary_notpi$summaryRsq,
  possum_summary_nowbc$summaryRsq
), Model = rep(c(
  "All",
  "Without Sex",
  "Without Age",
  "Without TPI",
  "Without WBC"
), each = 2500))

p_rsq <- ggplot(rsq_df) +
  geom_density(aes(x = Rsq, y = after_stat(density), fill = Model), color = 'white', alpha = 0.3) +
  xlab("Summary $R^2$") +
  ylab("Density") +
  theme_bw()

plot(p_rsq)


## Plot the baseline hazard function over time ---------------------------------

g <- c(lag(fitted_coxpe$time_grid)[-c(1:2)], 20)
adjusted_hazard <- fitted_coxpe$base_haz

haz_df <- data.frame(
  haz_hat = colMeans(adjusted_hazard),
  steps = fitted_coxpe$time_grid[-1],
  widths = g,
  lower = apply(adjusted_hazard, 2, \(x) quantile(x, 0.025)),
  upper = apply(adjusted_hazard, 2, \(x) quantile(x, 00.975))
)

p_1 <- ggplot(haz_df, aes(x = steps, y = haz_hat)) +
  geom_step() +
  geom_step(aes(y = lower), lty = 2) +
  geom_step(aes(y = upper), lty = 2) +
  xlab("Time") +
  ylab("") +
  xlim(1, 14) +
  ylim(0, 1.15) +
  theme_bw()

p_2 <- ggplot(haz_df, aes(x = steps, y = haz_hat)) +
  geom_step() +
  geom_step(aes(y = lower), lty = 2) +
  geom_step(aes(y = upper), lty = 2) +
  xlab("Time") +
  ylab("") +
  xlim(0, 1) +
  theme_bw()

p_3 <- ggplot(haz_df, aes(x = steps, y = haz_hat)) +
  geom_step() +
  geom_step(aes(y = lower), lty = 2) +
  geom_step(aes(y = upper), lty = 2) +
  xlab("Time") +
  ylab("Hazard") +
  theme_bw()

hazplot <- gridExtra::grid.arrange(p_3, p_2, p_1, nrow = 1)

## Make Pretty Plots -----------------------------------------------------------

tikzprint(add_sum, "01_add_sum", height = 4)
tikzprint(p_rsq, "01_p_rsq", height = 4)
tikzprint(hazplot, "01_hazplot", height = 3.5)
