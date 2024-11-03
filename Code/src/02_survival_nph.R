## Load ------------------------------------------------------------------------

source("load.R")

data("LeukSurv", package = "spBayesSurv")
data("leuk_data", package = "BART4RS")

SEED <- digest::digest2int("02_survival_nph")

## Fit the non-proportional hazards BART model ---------------------------------

filename <- "cache/02_my_fit.rds"
if (file.exists(filename)) {
  my_fit <- readRDS(filename)
} else {
  set.seed(SEED + 1)
  my_fit <- coxnph_bart(
    Surv(event_time, status) ~ age + sex + wbc + tpi,
    data = leuk_data,
    pop_haz = leuk_data$haz_rate,
    num_burn = 1000,
    num_save = 1000,
    num_thin = 1
  )
  saveRDS(my_fit, filename)
}

## Compute leave-one-out expected log predictive density (ELPD) ----------------

loo_nph <- loo(my_fit$loglik_obs)
print(loo_nph)

## Extract R(b | x) function ---------------------------------------------------

R <- extract_R(my_fit)

## Create visualization function for each bin using possum package -------------

possum_R <- function(j) {
  SE_j <- R[[j]]
  SE_j_hat <- colMeans(SE_j)
  
  possum_summary_SE1 <- possum::additive_summary(
    SE_j_hat ~ s(age) + factor(sex) + s(logwbc) + s(tpi), 
    fhatSamples = t(SE_j), 
    fhat = SE_j_hat, 
    df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
  )
  
  return(possum::additive_summary_plot(possum_summary_SE1) + theme_bw())
}

## Plot the results for bins 1, 10, and 20 -------------------------------------

plot_1 <- possum_R(1)
plot_10 <- possum_R(10)
plot_20 <- possum_R(20)

plot_1 / plot_10 / plot_20


## Create variable importance plot function for a given bin --------------------

possum_varimp <- function(j) {
  SE_j <- R[[j]]
  SE_j_hat <- colMeans(SE_j)
  
  possum_summary_SE1 <- possum::additive_summary(
    SE_j_hat ~ s(age) + factor(sex) + s(logwbc) + s(tpi), 
    fhatSamples = t(SE_j), 
    fhat = SE_j_hat, 
    df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
  )
  possum_summary_noage <- possum::additive_summary(
    SE_j_hat ~ factor(sex) + s(logwbc) + s(tpi), 
    fhatSamples = t(SE_j), 
    fhat = SE_j_hat, 
    df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
  )
  possum_summary_nosex <- possum::additive_summary(
    SE_j_hat ~ s(age) + s(logwbc) + s(tpi), 
    fhatSamples = t(SE_j), 
    fhat = SE_j_hat, 
    df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
  )
  possum_summary_nowbc <- possum::additive_summary(
    SE_j_hat ~ s(age) + factor(sex) + s(tpi), 
    fhatSamples = t(SE_j), 
    fhat = SE_j_hat, 
    df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
  )
  possum_summary_notpi <- possum::additive_summary(
    SE_j_hat ~ s(age) + factor(sex) + s(logwbc), 
    fhatSamples = t(SE_j), 
    fhat = SE_j_hat, 
    df = LeukSurv %>% mutate(logwbc = log(1 + wbc))
  )

  rsq_df <- data.frame(
    Rsq = c(possum_summary_SE1$summaryRsq, possum_summary_nosex$summaryRsq, possum_summary_noage$summaryRsq, possum_summary_notpi$summaryRsq, possum_summary_nowbc$summaryRsq),
    Model = rep(c("All", "Without Sex", "Without Age", "Without TPI", "Without WBC"), each = 1000)
  )

  p_rsq <- ggplot(rsq_df) +
    geom_density(aes(x = Rsq, y = after_stat(density), fill = Model), color = 'white', alpha = 0.3) +
    xlab("Summary $R^2$") +
    ylab("Density") + theme_bw()

  return(p_rsq)
}

## Plot the variable importance for bins 1, 10, and 20 -------------------------

pv1 <- possum_varimp(1)
pv10 <- possum_varimp(10)
pv20 <- possum_varimp(20)

pv1 / pv10 / pv20

## tikzPrint -------------------------------------------------------------------

tikzprint(plot_1 / plot_10 / plot_20, "02_survival_gam")
tikzprint(pv1 / pv10 / pv20, "02_survival_varimp")
