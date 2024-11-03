## Load ------------------------------------------------------------------------
library(devtools)

library(loo)
library(ggdist)
library(spBayesSurv)
library(mgcv)
library(rpart)
library(rpart.plot)
library(possum)
library(Batman)
library(BART4RS)
library(tidybayes)
library(splines)
library(survival)
library(rstan)
library(PWEXP)
library(zeallot)
library(bayestestR)
library(reshape)
data("leuk_data")

theme_set(theme_bw())
source("lib/nph_SimFunctions.R")

## Globals ---------------------------------------------------------------------

num_sim <- 200
set.seed(20983)

## Run Simulation --------------------------------------------------------------

run_sim_nph<- function(n_sim = 200) {
  out_data                          <- list(NA, n_sim)
  for(m in 1:n_sim){
    data_gen <- nph_data_gen()
    sim_data <- data_gen$sim_data
    truefit  <- data_gen$truefit
    true_s   <- data_gen$true_s
    fits     <- nph_pe_fits(sim_data = sim_data, truefit = truefit)
    
    s_hats   <- get_s_hats(fits = fits)
    cis      <- get_credintervals(true_s = true_s, s_hats = s_hats)
    
    out_data[[m]] <- list(true_s      = true_s, truefit = truefit,
                          s_hat_coxnph_bart            = s_hats$s_hat_coxnph_bart,
                          s_hat_coxpe_bart             = s_hats$s_hat_coxpe_bart,
                          
                          mean_s_hat_coxnph_bart       = s_hats$mean_s_hat_coxnph_bart,
                          mean_s_hat_coxpe_bart        = s_hats$mean_s_hat_coxpe_bart,
                          
                          cp_coxnph_bart               = cis$cp_coxnph_bart,
                          cp_coxpe_bart                = cis$cp_coxpe_bart,
                          
                          cp_indic_coxnph_bart         = cis$cp_indic_coxnph_bart,
                          cp_indic_coxpe_bart          = cis$cp_indic_coxpe_bart,
                          
                          ci_s_hat_coxnph_bart         = cis$ci_s_hat_coxnph_bart,
                          ci_s_hat_coxpe_bart          = cis$ci_s_hat_coxpe_bart,
                          
                          n_sim                        = n_sim)
    
    print(m)
  }
  return(out_data = out_data)
}

out_data_new <- run_sim_nph(n_sim = num_sim)
saveRDS(out_data, file = "cache/out_data_nph.RDS")

final_output_new <- output_fn(out_data = out_data_new)
saveRDS(final_output_nph, file = "cache/final_output_nph.RDS")
