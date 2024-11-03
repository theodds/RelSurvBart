## Load ------------------------------------------------------------------------

library(loo)
library(ggdist)
library(spBayesSurv)
library(mgcv)
library(rpart)
library(rpart.plot)
library(possum)
library(BART4RS)
library(tidybayes)
library(tidyverse)
library(splines)
library(survival)
library(rstan)
library(PWEXP)
library(bayestestR)
data("leuk_data")
data("LeukSurv")

theme_set(theme_bw())
source("lib/DataGenerationFunction.R")
source("lib/FittedModelFunction.R")
source("lib/Get_Eta_Hat_Fn.R")
source("lib/get_cred_inter_fn.R")
source("lib/Get_output_fn.R")

## Globals ---------------------------------------------------------------------

num_sim <- 200
set.seed(20983)

## Run Simulation --------------------------------------------------------------

run_sim <- function(n_sim = 100, true_model = "coxpe") {
  
  out_data                          <- list(NA, n_sim)
  
  for(m in 1:n_sim){
    data_gen             <- data_generation_fun(true_model = true_model)
    # saveRDS(data_gen, file = "data_gen.RDS")
    eta_centered         <- data_gen$eta - mean(data_gen$eta)
    fits                 <- model_fit(data = data_gen)
    # saveRDS(fits, file = "fits.RDS")
    eta_hats             <- get_eta_hats(fits = fits, data = data_gen)
    cis                  <- get_cis(fits = fits, true_eta_centered = eta_centered, 
                                    data = data_gen)
    
    ## 3rd item --------------------------
    out_data[[m]] <- data.frame(
      true_eta_centered   = eta_centered,
      etahat_centered_coxpebart = as.numeric(eta_hats$eta_hat_cenetered_coxpebart),
      etahat_centered_weibcoxph = as.numeric(eta_hats$eta_hat_cenetered_weibcoxph),
      etahat_centered_coxpe     = as.numeric(eta_hats$eta_hat_cenetered_coxpe),
      true_model                = true_model,
      cp_coxpebart              = cis$cov_prob_coxpebart,
      cp_weibcoxph              = cis$cov_prob_weibcoxph,
      cp_coxpe                  = cis$cov_prob_coxpe,
      lci_coxpebart             = cis$ci_etahat_centered_coxpebart[,"lci_coxpebart"],
      uci_coxpebart             = cis$ci_etahat_centered_coxpebart[,"uci_coxpebart"],
      lci_weibcoxph             = cis$ci_etahat_centered_weibcoxph[,"lci_weibcoxph"],
      uci_weibcoxph             = cis$ci_etahat_centered_weibcoxph[,"uci_weibcoxph"],
      lci_coxpe                 = cis$ci_etahat_centered_coxpe[,"lci_coxpe"],
      uci_coxpe                 = cis$ci_etahat_centered_coxpe[,"uci_coxpe"],
      n_sim                     = n_sim)
    
    rm(data_gen)
    print(m)
  }
  return(out_data = out_data)
}

out_data         <- run_sim(n_sim = num_sim, true_model = "weibcoxph_linear")
final_output     <- output_fn(out_data = out_data)
saveRDS(out_data,     "cache/out_data_weibcoxph_linear.RDS")
saveRDS(final_output, "cache/final_output_weibcoxph_linear.RDS")
