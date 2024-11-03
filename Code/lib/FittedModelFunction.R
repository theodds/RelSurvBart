model_fit <- function(data = data){
  ## COXPE_BART FIT ------------------------------------------------------------------------
  
  fit_coxpe_bart   <- BART4RS::coxpe_bart(formula = Surv(event_time, status) ~ age + sex + wbc + tpi,
                                          data = data,
                                          pop_haz = data$haz_rate,
                                          num_save = 250,
                                          num_thin = 10,
                                          num_burn = 100)
  
  
  ## Weibull Baseline Cox PH fit -----------------------------------------------------------
  
  fit_weibull_coxph <- BART4RS::weibull_coxph(formula = Surv(event_time, status) ~ age + sex + wbc + tpi,
                                              data = data,
                                              pop_haz = data$haz_rate)
  
  ## COXPE_Linear FIT ------------------------------------
  ## Preprocess data ----
  
  coxpe_formula <- Surv(event_time, status) ~ age + sex + wbc + tpi - 1
  
  
  train_frame <- model.frame(coxpe_formula, data)
  
  X_train <- model.matrix(coxpe_formula, data)
  s_train <- model.response(train_frame)
  Y_train <- s_train[,1]
  delta_train <- s_train[,2]
  
  grid <- c(0, quantile(Y_train, c(.1, .2, .3, .4, .5, .6, .7, .8, .9)))
  
  stan_data <- list(X = X_train,
                    y = Y_train,
                    delta = delta_train,
                    N = nrow(X_train),
                    P = ncol(X_train),
                    grid_size = length(grid),
                    grid = grid,
                    pop_haz = data$haz_rate)
  
  
  ## Fit CoxPE ----
  
  fitted_pe <- rstan::sampling(coxpe_stan,
                               data = stan_data,
                               chain = 2,
                               cores = 1,
                               iter = 1000)
  
  return(fits = list(fit_coxpe_bart = fit_coxpe_bart, 
         fit_weibull_coxph = fit_weibull_coxph,
         fitted_pe = fitted_pe))
}
