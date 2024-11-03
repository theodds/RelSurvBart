## Compile the STAN ----

coxpe_stan <- BART4RS::stanmodels$PEHaz


data_generation_fun <- function(true_model = true_model){
  
  
  if (true_model == "coxpe"){
    
    ## Preprocess data ----
    
    #coxpe_formula <- Surv(event_time, status) ~ ns(age, 5) + ns(wbc,5) + ns(tpi,5) + sex - 1
    coxpe_formula <- Surv(event_time, status) ~ age + wbc + tpi + sex - 1
    train_frame   <- model.frame(coxpe_formula, data = leuk_data)
    X_train       <- model.matrix(coxpe_formula, data = leuk_data)
    s_train       <- model.response(train_frame)
    Y_train       <- s_train[,1]
    delta_train   <- s_train[,2]
    grid          <- c(0, quantile(Y_train, c(.1, .2, .3, .4, .5, .6, .7, .8, .9)))
    stan_data     <- list(X = X_train,
                          y = Y_train,
                          delta = delta_train,
                          N = nrow(X_train),
                          P = ncol(X_train),
                          grid_size = length(grid),
                          grid = grid,
                          pop_haz = rep(0, nrow(X_train)))
    ## Fit ----
    
    fitted_pe     <- rstan::sampling(coxpe_stan,
                                     data = stan_data,
                                     chain = 2,
                                     cores = 1,
                                     iter = 1000)
    
    lambda_param     <- summary(fitted_pe, pars="lambda")$summary[,"mean"]
    linear_predictor <- summary(fitted_pe, pars="eta")$summary[,"mean"]
    
    
    X <- model.matrix(~ age + wbc + tpi + sex - 1 ,
                                           data = leuk_data)
                        
    ## We want to simulate from piecewise exponential proportional hazard model.
    
    C         <- 14 * runif(nrow(X))
    
    t_disease <- array(NA , nrow(X))
    for(j in 1:nrow(X)){
      
      t_disease[j] <- PWEXP::rpwexp(1,  rate = (lambda_param*exp(linear_predictor[j])), breakpoint = grid[-1])
    }
    
    t_pop <- rexp(n = nrow(X), rate = leuk_data$haz_rate)
    t_both <- pmin(t_disease, t_pop)
    
    delta <- ifelse(t_both < C, 1, 0)
    Y <- pmin(t_both, C)
    
    return(leuk_data %>% mutate(event_time = Y, status = delta, eta = linear_predictor))
    
  } else if (true_model == "weibcoxph"){
    
    weibull_formula <-Surv(event_time, status) ~ ns(age, 5) +
      ns(log(1 + wbc),5) + ns(tpi,5) + sex
    
    weibull_fit <-
      weibull_coxph(
        weibull_formula,
        data = leuk_data,
        pop_haz = ,
        test_data = leuk_data,
        verbose = FALSE
      )
    
    ## Get Coefficients and stuff --------------------------------------------------
    
    weibull_samples <- as.matrix(weibull_fit)
    shape_parameter <- mean(weibull_samples[,"weibull_shape"]) ## 0.574
    scale_parameter <- mean(weibull_samples[,"weibull_scale"]) ## 4.82
    beta_weibull    <- colMeans(weibull_samples[,grep("beta", colnames(weibull_samples))])
    
    X <- model.matrix(~ ns(age, 5) +
                        ns(log(1 + wbc),5) + ns(tpi,5) + sex - 1,
                      data = leuk_data)
    
    linear_predictor <- X %*% beta_weibull %>% as.numeric()
    
    ## We want to simulate from shape / scale * (t / scale)^(shape - 1) * exp(linear_predictor)
    ## shape * exp(linear_predictor) / (scale^(shape)) * t^(shape - 1)
    ## shape * (exp(linear_predictor / shape) / scale)^shape * t^(shape - 1)
    
    C         <- 14 * runif(nrow(X))
    t_disease <- rweibull(n = nrow(X),
                          shape = shape_parameter,
                          scale = scale_parameter /
                            exp(linear_predictor / shape_parameter))
    t_pop     <- rexp(n = nrow(X), rate = leuk_data$haz_rate)
    t_both    <- pmin(t_disease, t_pop)
    
    delta     <- ifelse(t_both < C, 1, 0)
    Y         <- pmin(t_both, C)
    
    return(leuk_data %>% mutate(event_time = Y, status = delta, eta = linear_predictor))
    
  } else if (true_model == "weibcoxph_linear"){
    weibull_formula <-Surv(event_time, status) ~ age + wbc + tpi + sex
    
    weibull_fit <-
      weibull_coxph(
        weibull_formula,
        data = leuk_data,
        pop_haz = ,
        test_data = leuk_data,
        verbose = FALSE
      )
    
    ## Get Coefficients and stuff --------------------------------------------------
    
    weibull_samples <- as.matrix(weibull_fit)
    shape_parameter <- mean(weibull_samples[,"weibull_shape"]) ## 0.574
    scale_parameter <- mean(weibull_samples[,"weibull_scale"]) ## 4.82
    beta_weibull    <- colMeans(weibull_samples[,grep("beta", colnames(weibull_samples))])
    
    X <- model.matrix(~ age + wbc + tpi + sex - 1,
                      data = leuk_data)
    
    linear_predictor <- X %*% beta_weibull %>% as.numeric()
    
    ## We want to simulate from shape / scale * (t / scale)^(shape - 1) * exp(linear_predictor)
    ## shape * exp(linear_predictor) / (scale^(shape)) * t^(shape - 1)
    ## shape * (exp(linear_predictor / shape) / scale)^shape * t^(shape - 1)
    
    C         <- 14 * runif(nrow(X))
    t_disease <- rweibull(n = nrow(X),
                          shape = shape_parameter,
                          scale = scale_parameter /
                            exp(linear_predictor / shape_parameter))
    t_pop     <- rexp(n = nrow(X), rate = leuk_data$haz_rate)
    t_both    <- pmin(t_disease, t_pop)
    
    delta     <- ifelse(t_both < C, 1, 0)
    Y         <- pmin(t_both, C)
    
    return(leuk_data %>% mutate(event_time = Y, status = delta, eta = linear_predictor))
  } else if (true_model == "coxpe_bart"){
    
    coxpe_BART_formula <-Surv(event_time, status) ~ age + sex + wbc + tpi 
    
    fit          <- BART4RS::coxpe_bart(formula = coxpe_BART_formula,
                                        data = leuk_data,
                                        pop_haz = rep(0, nrow(leuk_data)),
                                        num_save = 20,
                                        num_thin = 10, num_burn = 100)
    eta          <- colMeans(fit$lambda_train)
    lambda_param <- colMeans(fit$hazard)
    
    
    X            <- model.matrix(~ age + sex + wbc + tpi - 1,
                                 data = leuk_data)
    
    C         <- 14 * runif(nrow(X))
    t_disease <- array(NA , nrow(X))
    for(j in 1:nrow(X)){
      
      t_disease[j] <- PWEXP::rpwexp(1,  rate = (lambda_param*exp(eta[j])), breakpoint = fit$time_grid[-c(1, length(fit$time_grid))])
    }
    t_pop     <- rexp(n = nrow(X), rate = leuk_data$haz_rate)
    t_both    <- pmin(t_disease, t_pop)
    delta     <- ifelse(t_both < C, 1, 0)
    Y         <- pmin(t_both, C)
    
    return(leuk_data %>% mutate(event_time = Y, status = delta, eta = eta))
    
  }
  
}
