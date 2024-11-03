nph_data_gen <- function(){
  ## Fit the model ---------------------------------------------------------------
  
  coxnph_BART_formula <-Surv(event_time, status) ~ age + sex + wbc + tpi 
  
  #set.seed(123098)
  fit <- BART4RS::coxnph_bart(formula = coxnph_BART_formula,
                              data = leuk_data,
                              pop_haz = rep(0, nrow(leuk_data)),
                              num_save = 25,
                              num_thin = 10, num_burn = 10)
  
  ## Extract a ground truth and design matrix ------------------------------------
  
  X            <- model.matrix(~ age + sex + wbc + tpi - 1, data = leuk_data)
  grid_length  <- length(fit$time_grid)
  
  h_t <- function(ind,timepoint){
    mean(fit$hazard[,timepoint] * exp(fit$lambda_train[ ind , timepoint, ]))
  }
  
  tmpf <- function() {
    
    rate_matrix <- matrix(NA, nrow = nrow(X), ncol = (grid_length - 1))
    R_matrix    <- matrix(NA, nrow = nrow(X), ncol = (grid_length - 1))
    S_matrix    <- matrix(NA, nrow = nrow(X), ncol = (grid_length - 1))
    
    for(i in 1:nrow(rate_matrix)) {
      for(t in 1:ncol(rate_matrix)) {
        rate_matrix[i,t] <- h_t(i,t)
      }
      R_matrix[i, ] <-
        exp(-rate_matrix[i, ] * (fit$time_grid[-1] - fit$time_grid[-grid_length]))
    }
    S_matrix[,1] <- R_matrix[,1]
    for(t in 2:ncol(rate_matrix)) {
      S_matrix[,t] <- S_matrix[,t-1] * R_matrix[,t]
    }
    return(list(
      rate_matrix = rate_matrix,
      R_matrix = R_matrix,
      S_matrix = cbind(1, S_matrix)
    ))
  }
  
  c(rate_matrix, R_matrix, S_matrix) %<-% tmpf()
  true_s   <- S_matrix
  
  ## A Function for simulating from the ground truth -----------------------------
  
  gen_data <- function() {
    C         <- 14 * runif(nrow(X))
    t_disease <- array(NA , nrow(X))
    for(i in 1:nrow(X)) {
      t_disease[i] <- PWEXP::rpwexp(
        n = 1,
        rate = rate_matrix[i,],
        breakpoint = fit$time_grid[-c(1, length(fit$time_grid))]
      )
    }
    t_pop     <- rexp(n = nrow(X), rate = leuk_data$haz_rate)
    t_both    <- pmin(t_disease, t_pop)
    delta     <- ifelse(t_both < C, 1, 0)
    Y         <- pmin(t_both, C)
    data      <- leuk_data %>% mutate(event_time = Y, status = delta)
    return(list(data = data, t_pop = t_pop, t_disease = t_disease, C = C))
  }
  
  genned_data <- gen_data()
  sim_data <- genned_data$data
  
  return(data_gen = list(sim_data = sim_data, true_s = true_s, truefit = fit))
}

nph_pe_fits  <- function(sim_data, truefit = truefit){
  
  ## Fit COXPE-BART --------------------------------------------------------------
  
  fit_coxpe_bart <- BART4RS::coxpe_bart(
    formula   = Surv(event_time, status) ~ age + sex + wbc + tpi,
    data      = sim_data,
    pop_haz   = sim_data$haz_rate,
    num_save  = 250,
    num_thin  = 10,
    num_burn  = 1000,
    time_grid = truefit$time_grid
  )
  
  ## Fit COXNPH-BART -------------------------------------------------------------
  
  fit_coxnph <- BART4RS::coxnph_bart(
    formula   = Surv(event_time, status) ~ age + sex + wbc + tpi,
    data      = sim_data,
    pop_haz   = sim_data$haz_rate,
    num_save  = 250,
    num_thin  = 10, num_burn = 1000,
    time_grid = truefit$time_grid
  )
  return(fits = list(fit_coxnph = fit_coxnph, 
                     fit_coxpe_bart = fit_coxpe_bart))
}

get_s_hats   <- function(fits = fits){
  fit_coxnph              <- fits[["fit_coxnph"]]
  fit_coxpe_bart          <- fits[["fit_coxpe_bart"]]
  
  s_hat_coxnph_bart       <- extract_S(fit = fit_coxnph)
  s_hat_coxpe_bart        <- extract_S_PH(fit = fit_coxpe_bart)
                              
  mean_s_hat_coxnph_bart  <- apply(extract_S(fit = fit_coxnph), c(2,3), mean)
  mean_s_hat_coxpe_bart   <- apply(extract_S_PH(fit = fit_coxpe_bart), c(1,3), mean)
  
  return(s_hats = list(s_hat_coxnph_bart = s_hat_coxnph_bart,
                       s_hat_coxpe_bart  = s_hat_coxpe_bart,
                  mean_s_hat_coxnph_bart = mean_s_hat_coxnph_bart,
                  mean_s_hat_coxpe_bart  = mean_s_hat_coxpe_bart))
}

get_credintervals <- function(true_s = true_s, s_hats = s_hats){
  
  s_hat_coxnph_bart   <- s_hats$s_hat_coxnph_bart
  s_hat_coxpe_bart    <- s_hats$s_hat_coxpe_bart
  
  ci_s_hat_coxpe_bart <- array(NA, dim = c(1043, 3, ncol(true_s)))
  cp_indic_coxpe_bart <- matrix(NA, nrow = 1043, ncol = ncol(true_s))
  cp_coxpe_bart       <- array(NA, ncol(true_s))
    for(k in 1:ncol(true_s)){
    ci_s_hat_coxpe_bart[, , k] <- matrix(unlist(apply(s_hat_coxpe_bart[, , k], 1, 
                                                    function(x) bayestestR::ci(x, ci=0.9, method="HDI"))), 
                                       nrow=1043, byrow = T)
    cp_indic_coxpe_bart[,k]    <- ifelse(true_s[,k] > ci_s_hat_coxpe_bart[,2, k] & 
                                           true_s[,k] < ci_s_hat_coxpe_bart[,3, k], 1, 0)
    cp_coxpe_bart[k]           <- sum(cp_indic_coxpe_bart[,k])/1043
    }  
  
  
  
  ci_s_hat_coxnph_bart        <- array(NA, dim = c(1043, 3, ncol(true_s)))
  cp_indic_coxnph_bart        <- matrix(NA, nrow = 1043, ncol = ncol(true_s))
  cp_coxnph_bart              <- array(NA, ncol(true_s))
  for(k in 1:ncol(true_s)){
    ci_s_hat_coxnph_bart[, , k] <- matrix(unlist(apply(s_hat_coxnph_bart[, , k], 2, 
                                                      function(x) bayestestR::ci(x, ci=0.9, method="HDI"))), 
                                         nrow=1043, byrow = T)
    cp_indic_coxnph_bart[,k]    <- ifelse(true_s[,k] > ci_s_hat_coxnph_bart[,2, k] & 
                                           true_s[,k] < ci_s_hat_coxnph_bart[,3, k], 1, 0)
    cp_coxnph_bart[k]           <- sum(cp_indic_coxnph_bart[,k])/1043
  } 
  
  return(cp = list(cp_coxpe_bart  = cp_coxpe_bart,
                   cp_coxnph_bart = cp_coxnph_bart,
                   
                   cp_indic_coxpe_bart  = cp_indic_coxpe_bart,
                   cp_indic_coxnph_bart = cp_indic_coxnph_bart,
                   
                   ci_s_hat_coxpe_bart  = ci_s_hat_coxpe_bart[, 2:3, ],
                   ci_s_hat_coxnph_bart = ci_s_hat_coxnph_bart[ ,2:3, ]))
}

output_fn <- function(out_data = out_data){
  n_sim                  <- out_data[[1]]$n_sim
  rmse_coxnphbart        <- matrix(NA, ncol(out_data[[1]]$true_s), n_sim )
  rmse_coxpebart         <- matrix(NA, ncol(out_data[[1]]$true_s),  n_sim)
  
  cp_coxnph              <- matrix(NA, ncol(out_data[[1]]$true_s),  n_sim)
  cp_coxpe               <- matrix(NA, ncol(out_data[[1]]$true_s),  n_sim)
  
  e_coxnph_bart_all      <- array(NA, dim = c(1043, ncol(out_data[[1]]$true_s), n_sim))
  e_coxpe_bart_all       <- array(NA, dim = c(1043, ncol(out_data[[1]]$true_s), n_sim))
  
  cp_indic_coxnph_bart_all      <- array(NA, dim = c(1043, ncol(out_data[[1]]$true_s), n_sim))
  cp_indic_coxpe_bart_all       <- array(NA, dim = c(1043, ncol(out_data[[1]]$true_s), n_sim))
  
  ## Obtain error for each fitted model for each simulated data and adding it to the dataframe ---
  for(m in 1:n_sim){
    e_coxnph_bart_all[, , m]         <- out_data[[m]]$mean_s_hat_coxnph_bart - out_data[[m]]$true_s
    e_coxpe_bart_all[ , , m]         <- out_data[[m]]$mean_s_hat_coxpe_bart  - out_data[[m]]$true_s
    
    cp_indic_coxnph_bart_all[, , m]  <- out_data[[m]]$cp_indic_coxnph_bart
    cp_indic_coxpe_bart_all[ , , m]  <- out_data[[m]]$cp_indic_coxpe_bart
    
    rmse_coxnphbart[, m] <- colMeans((out_data[[m]]$mean_s_hat_coxnph_bart - out_data[[m]]$true_s)^2)
    rmse_coxpebart[, m]  <- colMeans((out_data[[m]]$mean_s_hat_coxpe_bart - out_data[[m]]$true_s)^2)
    
    cp_coxnph[ , m]      <- out_data[[m]]$cp_coxnph_bart
    cp_coxpe[  , m]      <- out_data[[m]]$cp_coxpe_bart
  }
  e_coxnph_bart_matrix <- apply(e_coxnph_bart_all, c(1,2), mean)
  e_coxpe_bart_matrix  <- apply(e_coxpe_bart_all, c(1,2), mean)
  
  cp_coxnph_bart_matrix <- apply(cp_indic_coxnph_bart_all, c(1,2), mean)
  cp_coxpe_bart_matrix  <- apply(cp_indic_coxpe_bart_all,  c(1,2), mean)
  
  rmse_array_nph <- rowMeans(rmse_coxnphbart)
  rmse_array_pe  <- rowMeans(rmse_coxpebart)
  cp_array_nph   <- rowMeans(cp_coxnph)
  cp_array_pe    <- rowMeans(cp_coxpe)
  
  ## Obtain Boxplot for error in s_hat--
  melt_error_nph  <- cbind(melt(e_coxnph_bart_matrix), Model = "nph")
  melt_error_pe   <- cbind(melt(e_coxpe_bart_matrix), Model = "pe")
  melt_error      <- rbind(melt_error_nph, melt_error_pe)
  boxplot_error_s <- ggplot(melt_error, aes(factor(X2), value, fill=Model)) 
  boxplot_error_s <- boxplot_error_s + geom_boxplot()
  
  ## Obtain Boxplot for cp in s_hat--
  melt_cp_nph  <- cbind(melt(cp_coxnph_bart_matrix), Model = "nph")
  melt_cp_pe   <- cbind(melt(cp_coxpe_bart_matrix), Model = "pe")
  melt_cp      <- rbind(melt_cp_nph, melt_cp_pe)
  boxplot_cp_s <- ggplot(melt_cp, aes(factor(X2), value, fill=Model)) 
  boxplot_cp_s <- boxplot_cp_s + geom_boxplot()
  
  ## Obtain plot of RMSE over time in s_hat--
  time             <- out_data[[1]]$truefit$time_grid
  RMSE_dataforplot <- data.frame(cbind(time, rmse_array_nph, rmse_array_pe)[-1, ])
  rmse_melt        <- melt(RMSE_dataforplot, id.vars = "time")
  rmse_plot        <- ggplot(rmse_melt, aes(x=time, y=value, group = variable, colour=variable))+geom_line()
  
  ## Obtain plot of CP over time in s_hat--
  CP_dataforplot   <- data.frame(cbind(time, cp_array_nph, cp_array_pe)[-1, ])
  cp_melt          <- melt(CP_dataforplot, id.vars = "time")
  cp_plot          <- ggplot(cp_melt, aes(x=time, y=value, group = variable, colour=variable))+geom_line()
  
  return(output = list(boxplot_error_s = boxplot_error_s,
                       boxplot_cp_s    = boxplot_cp_s,
                       rmse_plot       = rmse_plot,
                       cp_plot         = cp_plot,
                       rmse_array_nph  = rmse_array_nph,
                       rmse_array_pe   = rmse_array_pe,
                       cp_array_nph    = cp_array_nph,
                       cp_array_pe     = cp_array_pe))
}

