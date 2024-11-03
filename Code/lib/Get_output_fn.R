output_fn <- function(out_data = out_data){
  
  n_sim                  <- out_data[[1]]$n_sim[1]
  rmse_coxpebart         <- array(NA, n_sim)
  rmse_weibcoxph         <- array(NA, n_sim)
  rmse_coxpe             <- array(NA, n_sim)
  
  E_matrix_coxpebart            <- matrix(NA, n_sim, 1043)
  E_matrix_weibcoxph            <- matrix(NA, n_sim, 1043)
  E_matrix_coxpe                <- matrix(NA, n_sim, 1043)
  
  cp_coxpebart           <- array(NA, n_sim)
  cp_weibcoxph           <- array(NA, n_sim)
  cp_coxpe               <- array(NA, n_sim)
  
  for(m in 1: n_sim){
    
    ## Obtain error for each fitted model for each simulated data and adding it to the dataframe -------------
    
    out_data[[m]]$e_coxpebart <- out_data[[m]]$etahat_centered_coxpebart - out_data[[m]]$true_eta_centered
    out_data[[m]]$e_weibcoxph <- out_data[[m]]$etahat_centered_weibcoxph - out_data[[m]]$true_eta_centered
    out_data[[m]]$e_coxpe     <- out_data[[m]]$etahat_centered_coxpe     - out_data[[m]]$true_eta_centered
    
    
    ## Obtain rmse for each fitted model for each simulated data and adding it to the dataframe -------------
    rmse_coxpebart[m]         <- sqrt (mean( (out_data[[m]]$e_coxpebart)^2) )
    rmse_weibcoxph[m]         <- sqrt (mean( (out_data[[m]]$e_weibcoxph)^2) )
    rmse_coxpe[m]             <- sqrt (mean( (out_data[[m]]$e_coxpe)^2) )
    
    cp_coxpebart[m]           <- out_data[[m]]$cp_coxpebart[1]
    cp_weibcoxph[m]           <- out_data[[m]]$cp_weibcoxph[1]
    cp_coxpe[m]               <- out_data[[m]]$cp_coxpe[1]
  }
  
  ## Obtain avg rmse for each fitted model -------------
  avg_rmse_coxpebart          <- mean(rmse_coxpebart)
  avg_rmse_weibcoxph          <- mean(rmse_weibcoxph)
  avg_rmse_coxpe              <- mean(rmse_coxpe)
  
  
  ## Obtain coverage probability for each fitted model -------------
  avg_cp_coxpebart            <- mean(cp_coxpebart)
  avg_cp_weibcoxph            <- mean(cp_weibcoxph)
  avg_cp_coxpe                <- mean(cp_coxpe)
  
  
  ## Prepare data for boxplot and Generate Boxplot of errors of the fitted models -------------
  E <- lapply(out_data, function(x) x[, c("true_model", "e_coxpebart", "e_weibcoxph", "e_coxpe")])
  for(m in 1:n_sim){
    E_matrix_coxpebart[m, ] <- E[[m]]$e_coxpebart
    E_matrix_weibcoxph[m, ] <- E[[m]]$e_weibcoxph
    E_matrix_coxpe[m, ]     <- E[[m]]$e_coxpe
  }
  
  
  boxplot_data      <- data.frame(error_coxpebart = colMeans(E_matrix_coxpebart),
                                  error_weibcox   = colMeans(E_matrix_weibcoxph),
                                  error_coxpe     = colMeans(E_matrix_coxpe))
  
  boxplot_data_long <- reshape(boxplot_data, direction = "long", 
                               varying = names(boxplot_data),
                               v.names = "Value", 
                               timevar = "fitted.model", 
                               times = names(boxplot_data))
  p <- ggplot(data = boxplot_data_long, aes(x=fitted.model, y=Value)) + geom_boxplot()
  
  return(list( p = p, 
               avg_rmse = list(avg_rmse_coxpebart = avg_rmse_coxpebart, 
                               avg_rmse_weibcoxph = avg_rmse_weibcoxph, 
                               avg_rmse_coxpe     = avg_rmse_coxpe),
               avg_cp   = list(avg_cp_coxpebart   = avg_cp_coxpebart,   
                               avg_cp_weibcoxph   = avg_cp_weibcoxph,   
                               avg_cp_coxpe       = avg_cp_coxpe)) )
}
