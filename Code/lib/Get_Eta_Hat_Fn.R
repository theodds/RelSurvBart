get_eta_hats <- function(fits = fits, data = data_gen){
      
    fit_coxpe_bart     <- fits[["fit_coxpe_bart"]]
    fit_weibull_coxph  <- fits[["fit_weibull_coxph"]]
    fitted_pe          <- fits[["fitted_pe"]]
    
    beta_weibull_coxph <- summary(fit_weibull_coxph)$summary[, "mean"][1:4]
    
    
    
    ## Obtain centered eta_hat for fitted models -----------------------
    
    eta_hat_cenetered_coxpebart <- colMeans(fit_coxpe_bart$lambda_train - rowMeans(fit_coxpe_bart$lambda_train))
    
    eta_hat_weibcoxph            <- cbind(data$age, data$sex, data$wbc, data$tpi) %*% beta_weibull_coxph
    eta_hat_cenetered_weibcoxph  <- eta_hat_weibcoxph - mean(eta_hat_weibcoxph)
    
    eta_hat_coxpe                <- summary(fitted_pe, pars = "eta")$summary[, "mean"]
    eta_hat_cenetered_coxpe      <- eta_hat_coxpe - mean(eta_hat_coxpe)
    
    
    return(eta_hats = list(eta_hat_cenetered_coxpebart  = eta_hat_cenetered_coxpebart,
                           eta_hat_cenetered_weibcoxph  = eta_hat_cenetered_weibcoxph,
                           eta_hat_cenetered_coxpe      = eta_hat_cenetered_coxpe))
}