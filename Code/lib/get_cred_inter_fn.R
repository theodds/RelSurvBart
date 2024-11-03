get_cis <- function(fits = fits, true_eta_centered   = eta_centered, data = data_gen){
  
  fit_coxpe_bart     <- fits[["fit_coxpe_bart"]]
  fit_weibull_coxph  <- fits[["fit_weibull_coxph"]]
  fitted_pe          <- fits[["fitted_pe"]]
  
  ## Get CI & CP for coxpe_bart --------------------------------
  eta_hat_cent_coxpebart_post  <- fit_coxpe_bart$lambda_train - rowMeans(fit_coxpe_bart$lambda_train)
  ci_etahat_centered_coxpebart <- matrix(unlist(apply(eta_hat_cent_coxpebart_post, 2, 
                                                      function(x) bayestestR::ci(x, ci=0.9, method="HDI"))), 
                                         nrow=1043, byrow = T)
  colnames(ci_etahat_centered_coxpebart) <- c("ci_level", "lci_coxpebart", "uci_coxpebart")
  cov_prob_coxpebart           <- sum(true_eta_centered > ci_etahat_centered_coxpebart[,2] & 
                                      true_eta_centered < ci_etahat_centered_coxpebart[,3]) / 1043
  
  
  
  ## Get CI & CP for WeibcoxPH --------------------------------
  beta_post_weibcoxph          <- as.matrix(fit_weibull_coxph, "beta")
  eta_hat_post_weibcoxph       <- beta_post_weibcoxph %*% t(cbind(data$age, data$sex, data$wbc, data$tpi))
  eta_hat_cent_post_weibcoxph  <- eta_hat_post_weibcoxph - rowMeans(eta_hat_post_weibcoxph)
  ci_etahat_centered_weibcoxph <- matrix(unlist(apply(eta_hat_cent_post_weibcoxph, 2, 
                                                      function(x) bayestestR::ci(x, ci=0.9, method="HDI"))),
                                         nrow = 1043, byrow = T)
  colnames(ci_etahat_centered_weibcoxph) <- c("ci_level", "lci_weibcoxph", "uci_weibcoxph")
  cov_prob_weibcoxph           <- sum(true_eta_centered > ci_etahat_centered_weibcoxph[,2] & 
                                        true_eta_centered < ci_etahat_centered_weibcoxph[,3]) / 1043
  
  
  ## Get CI & CP for CoxPE --------------------------------
  beta_post_coxpe              <- as.matrix(fitted_pe, "beta")
  eta_hat_post_coxpe           <- beta_post_coxpe %*% t(cbind(data$age, data$sex, data$wbc, data$tpi))
  eta_hat_cent_post_coxpe      <- eta_hat_post_coxpe - rowMeans(eta_hat_post_coxpe)
  ci_etahat_centered_coxpe     <- matrix(unlist(apply(eta_hat_cent_post_coxpe, 2, 
                                                      function(x) bayestestR::ci(x, ci=0.9, method="HDI"))),
                                         nrow = 1043, byrow = T)
  colnames(ci_etahat_centered_coxpe) <- c("ci_level", "lci_coxpe", "uci_coxpe")
  cov_prob_coxpe               <- sum(true_eta_centered > ci_etahat_centered_coxpe[,2] & 
                                        true_eta_centered < ci_etahat_centered_coxpe[,3]) / 1043
  
  
  return(cp = list(cov_prob_coxpebart = cov_prob_coxpebart,
                   cov_prob_weibcoxph = cov_prob_weibcoxph,
                   cov_prob_coxpe     = cov_prob_coxpe,
                   ci_etahat_centered_coxpebart = ci_etahat_centered_coxpebart[,2:3],
                   ci_etahat_centered_weibcoxph = ci_etahat_centered_weibcoxph[,2:3],
                   ci_etahat_centered_coxpe     = ci_etahat_centered_coxpe[,2:3]))
}