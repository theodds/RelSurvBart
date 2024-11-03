weibsurv <- function(X, y, delta, X_test, t_test,
                     hypers, opts, num_thin, num_burn, num_save)
{

  ## Initialize forest
  X_obs <- X[delta == 1, ]
  t_obs <- y[delta == 1]
  t_mis <- y[delta == 0]
  idx   <- 1:nrow(X)
  idx   <- idx[delta == 1] - 1
  haz_est <- matrix(data = NA, nrow = num_save, ncol = nrow(X_test))
  params_weib <- vector("list", num_save)
  surv_samps <- array(NA, c(nrow(X_test), length(t_test), num_save))
  
  weib_forest <- MakeWeib(
    probs = hypers$probs,
    num_trees = hypers$num_trees,
    scale_lambda = hypers$scale_lambda,
    shape_lambda_0 = hypers$shape_lambda_0,
    rate_lambda_0 = hypers$rate_lambda_0,
    weibull_power = hypers$weibull_power,
    update_scale = hypers$update_scale
  )


  if(opts$do_ard) weib_forest$do_ard()

  for(i in 1:num_burn) {
    weib_forest$do_gibbs(X,y,t_obs,idx, X_test, 1)
    if(i %% 100 == 0)
      cat(paste("\rFinishing warmup iteration", i, "\t\t\t"))
  }
  cat("\n")

  for(i in 1:num_save) {
    for(j in 1:num_thin) {
      weib_forest$do_gibbs(X, y, t_obs, idx, X_test, 1)
    }
    haz_est[i,] <- exp(weib_forest$predict(X_test))
    params_weib[[i]] <- weib_forest$get_params()
    
    for(k in 1:nrow(X_test)) {
      surv_samps[k,,i] <- exp(-haz_est[i,k] * t_test^(hypers$weibull_power))
    }
    
    if(i %% 100 == 0)
      cat(paste("\rFinishing save", i, "\t\t\t"))
  }
  
  cat("\n")
  
  surv_hat <- apply(surv_samps,c(1,2),mean)
  
  return(list(
    weib_rate_test = haz_est,
    surv_samps = surv_samps,
    surv_hat = surv_hat,
    params_weib = params_weib
  ))
}
