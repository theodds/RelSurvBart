expand_x_test <- function(x, t) {
  t_test <- rep(t, nrow(x))
  n_t    <- length(t)
  f <- function(i) {
    matrix(rep(x[i,], n_t), nrow = n_t, byrow = TRUE)
  }
  x_test <- do.call(rbind, lapply(1:nrow(x), f))
  return(list(X_test = x_test, t_test = t_test))
}


update_trees_weib <- function(my_forest, alpha, t, X, t_obs, X_obs, idx_obs,
                              rate, shape, X_test, t_test)
{

  ## Augment
  out            <- WeibAugment(t, X, rate, shape)
  t_reject       <- out[["rejected_times"]]
  subject_id     <- out[["subject_id"]]
  X_reject       <- out[["X_reject"]]

  ## Do the thinning and get mu_hat for thinned
  mu_hat_reject  <- alpha + my_forest$predict(X_reject, t_reject)
  p_accept       <- pnorm(mu_hat_reject)
  N_reject       <- length(t_reject)
  accept_idx     <- ifelse(runif(N_reject) < p_accept, TRUE, FALSE)
  X_thin         <- X_reject[!accept_idx, ]
  t_thin         <- t_reject[!accept_idx]
  idx_thin       <- subject_id[!accept_idx]
  mu_hat_thin    <- mu_hat_reject[!accept_idx] 

  ## Get mu_hat for observed
  mu_hat_obs <- alpha + my_forest$predict(X_obs, t_obs)

  ## Get Z
  Z_obs <- rtruncnorm(n = length(mu_hat_obs),
                      mean = mu_hat_obs,
                      a = 0, b = Inf)
  if(length(mu_hat_thin) > 0) {
    Z_thinned <- rtruncnorm(n = length(mu_hat_thin),
                            mean = mu_hat_thin,
                            a = -Inf, b = 0)
  }
  else {
    Z_thinned <- numeric(0)
  }

  ## Update the trees
  X_aug <- rbind(X_thin, X_obs)
  Z_aug <- c(Z_thinned, Z_obs)
  t_aug <- c(t_thin, t_obs)
  idx_aug <- c(idx_thin, idx_obs)
  mu_hat_pred <- my_forest$do_gibbs(X_aug, Z_aug - alpha,
                                    t_aug, X_aug, t_aug, 1)

  R <- Z_aug - mu_hat_pred
  alpha <- rtruncnorm(1, 0, Inf,
                      sum(R) / (length(R) + 1), 1 / sqrt(length(R) + 1))
  ## alpha <- rnorm(1, (1 + sum(R)) / (length(R) + 1),
  ##                1 / sqrt(length(R)+1))
  return(list(X_aug = X_aug,
              t_aug = t_aug,
              mu_hat = alpha + mu_hat_pred,
              idx_aug = idx_aug,
              alpha = alpha
              ))

}


modsurv <- function(X, y, delta, X_test, t_test,
                    hypers_base, opts_base,
                    hypers_thin, opts_thin,
                    num_thin, num_burn, num_save)
{

  ## Initialize forest
  X_obs <- X[delta == 1, ]
  t_obs <- y[delta == 1]
  t_mis <- y[delta == 0]
  idx   <- 1:nrow(X)
  idx   <- idx[delta == 1] - 1
  tmp_test       <- expand_x_test(X_test, t_test)
  tX             <- X_test
  tt             <- t_test
  X_test         <- tmp_test$X_test
  t_test         <- tmp_test$t_test
  haz_est <- matrix(data = NA, nrow = num_save, ncol = length(t_test))
  alpha_probit_samps <- numeric(num_save)
  augment_size <- numeric(num_save)
  params_thin <- vector("list", num_save)
  params_base <- vector("list", num_save)
  
  weib_forest <- with(hypers_base, MakeWeib(probs = probs, 
                                            num_trees = num_trees, 
                                            scale_lambda = scale_lambda, 
                                            shape_lambda_0 = shape_lambda_0, 
                                            rate_lambda_0 = rate_lambda_0, 
                                            weibull_power = weibull_power,
                                            update_scale = update_scale))
  
  mod_forest   <- ModBart::MakeForest(hypers_thin, opts_thin)
  alpha_probit <- 1
  
  if(opts_base$do_ard) weib_forest$do_ard()
  
  for(i in 1:1000) {
    weib_forest$do_gibbs(X,y,t_obs,idx, X, 1)
    cat(paste("\rFinishing pre-warmup iteration", i, "\t\t\t"))
  }
  cat("\n")
  alpha_probit <- 1

  for(i in 1:num_burn) {
    rate <- exp(weib_forest$predict(X))
    tmp <- update_trees_weib(mod_forest, alpha_probit, y, X, t_obs, X_obs, idx,
                             rate, hypers_base$weibull_power, X_test, t_test)
    alpha_probit <- tmp$alpha
    tmp2 <- weib_forest$do_gibbs(X, y, tmp$t_aug, tmp$idx_aug, X_test, 1)
    
    if(i %% 100 == 0)
      cat(paste("\rFinishing warmup", i, "Augment size = ", length(tmp$t_aug), "\t\t\t"))
  }
  cat("\n")
  
  for(i in 1:num_save) {
    for(j in 1:num_thin) {
      rate <- exp(weib_forest$predict(X))
      tmp <- update_trees_weib(mod_forest, alpha_probit, y, X, t_obs, X_obs, idx,
                               rate, hypers_base$weibull_power, X_test, t_test)
      alpha_probit <- tmp$alpha
      tmp2 <- weib_forest$do_gibbs(X, y, tmp$t_aug, tmp$idx_aug, X, 1)
    }
    
    rate_test <- exp(weib_forest$predict(X_test))
    base_rate <- rate_test * hypers_base$weibull_power * t_test^(hypers_base$weibull_power - 1)
    thin_rate <- alpha_probit + mod_forest$predict(X_test, t_test)
    haz_test  <- base_rate * pnorm(thin_rate)
    haz_est[i,] <- haz_test
    alpha_probit_samps[i] <- alpha_probit
    augment_size[i] <- length(tmp$idx_aug)
    params_thin[[i]] <- mod_forest$get_params()
    params_base[[i]] <- weib_forest$get_params()

    if(i %% 100 == 0)
      cat(paste("\rFinishing save", i,  "Augment size = ", length(tmp$t_aug), "\t\t\t"))
    
    
  }
  cat("\n")

  Delta <- t_test[2] - t_test[1]
  get_surv <- function(i) {
    h <- matrix(haz_est[i,], nrow = nrow(tX), ncol = length(tt), byrow = TRUE)
    s_est <- apply(h, 1, function(x) {
      c(1, exp(-cumsum(x) * Delta))[-(length(tt)+1)]
    })
    return(s_est)
  }

  surv_samps <- array(sapply(1:num_save, get_surv),
                      c(length(tt),nrow(tX),num_save))
  surv_hat <- t(apply(surv_samps,c(1,2),mean))

  return(list(haz_test = haz_est, 
              # weib_forest = weib_forest, 
              # mod_forest = mod_forest,
              surv_samps = surv_samps,
              surv_hat = surv_hat,
              alpha_probit = alpha_probit_samps,
              augment_size = augment_size,
              params_base = params_base,
              params_thin = params_thin
              ))
  
}

