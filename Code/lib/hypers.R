WeibHypers <- function(X, Y,
                       update_scale_mu = FALSE,
                       do_ard = TRUE,
                       num_trees = 50,
                       num_burn = 5000,
                       num_thin = 1,
                       num_save = 5000) {
  hypers <- list()
  opts <- LVBart::Opts(num_burn = num_burn,
                       num_thin = num_thin,
                       num_save = num_save)

  opts$do_ard = do_ard
  hypers$probs     <- Matrix::Matrix(diag(ncol(X)))
  hypers$num_trees <- num_trees
  hypers$scale_lambda <- 1.5 * sd(log(Y)) / sqrt(num_trees)
  hypers$shape_lambda_0 <- 0.5
  hypers$rate_lambda_0 <- 0.5
  hypers$weibull_power <- 1
  hypers$update_scale <- update_scale_mu

  return(list(hypers = hypers, opts = opts))
  
}

ThinHypers <- function(X, Y, 
                       update_scale_mu = TRUE,
                       do_ard = TRUE,
                       num_burn = 5000, num_thin = 1, num_save = 5000,
                       ...) {


  opts <- LVBart::Opts(num_burn = num_burn,
                       num_thin = num_thin,
                       num_save = num_save)
  hypers <- LVBart::Hypers(X = X, Y = Y, ...)
  hypers$sigma_mu_hat <- 3 / 2/ sqrt(hypers$num_tree)
  hypers$k            <- 2
  hypers$sigma_hat    <- 1
  hypers$length_scale <- max(Y) / pi
  hypers$shape_length_scale <- 1
  hypers$rate_length_scale <- 1 / hypers$length_scale^2
  opts$update_sigma <- FALSE
  opts$update_sigma_mu <- update_scale_mu
  opts$do_ard <- do_ard

  return(list(hypers = hypers, opts = opts))

}

fit_coxbart <- function(X, Y, delta,
                        X_test = NULL,
                        k_lambda = 2,
                        num_trees = 50,
                        num_burn = 5000,
                        num_thin = 1,
                        num_save = 5000) {

  
  process_surv <- function(Y, delta, X) 
  {
    o     <- order(Y, 1-delta)
    Y     <- Y[o]
    delta <- delta[o]
    X     <- X[o,]
    survs <- Surv(Y, delta)
    k <- 1
    U <- numeric(length(unique(survs)))
    L <- numeric(length(unique(survs)))
    L[1] <- 1
    for(i in 2:length(survs))
    {
      if(!identical(survs[i], survs[i-1])) {
        U[k] <- i - 1
        k <- k + 1
        L[k] <- i
      }
    }
    U[length(unique(survs))] <- length(survs)
    return(list(o = 1:length(Y) - 1, Y = Y, L = L-1, U = U-1,
                delta = delta, X = X))
  }
  

  c(o,Y,L,U,delta,X) %<-% process_surv(Y = Y, delta = delta, X = X)
  if(is.null(X_test)) X_test <- X

  sigma_lambda <- 3 / k_lambda / sqrt(num_trees)
  probs <- Matrix::Matrix(diag(ncol(X)))

  fitted_coxbart <- CoxBart(X, Y, delta, o, L, U, probs, X_test, num_trees,
                            sigma_lambda, num_burn, num_thin, num_save)
  
  get_surv <- function(hazard, log_rates, iteration, times) {
    baseline_cumulative_hazard <- c(0,cumsum(hazard))
    baseline_survival <- exp(-baseline_cumulative_hazard)
    rate <- exp(log_rates)
    surv <- sapply(rate, function(x) baseline_survival^x)
    # N <- length(log_rates)
    # N_rel <- length(baseline_survival)
    # out_df <- tibble(survival = as.numeric(surv),
    #                  idx = rep(1:N, each = N_rel),
    #                  iteration = iteration,
    #                  times = rep(c(0,times), N))
    # return(out_df)
    return(surv)
  }

  get_surv_iteration <- function(i) {
    hazard <- fitted_coxbart$hazard[i,]
    log_rates <- fitted_coxbart$lambda_test[i,]
    get_surv(hazard, log_rates, i, sort(Y))
  }
    
  surv_mean <- get_surv_iteration(1)
  for(i in 2:num_save) {
    surv_mean <- surv_mean + get_surv_iteration(i)
  }
  surv_mean <- surv_mean / num_save
  
  fitted_coxbart$survival_mean <- t(surv_mean)
  
  # fitted_coxbart$survival_df <- 
  #   do.call(rbind, lapply(1:num_save, get_surv_iteration))
  # fitted_coxbart$survival_df_mean <- 
  #   fitted_coxbart$survival_df %>% group_by(times, idx) %>% 
  #   summarise(survival = mean(survival))


  return(fitted_coxbart)

}










