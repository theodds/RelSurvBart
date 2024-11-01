#' Additive Hazards Cox Proportional Hazards Model with Bayesian Additive Regression Trees (BART)
#'
#' This function fits a Cox proportional hazards regression models with a BART component.
#' The baseline hazard is modeled using a piecewise constant hazard function to obtain a
#' flexible model for the baseline hazard. The model is given by
#'
#' \deqn{h(t | X) = h_P(t | X) + h_0(t) e^{g(X)}}
#'
#' where:
#' \itemize{
#'  \item \eqn{h(t | X)}: is the hazard function at time \eqn{t} given covariate vector \eqn{X}.
#'  \item \eqn{h_P(t | x)}: is the population hazard rate, which is assumed known.
#'  \item \eqn{h_0(t)}: is the baseline hazard function at time \eqn{t} of the excess hazard, modeled as a piecewise constant function.
#'  \item \eqn{X}: is the covariate vector.
#'  \item \eqn{g(X)}: is the non-linear component estimated using BART.
#' }
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#'   a symbolic description of the model to be fitted. The response should be a `Surv` object
#'   as created by the `survival::Surv` function.
#' @param data A data frame, list, or environment (or object coercible by `as.data.frame`
#'   to a data frame) containing the variables in the model.
#' @param pop_haz (optional) A numeric vector specifying the population hazard. If `NULL`,
#'   a relative survival model is not fit.
#' @param test_data (optional) A data frame for the test set. If `NULL`, the training data
#'   is used for testing.
#' @param num_tree The number of trees to be used in the BART model. Defaults to 50.
#' @param k A parameter controlling the amount of regularization in the BART model.
#'   Defaults to 1.
#' @param num_burn The number of MCMC iterations to discard as burn-in. Defaults to 2500.
#' @param num_thin The thinning parameter for the MCMC sampler. Defaults to 1.
#' @param num_save The number of MCMC samples to save. Defaults to 2500.
#'
#' @return Returns a fitted Cox proportional hazards model object with BART, of class `"coxpe_bart"`. The object is a list containing the following components:
#' \itemize{
#' \item \code{r_test}: A numeric matrix representing the samples of \eqn{g(X)} for each observation in the test set, which has been centered to have mean 0 over the training set.
#' \item \code{r_train}: A numeric matrix representing the samples of \eqn{g(X)} for each observation in the training set, which has been centered to have mean 0 over the training set.
#' \item \code{base_haz}: A numeric matrix representing the samples of the parameters for the piecewise constant hazard.
#' \item \code{lambda_test}: A numeric matrix representing the samples of \eqn{g(X)} for each observation in the test set, which has not been centered.
#' \item \code{lambda_train}: A numeric matrix representing the samples of \eqn{g(X)} for each observation in the training set, which has not been.
#' \item \code{base_haz}: A numeric matrix representing the samples of the parameters for the piecewise constant hazard, which accounts for centering of \eqn{g(X)}.
#' \item \code{hazard}: A numeric matrix representing the samples of the parameters for the piecewise constant hazard, which does not account for centering of \eqn{g(X)}.
#' \item \code{cum_hazard}: A numeric matrix representing the samples of the parameters for the piecewise constant cumulative hazard.
#' \item \code{counts}: A matrix containing the count of the number of times each predictor is used in the trees of the BART model, across all saved MCMC samples.
#' \item \code{sigma_lambda}: A numeric vector containing the posterior samples of the leaf parameter standard deviation for the BART component of the model.
#' \item \code{loglik}: A numeric vector representing the log-likelihood of the model at each saved MCMC iteration.
#' \item \code{time_grid}: grid representing the bins along which the hazard is piecewise constant.
#' \item \code{loglik_obs}: matrix containing the unit-level log-likelihood evaluations, which can be used with the loo function in the loo package to compute leave-one-out expected log predictive densities for model comparison.
#' }
#'
#'
#' @examples
#' \dontrun{
#' data("leuk_data")
#' coxpe_bart(formula = Surv(event_time, status) ~ age + sex,
#'            data = leuk_data,
#'            pop_haz = leuk_data$haz_rate)
#' }
coxpe_bart <- function(formula,
                       data,
                       pop_haz = NULL,
                       test_data = NULL,
                       time_grid = NULL,
                       num_tree = 50,
                       k = 1,
                       num_burn = 2500,
                       num_thin = 1,
                       num_save = 2500) {

  if(is.null(test_data)) test_data <- data

  dv <- dummyVars(formula, data)
  terms <- attr(dv$terms, "term.labels")
  group <- dummy_assign(dv)
  suppressWarnings({
    X_train <- predict(dv, data)
    X_test  <- predict(dv, test_data)
  })
  train_frame <- model.frame(formula, data)
  test_frame <- model.frame(formula, test_data)
  s_train <- model.response(train_frame)
  s_test  <- model.response(test_frame)
  Y_train <- s_train[,1]
  Y_test  <- s_test[,1]
  delta_train <- s_train[,2]
  delta_test  <- s_test[,2]

  if(is.null(pop_haz)) pop_haz <- rep(0, length(Y_train))

  stopifnot(is.numeric(Y_train))

  ## Set up hypers
  scale_lambda <- 1 / sqrt(num_tree) / k

  ## Preprocess Data
  ab              <- auto_bin(Y_train, time_grid = time_grid)
  bin_to_obs_list <- ab$bin_to_obs
  obs_to_bin      <- ab$obs_to_bin
  time_grid       <- ab$time_grid
  bin_width       <- ab$bin_width
  base_haz_init   <- rgamma(ab$num_bins, 1)
  do_rel_surv     <- ifelse(is.null(pop_haz), FALSE, TRUE)
  probs           <- Matrix::Matrix(diag(length(terms)))

  ## Normalize!

  make_01_norm <- function(x) {
    a <- min(x)
    b <- max(x)
    return(function(y) (y - a) / (b - a))
  }

  ecdfs   <- list()
  for(i in seq_len(ncol(X_train))) {
    ecdfs[[i]] <- ecdf(X_train[,i])
    if(length(unique(X_train[,i])) == 1) ecdfs[[i]] <- identity
    if(length(unique(X_train[,i])) == 2) ecdfs[[i]] <- make_01_norm(X_train[,i])
  }
  for(i in seq_len(ncol(X_train))) {
    X_train[,i] <- ecdfs[[i]](X_train[,i])
    X_test[,i] <- ecdfs[[i]](X_test[,i])
  }


  ## Fit the model ----

  fitted_cox <-
    CoxPEBart(
      X_train,
      Y_train,
      delta_train,
      bin_to_obs_list,
      obs_to_bin,
      time_grid,
      bin_width,
      base_haz_init,
      probs,
      X_test,
      num_tree,
      scale_lambda,
      do_rel_surv,
      pop_haz,
      num_burn,
      num_thin,
      num_save
    )

  r     <- fitted_cox$lambda_train - rowMeans(fitted_cox$lambda_train)
  r_hat <- colMeans(r)
  fitted_cox$r_train      <- r
  fitted_cox$r_train_hat  <- r_hat

  r     <- fitted_cox$lambda_test - rowMeans(fitted_cox$lambda_train)
  r_hat <- colMeans(r)
  fitted_cox$r_test      <- r
  fitted_cox$r_test_hat  <- r_hat

  class(fitted_cox)    <- "coxpe_bart"
  fitted_cox$time_grid <- time_grid

  fitted_cox$base_haz <- fitted_cox$hazard * exp(rowMeans(fitted_cox$lambda_train))

  fitted_cox$loglik_obs <- eval_loglik(
    y = Y_train,
    delta = delta_train,
    lambda = exp(fitted_cox$lambda_train),
    grid = fitted_cox$time_grid,
    base_haz = fitted_cox$hazard,
    cum_base_haz = fitted_cox$cum_hazard,
    pop_haz = pop_haz,
    obs_to_bin = obs_to_bin
  )

  return(fitted_cox)
}
