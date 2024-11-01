#' Additive Non-Proportional Hazards Model with Bayesian Additive Regression Trees (BART)
#'
#' This function fits a non-proportional hazards regression model with a BART
#' component. The hazard function is modeled as piecewise constant on time
#' intervals \eqn{0 = t_0 < t_1 < \cdots < t_B = \infty}, and is given by
#'
#' \deqn{h(t \mid X) = h_P(t \mid X) + h_0(t) e^{g(X,t)}}
#'
#' where:
#' \itemize{
#'  \item \eqn{h(t | X)}: is the hazard function at time \eqn{t} given covariate vector \eqn{X}.
#'  \item \eqn{h_P(t | x)}: is the population hazard rate, which is assumed known.
#'  \item \eqn{h_0(t)}: is the baseline hazard function at time \eqn{t} of the excess hazard, modeled as a piecewise constant function.
#'  \item \eqn{X}: is the covariate vector.
#'  \item \eqn{g(X,t)}: is the non-linear component estimated using BART, which is also assumed to be piecewise-constant as a function of $t$.
#' }
#'
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
#' @return Returns a fitted non-proportional hazards model object of class "coxnph_bart". The object is a list containing the following components.
#' \itemize{
#'   \item \code{lambda_test}: A 3-dimensional array where lambda_test[i,j,k] represents the value of \eqn{g(X_i,T_j)} at iteration \eqn{k} of the chain for observations in the test set.
#'   \item \code{lambda_train}: A 3-dimensional array where lambda_test[i,j,k] represents the value of \eqn{g(X_i,T_j)} at iteration \eqn{k} of the chain for observations in the training set.
#'   \item \code{hazard}: A matrix where hazard[i,j] is the baseline hazard at iteration \eqn{i} of the chain in bin \eqn{j}.
#'   \item \code{counts}: Matrix consisting of the number of times that predictor \eqn{j} was used at iteration \eqn{i}.
#'   \item \code{loglik_obs}: Individual level log-likelihoods, such that loglik_obs[i,j] is the log-likelihood of individual \eqn{j} at iteration \eqn{i} of the chain.
#'   \item \code{time_grid}: The grid of times used by the NPH model.
#' }
#'
#'
#' @examples
#' \dontrun{
#' data("leuk_data")
#' }
coxnph_bart <- function(formula,
                        data,
                        pop_haz = NULL,
                        test_data = NULL,
                        debug_ph = FALSE,
                        num_tree = 50,
                        k = 1,
                        time_grid = NULL,
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

  stopifnot(is.numeric(Y_train))

  ## Set up hypers
  scale_lambda <- 1 / sqrt(num_tree) / k

  ## Preprocess Data
  ab <- auto_bin(Y_train, time_grid = time_grid)
  bin_to_obs_list <- ab$bin_to_obs
  obs_to_bin      <- ab$obs_to_bin
  time_grid       <- ab$time_grid
  bin_width       <- ab$bin_width
  base_haz_init   <- rgamma(ab$num_bins, 1)
  do_rel_surv     <- ifelse(is.null(pop_haz), FALSE, TRUE)

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
  X_train <- cbind(X_train, 0)
  X_test <- cbind(X_test, 0)
  probs           <- Matrix::Matrix(diag(ncol(X_train)))

  if(debug_ph) {
    probs[,ncol(probs)] <- c(rep(1/(ncol(X_train) - 1), ncol(X_train) - 1), 0)
  }

  fitted_nph <-
    CoxNPHBart(
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

  ## Compute the log-likelihood for model evaluation ----

  eval_loglik_nph <- function(Z, r, i, delta, bin, haz, pop_haz) {
    bin <- bin[i]
    A <- delta[i] * log(pop_haz[i] + r[i,bin,] * haz[,bin])
    num_bin <- length(haz[1,])
    B <- Z[i,1] * r[i,1,] * haz[,1]
    for(b in 2:(num_bin)) {
      B <- B + Z[i,b] * r[i,b,] * haz[,b]
    }
    return(A - B)
  }

  bin <- obs_to_bin
  haz <- fitted_nph$hazard
  r <- exp(fitted_nph$lambda_train)

  Z <- matrix(nrow = length(Y_train), ncol = ncol(haz))
  for(i in 1:length(Y_train)) {
    for(j in 1:ncol(haz)) {
      if(Y_train[i] < time_grid[j]) {
        Z[i,j] <- 0
      }
      if(Y_train[i] >= time_grid[j]) {
        if(Y_train[i] < time_grid[j + 1]) {
          Z[i,j] <- Y_train[i] - time_grid[j]
        } else {
          Z[i,j] <- time_grid[j+1] - time_grid[j]
        }
      }
    }
  }

  out <- sapply(1:length(Y_train),
                \(i) eval_loglik_nph(Z = Z, r = r, i = i,
                                     delta = delta_train, bin = bin + 1,
                                     haz = haz, pop_haz =pop_haz))
  
  fitted_nph$loglik_obs <- out
  fitted_nph$time_grid <- time_grid

  return(fitted_nph)
}

extract_R <- function(fit) {
  get_j <- function(j) {
    exp(-fit$hazard[,j] * (fit$time_grid[j+1] - fit$time_grid[j]) *
        t(exp(fit$lambda_train[,j,])))
  }
  J <- length(fit$time_grid) - 1
  out <- lapply(1:J, get_j)
  return(out)
}
