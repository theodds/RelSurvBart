#' Weibull Cox Proportional Hazards Regression Model
#'
#' This function fits an additive/excess hazard regression model to survival
#' data using a Cox proportional hazards model with a Weibull baseline hazard.
#' The model can be mathematically formulated as:
#'
#' \deqn{h(t | X) = h_P(t | X) + h_E(t | X)}
#'
#' where:
#' \itemize{
#' \item \eqn{h(t | X)}: is the hazard function at time \eqn{t} given covariate vector \eqn{X}.
#' \item \eqn{h_P(t | X)}: is the population hazard, which is assumed known.
#' \item \eqn{h_E(t | X)}: is the excess hazard, which follows a Cox proportional hazards model with Weibull baseline hazard \eqn{h_0(t | X) \, e^{X^\top\beta}}.
#' \item \eqn{X}: is the covariate vector.
#' \item \eqn{\beta}: is the vector of coefficients.
#' }
#'
#' @param formula An object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted.
#' @param data A data frame, list, or environment (or object coercible by
#' as.data.frame to a data frame) containing the variables in the model.
#' @param pop_haz A numeric vector specifying the population hazard. If NULL,
#' a vector of zeros with length equal to the number of rows in the training
#' data is used, indicating that a standard survival model will be fit.
#' @param test_data A data frame for the test set. If NULL, the training data
#' is used for testing.
#' @param verbose A logical value. If TRUE, it prints out additional
#' information including a warning message for unimplemented functionality
#' related to test data evaluations.
#' @param ... Additional arguments to be passed to the `rstan::sampling`
#' function.
#'
#' @return Returns a fitted Weibull regression model object from the `rstan`
#' package.
#'
#' @seealso
#' \code{\link[rstan]{sampling}}
#'
#' @examples
#' \dontrun{
#' data("leuk_data")
#' weibull_coxph(formula = Surv(event_time, status) ~ age + sex,
#'               pop_haz = leuk_data$haz_rate,
#'               data = leuk_data)
#' }
weibull_coxph <- function(formula,
                          data,
                          pop_haz = NULL,
                          test_data = NULL,
                          verbose = TRUE,
                          ...) {

  ## Preprocess data

  if(is.null(test_data)) test_data <- data

  no_intercept_formula <- update(formula, . ~ . - 1)

  train_frame <- model.frame(no_intercept_formula, data)
  test_frame  <- model.frame(no_intercept_formula, test_data)

  X_train <- model.matrix(no_intercept_formula, data)
  X_test  <- model.matrix(no_intercept_formula, test_data)
  s_train <- model.response(train_frame)
  s_test  <- model.response(test_frame)
  Y_train <- s_train[,1]
  Y_test  <- s_test[,1]
  delta_train <- s_train[,2]
  delta_test  <- s_test[,2]

  if(is.null(pop_haz)) pop_haz <- rep(0, nrow(X_train))

  ## Make stan list

  stan_data <- list(X = X_train,
                    y = Y_train,
                    delta = delta_train,
                    N = nrow(X_train),
                    P = ncol(X_train),
                    pop_haz = pop_haz)

  ## Fitted model

  fitted_weibull <- rstan::sampling(BART4RS::stanmodels$WeibHaz,
                                    data = stan_data,
                                    chain = 4,
                                    cores = 4,
                                    ...)

  ## Return

  if(verbose) {
    warning("TODO: add functionality for returning evaluations on test data")
  }

  return(fitted_weibull)

}
