## Load Libraries and Data -----------------------------------------------------

source("load.R")

data("leuk_data", package = "BART4RS")
data("LeukSurv", package = "spBayesSurv")

log_sum_exp <- function(x) max(x) + log(sum(exp(x - max(x))))

## Global Variables ------------------------------------------------------------

NUM_SPLIT <- 10
NUM_FOLD  <- 5
SEED      <- digest::digest2int("CV-leukemia-4")
seeds     <- sample.int(.Machine$integer.max, NUM_SPLIT)

## Evaluation Function ---------------------------------------------------------

eval_split <- function(seed, fold) {
  set.seed(seed)
  folds <- caret::createFolds(seq_len(nrow(leuk_data)), k = NUM_FOLD)
  test_idx <- folds[[fold]]
  train_idx <- setdiff(seq_len(nrow(leuk_data)), test_idx)

  leuk_train <- leuk_data[train_idx,]
  leuk_test  <- leuk_data[test_idx,]

  ## Fit CoxBart ---------------------------------------------------------------
  fitted_coxpe <- BART4RS::coxpe_bart(
    formula = Surv(event_time, status) ~ age + sex + wbc + tpi,
    data = leuk_train,
    test_data = leuk_test,
    pop_haz = leuk_train$haz_rate,
    num_save = 4000,
    num_thin = 1,
    num_burn = 4000
  )

  ## Get Heldout Log-likelihood for CoxBart ------------------------------------
  formula <- Surv(event_time, status) ~ age + sex + wbc + tpi
  dv <- caret::dummyVars(formula, leuk_data)
  terms <- attr(dv$terms, "term.labels")
  group <- dummy_assign(dv)
  suppressWarnings({
    X_train <- predict(dv, leuk_train)
    X_test  <- predict(dv, leuk_test)
  })
  train_frame <- model.frame(formula, leuk_train)
  test_frame <- model.frame(formula, leuk_test)
  s_train <- model.response(train_frame)
  s_test  <- model.response(test_frame)
  Y_train <- s_train[,1]
  Y_test  <- s_test[,1]
  delta_train <- s_train[,2]
  delta_test  <- s_test[,2]

  obs_to_bin_test <- sapply(Y_test, find_bins, time_grid = fitted_coxpe$time_grid) - 1
  obs_to_bin_train <- sapply(Y_train, find_bins, time_grid = fitted_coxpe$time_grid) - 1

  test_loglik <- eval_loglik(y = Y_test,
                             delta = delta_test,
                             lambda = exp(fitted_coxpe$lambda_test),
                             grid = fitted_coxpe$time_grid,
                             pop_haz = leuk_test$haz_rate,
                             base_haz = fitted_coxpe$hazard,
                             cum_base_haz = fitted_coxpe$cum_hazard,
                             obs_to_bin = obs_to_bin_test)

  ## Fit Weibull ---------------------------------------------------------------
  weibform <- Surv(event_time, status) ~ ns(age, 5) + ns(log(1 + wbc), 5) + ns(tpi, 5) + sex
  dv <- caret::dummyVars(weibform, leuk_data)
  terms <- attr(dv$terms, "term.labels")
  group <- dummy_assign(dv)
  suppressWarnings({
    Xw_train <- predict(dv, leuk_train)
    Xw_test  <- predict(dv, leuk_test)
  })

  weibull_fit <- weibull_coxph(
    weibform,
    data = leuk_train,
    pop_haz = leuk_train$haz_rate,
    test_data = leuk_test,
    verbose = FALSE
  )

  ## Compute Weibull Heldout Log-likelihood ------------------------------------
  beta_samples  <- as.matrix(weibull_fit, "beta")
  shape         <- as.matrix(weibull_fit, "weibull_shape")
  ws            <- as.matrix(weibull_fit, "weibull_scale")
  ll            <- as.matrix(weibull_fit, "loglik")

  loglik_iter <- function(y, delta, X, beta, shape, scale, pop) {
    eta <- X %*% beta
    base <- (y / scale)^(shape - 1) * shape / scale
    cum <- (y / scale)^shape
    loglik <- delta * log(base * exp(eta) + pop) - cum * exp(eta)
    return(sum(loglik))
  }

  weibull_test_loglik <- sapply(
    1:4000,
    \(i) loglik_iter(
      y = Y_test,
      delta = delta_test,
      X = Xw_test,
      beta = beta_samples[i, ],
      shape = shape[i],
      scale = ws[i],
      pop = leuk_test$haz_rate
    )
  )

  ## Fit CoxLinear -------------------------------------------------------------
  formula_cox <- Surv(event_time, status) ~ ns(age, 5) + ns(wbc, 5) + ns(tpi, 5) + sex
  data_cox    <- leuk_data[train_idx,]
  pop_haz_cox <- leuk_data$haz_rate[train_idx]

  no_intercept_formula <- update(formula_cox, . ~ . - 1)
  train_frame     <- model.frame(no_intercept_formula, data_cox)
  X_train_cox     <- model.matrix(no_intercept_formula, data_cox)
  s_train_cox     <- model.response(train_frame)
  Y_train_cox     <- s_train[,1]
  delta_train_cox <- s_train[,2]
  grid_cox        <- fitted_coxpe$time_grid

  test_frame     <- model.frame(no_intercept_formula, leuk_data[test_idx,])
  X_test_cox     <- model.matrix(no_intercept_formula, leuk_data[test_idx,])

  stan_data <- list(
    X = X_train_cox,
    y = Y_train_cox,
    delta = delta_train_cox,
    N = nrow(X_train_cox),
    P = ncol(X_train_cox),
    grid_size = length(grid_cox) - 1,
    grid = grid_cox[-length(grid_cox)],
    pop_haz = leuk_data$haz_rate[train_idx]
  )

  fitted_cox_linear <- rstan::sampling(
    coxpe_stan,
    data = stan_data,
    chain = 4,
    cores = 4,
    iter = 2000
  )

  ## Getting Log-likelihood for CoxLinear --------------------------------------
  base_haz_cox     <- as.matrix(fitted_cox_linear, "lambda")
  cum_base_haz_cox <- as.matrix(fitted_cox_linear, "base_cum_haz")
  beta_cox         <- as.matrix(fitted_cox_linear, "beta")
  lambda_cox       <- beta_cox %*% t(X_test_cox)

  test_loglik_cox <- eval_loglik(
    y = Y_test,
    delta = delta_test,
    lambda = exp(lambda_cox),
    grid = fitted_coxpe$time_grid,
    pop_haz = leuk_test$haz_rate,
    base_haz = base_haz_cox,
    cum_base_haz = cum_base_haz_cox,
    obs_to_bin = obs_to_bin_test
  )

  log_evidence_bart_cox <- log_sum_exp(rowSums(test_loglik)) -
    log_sum_exp(rowSums(test_loglik_cox))
  log_evidence_bart_weib <- log_sum_exp(rowSums(test_loglik)) -
    log_sum_exp((weibull_test_loglik))

  out <- data.frame(
    seed = seed,
    fold = fold,
    log_evidence_bart_cox = log_evidence_bart_cox,
    log_evidence_bart_weib = log_evidence_bart_weib
  )

  return(out)
}

## Run Experiment or Load Results ----------------------------------------------

sim_settings <- expand.grid(seed = seeds, fold = 1:NUM_FOLD) %>% arrange(seed)
tmpf <- function(i) eval_split(sim_settings$seed[i], sim_settings$fold[i])

filename <- "cache/03_CV.rds"
if (!file.exists(filename)) {
  coxpe_stan <- stan_model("lib/PEHaz.stan")
  results <- parallel::mclapply(
    1:nrow(sim_settings),
    tmpf,
    mc.preschedule = FALSE,
    mc.cores = 5
  )
  results_df <- results %>% do.call(rbind, .)
  results_df %>% group_by(seed) %>% summarise_all(sum)
  saveRDS(results_df, file = filename)
}

results_df <- readRDS(file = filename)

## Summarize and Display Results -----------------------------------------------

results_sum <- results_df %>% group_by(seed) %>% summarise_all(sum) %>%
  rename(`Cox Linear` = log_evidence_bart_cox,
         Weibull = log_evidence_bart_weib) %>%
  pivot_longer(cols = c(`Cox Linear`, Weibull),
               names_to = "Model",
               values_to = "Log-Loss")

results_sum %>%
  group_by(Model) %>%
  summarise_all(mean) %>%
  mutate(Deviance = -2 * `Log-Loss`)
