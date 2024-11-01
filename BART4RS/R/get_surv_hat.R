get_surv_hat <- function(j, time_grid, surv_grid) {
  base <- fitted_cox$cum_hazard[j,] * exp(mean(fitted_cox$lambda_train[j,]))
  Lambda_hat <- approx(x = time_grid, y = c(0, base), xout = surv_grid)
  surv_hat <- exp(-Lambda_hat$y)
  return(surv_hat)
}

get_surv_samples <- function(fit, surv_grid = NULL) {

  if(is.null(surv_grid)) {
    surv_grid <- seq(from = 0, to = max(fit$time_grid), length = 400)
  }

  out <-
    sapply(
      seq_len(length(fit$sigma_lambda)),
      get_surv_hat,
      time_grid = fit$time_grid,
      surv_grid = surv_grid
    )

  return(list(samples = t(out), grid = surv_grid))

}
