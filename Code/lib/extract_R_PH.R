extract_R_PH <- function(fit = fit_coxpe_bart){
  get_j      <- function(j){
                exp(-fit$hazard[, j] * (fit$time_grid[j + 1] - fit$time_grid[j]) * 
                 t(exp(fit$lambda_train)))
  }
  J          <- length(fit$time_grid) - 1
  out        <- lapply(1:J, get_j)
  return(out)
  
}

extract_S_PH <- function(fit = fit_coxpe_bart){
           R <- extract_R_PH(fit)
    num_time <- length(R)
           S <- list()
      S[[1]] <- R[[1]] * 0 + 1
  for (j in 2:(num_time + 1)) {
    S[[j]]   <- S[[j - 1]] * R[[j - 1]]
  }
  num_row    <- nrow(R[[1]])
  num_col    <- ncol(R[[1]])
  num_slice  <- length(R) + 1
  S          <- array(unlist(S), dim = c(num_row, num_col, num_slice))
  return(S)
}

