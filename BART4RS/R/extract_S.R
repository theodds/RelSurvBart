extract_S <- function(fit) {
  R <- extract_R(fit)
  num_time <- length(R)
  S <- list()
  S[[1]] <- R[[1]] * 0 + 1
  for(j in 2:(num_time + 1)) {
    S[[j]] <- S[[j-1]] * R[[j-1]]
  }
  num_row <- nrow(R[[1]])
  num_col <- ncol(R[[1]])
  num_slice <- length(R) + 1
  S <- array(unlist(S), dim = c(num_row, num_col, num_slice))
  return(S)
}
