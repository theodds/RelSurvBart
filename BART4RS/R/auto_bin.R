auto_bin <- function(y, factor = 0.5, time_grid = NULL) {
  N <- length(y)

  if(is.null(time_grid)) {
    obs_per_bin <- factor * N^(2/3)
    num_bins <- floor(N / obs_per_bin)
    time_grid <- quantile(y, probs = seq(1 / num_bins, 1, length = num_bins))
    time_grid <- c(0, time_grid)
  } else {
    stopifnot(time_grid[1] == 0)
    num_bins <- length(time_grid) - 1
  }

  obs_to_bin <- sapply(y, find_bins, time_grid = time_grid)
  bin_to_obs <- bin_to_obs_f(obs_to_bin, num_bins)

  bin_width <- time_grid[-1] - time_grid[-(num_bins + 1)]

  return(list(
    time_grid = time_grid,
    obs_to_bin = obs_to_bin - 1,
    bin_to_obs = sapply(bin_to_obs, \(i) i - 1),
    num_bins = num_bins,
    bin_width = bin_width
  ))

}
