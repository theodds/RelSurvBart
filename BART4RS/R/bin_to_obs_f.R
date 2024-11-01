bin_to_obs_f <- function(y_bin, num_bins) {
  out <- vector(num_bins, mode = 'list')
  for(k in seq_len(num_bins)) {
    out[[k]] <- which(y_bin == k)
  }
  return(out)
}
