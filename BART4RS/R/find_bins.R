find_bins <- function(y, time_grid) {
  for(i in seq_len(length(time_grid) - 1)) {
    if((y > time_grid[i]) & (y <= time_grid[i+1])) return(i)
  }
  return(length(time_grid)-1)
}
