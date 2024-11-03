modbart_density_approx <- function(y, delta, tgrid, sgrid) {
  sfun <- approxfun(x = tgrid, y = sgrid)
  epsilon <- 0.001
  if(delta == 0) {
    return(sfun(y))
  }
  approx_deriv <- (sfun(y + epsilon) - sfun(y)) / epsilon
  return(-approx_deriv)
}