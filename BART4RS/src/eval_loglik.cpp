#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix eval_loglik(NumericVector& y,
                          NumericVector& delta,
                          NumericMatrix& lambda,
                          NumericVector& grid,
                          NumericMatrix& base_haz,
                          NumericMatrix& cum_base_haz,
                          NumericVector& pop_haz,
                          std::vector<int>& obs_to_bin) {

  int num_obs = y.size();
  int grid_size = grid.size();
  int num_iter = lambda.nrow();

  NumericMatrix loglik(num_iter, num_obs);
  for(int t = 0; t < num_iter; t++) {
    for(int i = 0; i < num_obs; i++) {
      int bin = obs_to_bin[i];
      double cum_base_haz_Y = (y(i) - grid(bin)) * base_haz(t, bin);
      if(bin > 0) {
        cum_base_haz_Y += cum_base_haz(t, bin - 1);
      }
      loglik(t,i) = delta(i) *
        log(base_haz(t, bin) * lambda(t,i) + pop_haz(i))
        - lambda(t,i) * cum_base_haz_Y;
    }
  }

  return loglik;

}
