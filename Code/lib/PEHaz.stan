data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> grid_size;
  vector[N] pop_haz;
  matrix[N,P] X;
  vector[N] y;
  vector[N] delta;
  vector[grid_size] grid;
}

parameters {
  vector[P] beta;
  vector<lower=0>[grid_size] lambda;
  real<lower=0> shape_lambda;
  real<lower=0> scale_lambda;
}

transformed parameters {
  real loglik;
  vector[N] eta;
  vector[grid_size - 1] base_cum_haz;
  {
    vector[N] base_haz;
    vector[N] cum_haz;

    eta = X * beta;
    for(i in 1:N) {
      base_haz[i] = 0;
      cum_haz[i] = 0;
      for(g in 1:(grid_size - 1)) {
        if(y[i] < grid[g]) {

        }
        if(y[i] >= grid[g] && y[i] < grid[g+1]) {
          cum_haz[i] += lambda[g] * (y[i] - grid[g]);
          base_haz[i] = lambda[g];
        }
        if(y[i] >= grid[g+1]) {
          cum_haz[i] += (grid[g+1] - grid[g]) * lambda[g];
        }
      }
      if(y[i] >= grid[grid_size]) {
        base_haz[i] = lambda[grid_size];
        cum_haz[i] += (y[i] - grid[grid_size]) * lambda[grid_size];
      }
    }
    loglik = 0;
    for(i in 1:N) {
      loglik += delta[i] * log(base_haz[i] * exp(eta[i]) + pop_haz[i]) -
        cum_haz[i] * exp(eta[i]);
    }
  }
  base_cum_haz[1] = (grid[2] - grid[1]) * lambda[1];
  for(g in 2:(grid_size - 1)) {
    base_cum_haz[g] = base_cum_haz[g-1] + (grid[g+1] - grid[g]) * lambda[g];
  }
}

model {
  target += loglik;
  lambda ~ gamma(shape_lambda, scale_lambda);
  beta ~ student_t(3, 0, 2);
}
