data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] pop_haz;
  matrix[N,P] X;
  vector[N] y;
  vector[N] delta;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[P] beta;
  real<lower=0> weibull_shape;
  real<lower=0> weibull_scale;
}

transformed parameters {
  real loglik;
  {
    vector[N] eta;
    vector[N] base_haz;
    vector[N] cum_haz;

    eta = X  * beta;
    for(i in 1:N) {
      base_haz[i] = pow(y[i] / weibull_scale, weibull_shape- 1) *
        weibull_shape / weibull_scale;
      cum_haz[i] = pow(y[i] / weibull_scale, weibull_shape);
    }

    loglik = 0;
    for(i in 1:N) {
      loglik += delta[i] * log(base_haz[i] * exp(eta[i]) + pop_haz[i]) -
        cum_haz[i] * exp(eta[i]);
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += loglik;
  beta ~ student_t(3, 0, 2);
  weibull_shape ~ uniform(0.1, 10);
  weibull_scale ~ gamma(1, 1);
}

