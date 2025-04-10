// functions{
//   real R_hat(real alpha, real beta, real Mass) {
//     R_hat = alpha * Mass ^ beta;
//     return R_hat;
//   }
// }

data {
  int N_obs;
  vector[N_obs] Mass;
  vector[N_obs] Resp;
}

parameters {
  real<lower = 0, upper = 1> beta;
  real<lower = 0> alpha;
  real<lower = 0> lambda;
}

transformed parameters {
  vector[N_obs] log_R_hat;
  vector[N_obs] R_hat;
  real log_alpha;
  log_alpha = log(alpha);
  log_R_hat = log_alpha + beta * log(Mass);
  R_hat = exp(log_R_hat);
  // R_hat = alpha * Mass ^ beta;
}

model {
  alpha ~ normal(0,1);
  beta ~ normal(0.75,0.2);
  lambda ~ normal(100,100);
  target += gamma_lpdf(Resp | R_hat * lambda, lambda);
}

