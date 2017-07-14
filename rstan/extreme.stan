data {
  int<lower=0> N;
  vector[N] exceedances;
}
parameters {
  real ksi_raw; 
  real<lower=0> sigma_raw;
}
transformed parameters {
    vector[2] sig_adj;
    real ksi;
    real sigma;
    vector[N] ex_scaled;

    ksi = 0.01;
    sig_adj[1] = -ksi*max(exceedances);
    sig_adj[2] = 0.0;
    sigma = sigma_raw + max(sig_adj);
    ksi = ksi_raw - sigma / max(exceedances);
    ex_scaled = -ksi/sigma * exceedances;
}
model {
  ksi ~ normal(0, 1);
  sigma ~ inv_gamma(1, 1);
  if (ksi > 0)
    target += pareto_type_2_lpdf(exceedances | 0, sigma/ksi, 1/ksi);
  else if (ksi == 0)
    target += exponential_lpdf(exceedances | 1/sigma);
  else if (ksi < 0)
    target += beta_lpdf(ex_scaled | 1, -1/ksi);
}
