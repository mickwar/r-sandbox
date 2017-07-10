data {
  int<lower=0> N;
  real exceedances[N];
  real threshold;
}
parameters {
  real ksi; 
  real<lower=0> sigma;
}
transformed parameters {
  real z[N];
  for (i in 1:N)
    z[i] = -exceedances[i]*ksi/sigma;
}
model {
  ksi ~ normal(0, 1);
  sigma ~ inv_gamma(1, 1);
  if (ksi > 0)
    target += pareto_type_2_lpdf(exceedances | 0, sigma/ksi, 1/ksi);
  else if (ksi == 0)
    target += exponential_lpdf(exceedances | 1/sigma);
  else
    target += beta_lpdf(z | 1, -1/ksi);
}
