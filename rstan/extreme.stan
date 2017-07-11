data {
  int<lower=0> N;
  real exceedances[N];
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
    exceedances ~ pareto_type_2(0, sigma/ksi, 1/ksi);
  else if (ksi == 0)
    exceedances ~ exponential(1/sigma);
  else
    z ~ beta(1, -1/ksi);
}
