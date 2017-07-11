data {
  int<lower=0> N;
  real exceedances[N];
}
parameters {
  real ksi; 
  real<lower=0> sigma;
}
model {
    ksi ~ normal(0, 1);
    sigma ~ inv_gamma(1, 1);
    if ((ksi < 0) && (1+ksi*max(exceedances)/sigma < 0))
        target += -1e10;
    else {
        target += -N*log(sigma);
        for (i in 1:N)
            target += -(1/ksi+1)*log(1+ksi*exceedances[i]/sigma);
    }
}
