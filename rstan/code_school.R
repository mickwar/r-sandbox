library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

schools = list(J = 8,
    y = c(28, 8, -3, 7, -1, 1, 18, 12),
    sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

nburn = 20000
nmcmc = 30000
nchain = 8
fit = stan(file = "./8schools.stan", data = schools,
    iter = nburn + nmcmc, warmup = nburn, chains = nchain, control = list("adapt_delta" = 0.95))

# Number of divergent iterations, proportion of divergences
sapply(get_sampler_params(fit, inc_warmup=FALSE), function(x) sum(x[,'divergent__']))
sapply(get_sampler_params(fit, inc_warmup=FALSE), function(x) mean(x[,'divergent__']))

params = extract(fit, permuted = FALSE)
plot(params$theta[,1])
plot(params[,4,11])
tmp = c(params[,1,11],
    params[,2,11],
    params[,3,11],
    params[,4,11])

sum(duplicated(params[,3,1]))


plot(params$theta[,1], xlim = c(1, 100))
plot(tmp, xlim = c(1, 100))
plot(c(t(matrix(tmp, ncol = 4))), xlim = c(1, 100))

print(fit)
plot(fit, pars = c("mu", "tau", "eta", "theta"))
pairs(fit)

str(fit)

plot(fit, pars = c("mu", "tau", "eta", "theta"))
