rgp = function(n, mu, sigma, ksi)
    mu + sigma/ksi*(runif(n)^(-ksi)-1)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

theta = 0.5
n = 10000
ww = -1/log(runif(n))
y = double(n)
y[1] = ww[1] / theta
for (i in 2:n)
    y[i] = max((1-theta)*y[i-1], ww[i])
plot(y)


y = rnorm(n)
u = quantile(y, 0.95)



ex = y[y > u] - u
N = length(ex)



set.seed(1)
N = 500
ksi = 0.0
sigma = 1
ex = rgp(N, 0, sigma, ksi)
ex = sort(ex, decreasing = TRUE)
dat = list("N" = N, "exceedances" = ex)


nburn = 1000
nmcmc = 1000
nchain = 1
fit = stan(file = "./extreme.stan", data = dat,
    iter = nburn + nmcmc, warmup = nburn, chains = nchain,
    init = list(list("ksi" = 0.0, "sigma" = 1)),
    control = list("adapt_delta" = 0.80))

# Number of divergent iterations, proportion of divergences
sapply(get_sampler_params(fit, inc_warmup=FALSE), function(x) sum(x[,'divergent__']))
sapply(get_sampler_params(fit, inc_warmup=FALSE), function(x) mean(x[,'divergent__']))

ind = c(sapply(get_sampler_params(fit, inc_warmup=FALSE), function(x) which(x[,'divergent__'] == 1)))
ind = 1

params = extract(fit, permuted = TRUE)

plot(params$ksi[-ind], params$sigma[-ind], col = 'blue',
    xlim = range(params$ksi), ylim = range(params$sigma))
points(params$ksi[ind], params$sigma[ind], col = 'blue')
points(mean(params$ksi), mean(params$sigma), col = 'red', pch = 15)
points(ksi, sigma, pch = 16, col = 'green')

plot(cumsum(params$ksi)/(1:nmcmc))
abline(h = -0.5)

plot(cumsum(params$sigma)/(1:nmcmc))
abline(h = 1)

mean(params$ksi)
mean(params$sigma)
plot(params


plot(params$theta[,1], xlim = c(1, 100))
plot(tmp, xlim = c(1, 100))
plot(c(t(matrix(tmp, ncol = 4))), xlim = c(1, 100))

print(fit)
plot(fit, pars = c("mu", "tau", "eta", "theta"))
pairs(fit)

str(fit)

plot(fit, pars = c("mu", "tau", "eta", "theta"))
