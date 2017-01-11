library(mwBASE)
library(ismev)

set.seed(2)
true.theta = 1.0
ny = 50
n = 365 * ny
ww = -1/log(runif(n))
x = double(n)
x[1] = ww[1] / true.theta
for (i in 2:n)
    x[i] = max((1-true.theta)*x[i-1], ww[i])

plot(log(x), xlim = c(1, 100))

u = 250
Tu = diff(which(x > u))
N = sum(x > u)
m1 = sum(Tu == 1)

m1 / N

B = 100000
rs = rbeta(B, N - m1 + 1, m1 + 2)
qs = rbeta(B, sum(Tu) -2*N+2+m1+1,N-m1)

theta = rs / qs
plot(density(theta))


y = x[x > u]
na = length(y)

na / n

calc.post = function(y, param){
    mu = param[1]
    sig = param[2]
    ksi = param[3]

    if (sig <= 0)
        return (-Inf)

    if (any(1 + ksi*(y - mu)/sig <= 0))
        return (-Inf)

    na = length(y)

    # likelihood
#   if (ksi != 0){
        out = -ny * (1 + ksi*(u - mu)/sig)^(-1/ksi) - na * log(sig) +
            (-1/ksi-1)*sum(log(1 + ksi*(y - mu)/sig))
#   } else {

#       }

    if (is.na(out))
        return (-Inf)

    # priors
#   out = out + dnorm(mu, 0, 10, log = TRUE)
#   out = out - log(sig)
#   out = out + dnorm(ksi, 0, 1, log = TRUE)

    return (out)
    }

out = mcmc_sampler(y, calc.post, 3, nburn = 100000, nmcmc = 100000)

mean(out$accept)
plot(out$params[,1], type='l')
plot(out$params[,2], type='l')
plot(out$params[,3], type='l')

mu = out$params[,1]
sig = out$params[,2]
ksi = out$params[,3]

mu.star = mu - sig/ksi*(1-theta^(-ksi))
sig.star = sig*theta^ksi

par(mfrow = c(2,2))
plot(density(mu))
plot(density(sig))
plot(density(mu.star))
plot(density(sig.star))

plot(density(ksi))
plot(density(theta))
