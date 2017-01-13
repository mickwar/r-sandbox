library(mwBASE)
library(mwEVT)
library(ismev)

set.seed(2)
true.theta = 0.5
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

de.x = decluster(x, u = u, r = 1, xlim = c(0, 1000), log = "y")

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


calc.post.declust = function(y, param){
    sig = param[1]
    ksi = param[2]

    na = length(y)

    if (sig <= 0)
        return (-Inf)

    if (any(1 + ksi * y/sig < 0)) 
        return(-Inf)

    # likelihood
    if (ksi != 0) {
        out = -na * log(sig) - (1 + 1/ksi) * sum(log(1 + 
            ksi * y/sig))
    } else {
        out = -na * log(sig) - sum(y)/sig
        }

    # piors
#   out = out - log(sig)
#   out = out + dnorm(ksi, 0, 1, log = TRUE)

    return(out)
    }

out = mcmc_sampler(y, calc.post, 3, nburn = 100000, nmcmc = 100000)
out2 = mcmc_sampler(de.x$max_cluster_excess, calc.post.declust, 2, nburn = 100000, nmcmc = 100000)
out3 = mcmc_sampler(de.x$max_cluster_excess + u, calc.post, 3, nburn = 100000, nmcmc = 100000)

mean(out$accept)
plot(out$params[,1], type='l')
plot(out$params[,2], type='l')
plot(out$params[,3], type='l')

mean(out2$accept)
plot(out2$params[,1], type='l')
plot(out2$params[,2], type='l')

mu = out$params[,1]
sig = out$params[,2]
ksi = out$params[,3]

mu.star = mu - sig/ksi*(1-theta^(-ksi))
sig.star = sig*theta^ksi

mu.star = mu + sig*theta^(-ksi)/ksi*(1-theta^(-ksi))

par(mfrow = c(2,2))
plot(density(mu))
plot(density(sig))
plot(density(mu.star))
plot(density(sig.star))

plot(density(ksi))
plot(density(theta))

plot(density(mu.star))
lines(density(out3$params[,1]), col = 'blue')

plot(density(out$params[,3]))
lines(density(out2$params[,2]), col = 'red')
lines(density(out3$params[,3]), col = 'blue')

plot(density(sig.star))
lines(density(out$params[,2]), col = 'black')
lines(density(out2$params[,1]), col = 'red')
lines(density(out3$params[,2]), col = 'blue')
lines(density(out3$params[,2]*theta^out3$params[,3]), col = 'blue')
