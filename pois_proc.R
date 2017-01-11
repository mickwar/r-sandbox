library(mwBASE)
library(ismev)

ny = 1
n = 3650 * ny
x = rexp(n, 1)

u = 2.2

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

colMeans(out$params)

tmp = apply(out$params, 1, function(p) calc.post(y, p))
out3 = pp.fit(x, 4.6)

out$params[which.max(tmp),]
out3$mle

 -max(tmp)
out3$nllh

pp.diag(out3)


calc.post2 = function(y, param){
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

out2 = mcmc_sampler(y - u, calc.post2, 2, nburn = 100000, nmcmc = 100000)

mean(out$accept)
plot(out$params[,1], type='l')
plot(out$params[,2], type='l')
plot(out$params[,3], type='l')

mean(out2$accept)
plot(out2$params[,1], type='l')
plot(out2$params[,2], type='l')

plot(density(out$params[,3]))
lines(density(out2$params[,2]), col = 'red')

plot(density(out$params[,2] + out$params[,3]*(u - out$params[,1])))
lines(density(out2$params[,1]), col = 'red')

