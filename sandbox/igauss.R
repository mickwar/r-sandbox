# Normal-inverse Gaussian distribution

digauss = function(x, mu, alpha, beta, delta, log = TRUE){
    gamma = sqrt(alpha^2 - beta^2)
    out = log(alpha) + log(delta) + log(besselK(alpha*sqrt(delta^2+(x-mu)^2) , 1)) -
        log(pi) - 0.5*log(delta^2+(x-mu)^2) + delta*gamma + beta*(x-mu)
    if (log)
        return (out)
    return (exp(out))
    }

# standard normal
x = seq(-10, 10, length=100)
plot(x, digauss(x, 0, 4, 4, 2, FALSE), type='l')
curve(dnorm(x, 0, 1), add = TRUE, col='red')


