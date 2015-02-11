##### MCMC on the beetles data
# I tried considering the log base for the link function (having logit form)
# as a model parameter. Bad things happened.


dose = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
n.beetles = c(59, 60, 62, 56, 63, 59, 62, 60)
y.beetles = c(6, 13, 18, 28, 52, 53, 61, 60)

g.inv = function(x, beta, base = exp(1))
    base^(x %*% beta) / (1 + base^(x %*% beta))

lik = function(p, n, y)
    prod(p^y * (1-p)^(n-y))

loglik = function(p, n, y)
    sum(log(p)*y + log(1-p)*(n-y))

prior = function(beta, delta)
    prod(pnorm(beta, 0, 5))#*dgamma(delta, 1, 1)

logprior = function(beta, delta)
    sum(pnorm(beta, 0, 5, log.p = TRUE))# + dgamma(delta, 1, 1, log = TRUE)

x = cbind(1, dose)
beta = c(1, 1)
delta = exp(1) # logarithm base

# center the x
#x = x - matrix(rep(apply(x, 2, mean), each = length(dose)), ncol = 2)
x[,2] = x[,2] - mean(x[,2])

full.loglik = sum((log(y.beetles/n.beetles)*y.beetles + 
    log(1 - y.beetles/n.beetles)*(n.beetles - y.beetles))[1:7])
c.sq = c(0.6, 10, 0.01)

ITERS = 100000
beta.save = matrix(0, 2, ITERS)
delta.save = double(ITERS)
accept = rep(0, 3)

beta.save[,1] = beta
delta.save[1] = delta
delta.cur = delta.save[1]

for (i in 2:ITERS){
    for (j in 1:2){
        delta.cur = delta
        beta.cur = beta
        beta.samp = rnorm(1, beta.cur[j], c.sq[j])
        beta.star = beta.cur
        beta.star[j] = beta.samp

        logalpha = logprior(beta.star, delta.cur) +
            loglik(g.inv(x, beta.star, delta.cur), n.beetles, y.beetles) -
            logprior(beta.cur, delta.cur) -
            loglik(g.inv(x, beta.cur, delta.cur), n.beetles, y.beetles)

#       logprior(beta.star, delta.cur)
#       loglik(g.inv(x, beta.star, delta.cur), n.beetles, y.beetles)
#       logprior(beta.cur, delta.cur)
#       loglik(g.inv(x, beta.cur, delta.cur), n.beetles, y.beetles)

        if (log(runif(1)) < logalpha){
            beta = beta.star
            accept[j] = accept[j] + 1
        } else {
            beta = beta.cur
            }
        }

    delta.cur = delta
    delta.star = rnorm(1, delta.cur, c.sq[3])
    logalpha = logprior(beta.cur, delta.star) +
        loglik(g.inv(x, beta.cur, delta.star), n.beetles, y.beetles) -
        logprior(beta.cur, delta.cur) -
        loglik(g.inv(x, beta.cur, delta.cur), n.beetles, y.beetles)
    
    if (log(runif(1)) < logalpha){
        delta = delta.star
        accept[3] = accept[3] + 1
    } else {
        delta = delta.cur
        }

    beta.save[,i] = beta
    delta.save[i] = delta
    }

par(mfrow=c(1,3))
plot(beta.save[1,], type='l')
plot(beta.save[2,], type='l')
plot(delta.save, type='l')
accept / ITERS
