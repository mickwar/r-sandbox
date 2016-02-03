digamma = function(x, alpha, beta)
    alpha*log(beta) - lgamma(alpha) -(alpha + 1)*log(x) - beta / x

set.seed(1)
n = 500
rho = 0.95
sig = 2
eps = rnorm(n, 0, sig)

y = double(n)
y[1] = eps[1]
for (i in 2:n)
    y[i] = rho*y[i-1] + eps[i]

plot(y, type='o', pch = 20)

sig^2 * 1/(1-0.95^2)
var(y)

sum(dnorm(y[2:n], 0.95*y[1:(n-1)], sig, log = TRUE))

# z[1] "=" y_{-1}
z = c(0, y)


alpha = n/2
beta = 1/2*sum((z[2:(n+1)] - rho*z[1:n])^2)

ndraws = 50000
post.sig2 = 1/rgamma(ndraws, alpha, beta)
post.rho = rnorm(ndraws, sum(z[2:(n+1)] * z[1:n]) / sum(z[1:n]^2),
    sqrt(post.sig2 / (sum(z[1:n]^2))))

par(mfrow=c(2,1), mar = c(4.1, 2.1, 2.1, 1.1))
xx = seq(0.001, 10, length = 500)
plot(xx, exp(digamma(xx, alpha, beta)), type='l', xlab = expression(sigma^2),
    ylab = "")
lines(density(post.sig2), col = 'red')

xx = seq(0, 1, length = 500)
plot(xx, dnorm(xx, sum(z[2:(n+1)] * z[1:n]) / sum(z[1:n]^2),
    sqrt(alpha/beta / (sum(z[1:n]^2)))), type='l', xlab = expression(rho))
lines(density(post.rho), col = 'red')

beta/(alpha - 1)
