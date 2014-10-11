get.a.star = function(n, test.mu, size = 1000, p = 6, mu, sig2){
    require(MASS)
    out = matrix(0, n, p)
    if (missing(mu))
        mu = rep(0, p)
    if (missing(sig2))
        sig2 = diag(p)
    for (i in 1:n){
        x = mvrnorm(size, mu, sig2)

        # discriminant function
        xbar = apply(x, 2, mean)
        S = var(x)
        a = solve(S) %*% (xbar - test.mu)

        out[i,] = diag(sqrt(diag(S))) %*% a
        }
    return (out)
    }

n = 10000
size = 500
p = 7
test.mu = c(1, 0.5, -0.5, 0, 0, 0, 0)
mu = rep(0, p)
sig2 = diag(p)
x = get.a.star(n, test.mu, size, p, mu, sig2)


plot(density(x[,1]), ylim=c(0, 10), xlim=c(-1.3, .3))
for (i in 1:ncol(x))
    points(density(x[,i]), type='l', col=i)

(ms = apply(x, 2, mean))
(ss = apply(x, 2, sd))

for (i in 1:ncol(x))
    curve(dnorm(x, ms[i], ss[i]), add = TRUE, lwd=3, col=i)

