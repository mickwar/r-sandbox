dose = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
n.beetles = c(59, 60, 62, 56, 63, 59, 62, 60)
y.beetles = c(6, 13, 18, 28, 52, 53, 61, 60)

g.inv = function(x, beta, base = exp(1))
    base^(x %*% beta) / (1 + base^(x %*% beta))

lik = function(p, n, y)
    prod(p^y * (1-p)^(n-y))

loglik = function(p, n, y)
    sum(log(p)*y + log(1-p)*(n-y))
