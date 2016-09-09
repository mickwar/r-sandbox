
# S = [0, 1]
fac = 0
rv = function(n, s)
    rbeta(n, 1+fac*(1-s), 1+fac*s)
rpareto = function(n)
    1/runif(n)

for (i in seq(0, 1, length = 101)){
    curve(dbeta(x, 1 + 3*(1-i), 1 + 3*i), ylim = c(0, 4))
    readline()
    }

n = 10000
y = rpareto(n)
v = rv(n, 0.5)

w0 = dbeta(1, 1+fac, 1)
w = y*v

w

n = 1000000
fac = 1000000
plot(density(rv(n, runif(n))))
lines(density(runif(n)), col = 'red')
(rv(n, runif(n)))

plot(density(y))

plot(density(y*v

mean(w/max(w))
mean(max(w)/w0 > 1)




bn = functio(n) qnorm(1 - 1/n)
an = function(n, ksi)
    ifelse(ksi != 0, ksi*(bn(2*n) - bn(n)) / (2^ksi - 1), (bn(2*n) - bn(n)) / log(2))

n = 100
m = 10000
x = apply(matrix(rnorm(n*m), n, m), 2, max)
y = (x - bn(n)) / an(n, 0)

(beta = sqrt(var(y)*6/(pi^2)))
(mu = mean(y) - beta * (-digamma(1)))
z = mu - beta*log(-log(runif(m)))

plot(density(y))
lines(density(z), col = 'red')

