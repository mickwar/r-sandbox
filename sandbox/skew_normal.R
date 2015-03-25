### random draws from a skew normal distribution

n = 100000
delta = 0.0
u_0 = rnorm(n, 0, 1)
v = rnorm(n, 0, 1)

u_1 = delta*u_0 + sqrt(1-delta^2)*v

z = ifelse(u_0 >= 0, u_1, -u_1)

plot(density(z))

alpha = delta / sqrt(1-delta^2)

dskewnormal = function(x, ksi, omega, alpha)
    2/omega * dnorm(x, ksi, omega) * pnorm(alpha*(x-ksi)/omega)

curve(dskewnormal(x, 0, 1, alpha), add=TRUE, from = -5, to = 5, col='red')
