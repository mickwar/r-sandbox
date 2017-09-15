### The integrals

library(rgl)

z = seq(1e-6, 1 - 1e-6, length = 101)
theta = pi * z / 2

alpha = 2
beta = 3

plot(z, dbeta(z, alpha, beta))
plot(theta, dbeta(2*theta/pi, alpha, beta) * 2/pi)

f = function(alpha, beta){
    g = function(x) ifelse(x < pi/4, tan(x)*dbeta(2*x/pi, alpha, beta)*2/pi,
        1/tan(x)*dbeta(2*x/pi, alpha, beta)*2/pi)
    integrate(g, 1e-12, pi/2-1e-12, subdivisions = 500L)$value
    }



ab = expand.grid(seq(-7, 5, length = 50), seq(-7, 5, length = 50))
z = apply(ab, 1, function(x) f(exp(x[1]), exp(x[2])))

ab = expand.grid(seq(2, 8, length = 50), seq(2, 8, length = 50))
z = apply(ab, 1, function(x) f(exp(x[1]), exp(x[2])))

plot3d(exp(ab[,1]), exp(ab[,2]), z, zlim = c(0,1), xlim = c(0, 10), ylim = c(0, 10))

plot3d(exp(ab[,1]), exp(ab[,2]), z, zlim = c(0,1))

plot3d(exp(ab[,1]), exp(ab[,2]), log(z), xlim = c(0, 10), ylim = c(0, 10))

plot3d(ab[,1], ab[,2], log(z))

