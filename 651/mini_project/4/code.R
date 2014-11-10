mw.pairs = function(x){
    require(MASS)
    require(fields)
    par(mfrow=rep(ncol(x), 2), mar=rep(0, 4))
    for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
            if (i == j){
                #hist(x[,i], axes = FALSE, main = "", col = "gray")
                plot(density(x[,i]), axes = FALSE, main = "")
#               legend("topright", legend = i, box.lty = 0, cex = 1.5)
            } else {
                if (i > j){
                    z = kde2d(x[,j], x[,i])
                    plot(NA, xlim = range(z$x), ylim = range(z$y), axes = FALSE)
#                   .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20),
#                       col=gray(seq(0.0, 1.0, length=20)))
                    .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20),
                        col = tim.colors(20))
#                   .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20),
#                       col=two.colors(20, start="dodgerblue"))
                } else {
                    plot(x[,j], x[,i], type = 'l', axes = FALSE)
                    }
                }
            box()
            }
        }
    }

y = scan("~/files/R/651/data/faculty.dat")

k = length(y)
n = k
ybar = mean(y)

# The model:
# Y_i ~ N(theta_i, sigma^2)
# theta_i ~ N(mu, tau^2)
# sigma^2 ~ IG(a, b)
# mu ~ N(m, s^2)
# tau^2 ~ IG(c, d)

# priors
# for mu
m = 6       # prior mean
s2 = 0.5^2  # prior variance on the prior mean

# for sigma^2
a = 3
b = 2

# for tau^2
c = 3
d = 2

nburn = 1000
nmcmc = 10000
p.theta = matrix(0, nburn + nmcmc, k)
p.mu = double(nburn + nmcmc)
p.sig2 = double(nburn + nmcmc)
p.tau2 = double(nburn + nmcmc)

# init
p.theta[1,] = m
p.mu[1] = m
p.sig2[1] = 1/(b * (a+1))   # mode for IG
p.tau2[1] = 1/(d * (c+1))

for (i in 2:(nburn+nmcmc)){
    # update thetas
    p.theta[i,] = rnorm(n, (y*p.tau2[i-1] + p.mu[i-1]*p.sig2[i-1]) /
        (p.tau2[i-1] + p.sig2[i-1]), sqrt(p.tau2[i-1]*p.sig2[i-1] /
        (p.tau2[i-1] + p.sig2[i-1])))

    # update mu
    p.mu[i] = rnorm(1, (mean(p.theta[i,])*k*s2 + m*p.tau2[i-1]) /
        (k*s2+p.tau2[i-1]), sqrt(s2*p.tau2[i-1] / (k*s2+p.tau2[i-1])))

    # update sigma^2
    p.sig2[i] = 1/rgamma(1, shape = a + n/2, scale = 1/(1/b + 0.5 *
        sum((y - p.theta[i,])^2)))

    # update tau^2
    p.tau2[i] = 1/rgamma(1, shape = c + k/2, scale = 1/(1/d + 0.5 *
        sum((p.theta[i,] - p.mu[i])^2)))
    }


for (i in 1:k)
    plot(p.theta[,i], type='l')

plot(p.mu, type='l')
plot(p.sig2, type='l')
plot(p.tau2, type='l')

mw.pairs(cbind(p.mu, p.sig2, p.tau2))
mw.pairs(p.theta[,1:5])

rr = cor(cbind(p.theta, p.mu, p.sig2, p.tau2))

library(fields)
image.plot(abs(rr))

preds = double(nburn + nmcmc)
for (i in 1:(nburn+nmcmc)){
    temp.theta = rnorm(1, p.mu[i], sqrt(p.tau2[i]))
    preds[i] = rnorm(1, temp.theta, sqrt(p.sig2[i]))
    }

hist(y, col='gray', freq=FALSE, breaks=6, ylim=c(0, 0.8), xlim=c(3.5,7.5))
points(density(y), col='red', type='l', lwd=3)
points(density(preds), col='black', type='l', lwd=3)
plot(density(preds), col='red', type='l')


