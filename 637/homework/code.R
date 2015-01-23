n = 50;

### generate from a normal with a covariate
set.seed(1)    # seed

x = cbind(1, runif(n))
beta = c(1.5, 0.75)
xbeta = x %*% beta;
sigma = 0.10; # standard deviation

y = rnorm(n, xbeta, sigma)
z = rnorm(n, sum(c(1, 0.5) * beta), sigma)

mod11 = glm(y ~ 1)
mod12 = glm(y ~ 0 + x)

mod21 = glm(z ~ 1)
mod22 = glm(z ~ 0 + x)

s2 = sum((z - x %*%coef(mod22))^2)/n

# don't divide by s2
deviance(mod11) / s2
deviance(mod12) / s2

deviance(mod21) / s2
deviance(mod22) / s2
