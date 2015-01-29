set.seed(1)
n = 100
x = cbind(runif(n, -2, 4), runif(n, 0, 4))
x = cbind(x, x[,1] * x[,2], x[,1]^2)
y = double(n)
for (i in 1:n)
    y[i] = 1.5 - 4*x[i,1] + 2*x[i,2] - 0.5*x[i,3] + 1.0*x[i,4] + rnorm(1, 0, 1)

pairs(cbind(x, y))
mod1 = lm(y ~ x)

ystar = (y - mean(y))/1
mod2 = lm(ystar ~ x)

summary(mod1)
summary(mod2)

coef(mod1)
coef(mod2)
