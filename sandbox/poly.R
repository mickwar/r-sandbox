set.seed(2)

gen.y = function(x)
    sin(x) + rnorm(length(x), 0, 0.25)

n = 12
x = sort(runif(n, 0, 3*pi))
x = sort(seq(0, 3*pi, length=n))

y = gen.y(x)

# data
plot(x, y, pch=20, xlim = c(0, 3*pi), ylim=c(-1.5, 1.5))

# truth
xx = seq(0, 3*pi, length=1000)
lines(xx, sin(xx), pch=20, lty=2, lwd=3)

mod1 = lm(y ~ x)
pred1 = cbind(1, xx) %*% mod1$coef
lines(xx, pred1, lwd=2)

newD = x
newX = cbind(xx)
for (i in 2:(length(x)-1)){
    readline(paste(i,"/",length(x)))
    newD = cbind(newD, x^i)
    newX = cbind(newX, xx^i)
    mod = lm(y ~ 0 + newD)
    pred = newX %*% mod$coef
    lines(xx, pred, lwd=2, col=i)
    }
points(x, y, pch=20, lwd=5)
