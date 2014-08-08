f = function(x)
    x^3-(x-4)^2+x + rnorm(length(x), 0, 2)
metric = function(x, y, type){
    if (type == "mse")
        return (mean((x-y)^2))
    if (type == "abs")
        return (mean(abs(x-y)))
    }

yy = f(xx)

xx = sort(runif(15, -3, 3))

mod = lm(yy ~ xx)
lin = function(x)
    coef(mod)[1] + coef(mod)[2] * x

type = "mse"
type = "abs"
met.line = metric(yy, lin(xx), type)
est1 = yy + met.line
# lin and est1 have same MSE
metric(yy, lin(xx), type)
metric(yy, est1, type)

plot(xx, yy)
lines(xx, lin(xx), col='red')
lines(xx, est1, col='green')

# get differences from data (f) and predictions
(dif.line = yy - lin(xx))
(dif.curv = yy - est1)
mean(dif.line)
mean(dif.curv)
var(dif.line)
var(dif.curv)
