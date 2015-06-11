set.seed(1)
n = 10
x = seq(0, 1, length = n)
y = 2-(x - 0.7)^2

# simple regression model
pred.y = predict(lm(y ~ x))

# offset model
const.y = y + 0.0801454

# mistake at one location
miss.y = y
at = 1
miss.y[at] = miss.y[at] + mean(abs(pred.y - y) / y)*miss.y[at]*n

plot(x, y, pch=20, ylim=range(c(y, pred.y, const.y, miss.y)), lwd=3)
lines(x, pred.y, lty=2, col='red')
lines(x, const.y, lty=2, col='blue')
lines(x, miss.y, lty=2, col='green')

# goodness of fit measure
mean(abs(pred.y - y) / y)
mean(abs(const.y - y) / y)
mean(abs(miss.y - y) / y)
