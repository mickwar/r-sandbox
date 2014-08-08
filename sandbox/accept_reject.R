# Example 5.6.9 (p. 254) of Casella and Berger
a1 = 2
b1 = 6
a2 = 2.7
b2 = 6.3


xx = seq(0, 1, length=10000)

# calculate M ( sup_y f_Y(y) / f_V(y) )
# M must be finite for Theorem 5.6.8 to hold
yy1 = dbeta(xx, a1, b1)
yy2 = dbeta(xx, a2, b2)
mm = yy2/yy1
M = max(mm, na.rm=TRUE)
plot(xx, mm, type='l')

curve(dbeta(x, a1, b1), ylim=c(0, M*max(c(yy1, yy2))))
curve(dbeta(x, a2, b2), add=TRUE, col='red')
curve(M*dbeta(x, a1, b1), add=TRUE, col='blue')

n = 100000
U = runif(n, 0, 1)
V = rbeta(n, a1, b1)

# target
Y = V[U < 1/M * dbeta(V, a2, b2) / dbeta(V, a1, b1)]
points(density(Y), col='green', type='l', lwd=3)
