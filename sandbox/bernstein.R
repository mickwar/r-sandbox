make.b = function(x, n){
    out = matrix(0, length(x), n+1)
    for (v in 0:n)
        out[,v+1] = choose(n, v) * x^v * (1-x)^(n-v)
    return (out)
    }


xx = seq(0, 1, length = 100)
X = make.b(xx, 4)
Y = make.b(xx, 10)

ex = svd(X)
ey = svd(Y)

par(mfrow = c(2,1), mar = c(2.1,2.1,2.1,2.1))
matplot(xx, X, type='l')
matplot(xx, Y, type='l')

matplot(ex$v, type='l')
matplot(ey$v, type='l')

matplot(ex$u)
matplot(ey$u)

head(ex$u)
head(ey$u)

(ex$u %*% diag(ex$d))[,1:3] - (ey$u %*% diag(ey$d))[,1:3]

hist((ex$u %*% diag(ex$d))[,3] - (ey$u %*% diag(ey$d))[,3])
hist(ex$u[,3] - ey$u[,3])

beta = sort(rnorm(NCOL(X), 5, 3))
y = X %*% beta
plot(xx, y, pch = 20)
