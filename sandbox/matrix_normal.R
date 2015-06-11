dmatnorm = function(X, M, U, V, log = FALSE){
    tr = function(x) sum(diag(x))
    n = NROW(X)
    p = NCOL(X)
    numer = -0.5 * tr(solve(V) %*% t(X - M) %*% solve(U) %*% (X - M))
    denom = (n*p/2)*log(2*pi) + (n/2)*determinant(V)$modulus[1] +
        (p/2)*determinant(U)$modulus[1]
    if (log){
        return (numer - denom)
    } else {
        return (exp(numer - denom))
        }
    }

rmatnorm = function(M, U, V){
    n = NROW(M)
    p = NCOL(M)
    A = t(chol(U))
    B = chol(V)
    X = matrix(rnorm(n*p, 0, 1), n, p)
    return (M + A %*% X %*% B)
    }

n = 100
p = 2
M = matrix(0, n, p)
M[1,] = c(0, 2)
#U = matrix(0.5, n, n) + 0.5*diag(n)
U = diag(n)
V = matrix(-0.8, p, p)
diag(V) = 1.0

plot(rmatnorm(M, U, V), pch=20, xlim=c(-6,6), ylim=c(-6,6),
    col = 1:1000)
for (i in 1:10000)
    points(rmatnorm(M, U, V), pch=20, col= 1:1000)

