lu_decomp = function(x){
    n = nrow(x)
    L = diag(n)
    U = matrix(0, n, n)
    
    U[1,] = x[1,]
    for (i in 2:n){
        L[i:n, i-1] = (x[i:n,i-1] - L[i:n,-(i-1)] %*% U[-(i-1),i-1]) / U[i-1,i-1]
        U[i, i:n] = x[i,i:n] - L[i,-i] %*% U[-i,i:n]
        }

    return (list("L"=L, "U"=U))
    }

x = matrix(c(2,1,5,3, -3, 4, 3, -6, 1, -3, -1, -3, 3, -3, -1, 1), 4, 4)

lu_decomp(x)

