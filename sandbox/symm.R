### Decompose a square matrix into its symmetric and skew-symmetric parts
sks_decomp = function(x){
    n = NROW(x)
    symm = matrix(0, n, n)
    skew = matrix(0, n, n)
    A = solve(matrix(c(1,1,1,-1),2,2))

    diag(symm) = diag(x)
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            y = A %*% c(x[i,j], x[j,i])
            symm[i,j] = y[1]
            skew[i,j] = y[2]
            }
        }

    # Copy upper triangles to lower triangles
    symm[lower.tri(symm)] = t(symm)[lower.tri(symm)]
    skew[lower.tri(skew)] = -t(skew)[lower.tri(skew)]

    return (list("symmetric"=symm, "skew"=skew))
    }

### Alternate vresion (slower)
# sks_decomp2 = function(x){
#     n = NROW(x)
#     symm = matrix(0, n, n)
#     skew = matrix(0, n, n)
#     A = solve(matrix(c(1,1,1,-1),2,2))
# 
#     diag(symm) = diag(x)
#     for (i in 1:(n-1)){
#         for (j in (i+1):n){
#             y = A %*% c(x[i,j], x[j,i])
#             symm[i,j] = symm[j,i] = y[1]
#             skew[i,j] = y[2]
#             skew[j,i] = -y[2]
#             }
#         }
# 
#     return (list("symmetric"=symm, "skew"=skew))
#     }
