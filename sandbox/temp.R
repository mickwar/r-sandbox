cross = function(...){
    vec = list(...)
    d = length(vec)
    N = double(d+2) + 1
    for (i in 1:d)
        N[i+1] = length(vec[[i]])
    out = matrix(0, prod(N), d)
    for (i in 1:d){
        out[,i] = rep(vec[[i]], times=prod(N[1:i]),
            each=prod(N[(i+2):(d+2)]))
        }
    return(out)
    }

xx = seq(0, 1, length=51)
yy = 4*seq(-1, 1, length=101)^2

A = cross(xx, yy)
f = function(x)
    x[1] ^ x[2]

B = cbind(A, apply(A, 1, f))

library(rgl)

plot3d(B)
