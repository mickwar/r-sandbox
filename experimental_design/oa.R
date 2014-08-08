# see Tang (1993)


oa = function(symbols, columns, strength=2, lambda=1){
    s = symbols  # v
    if (s < 2)
        stop("s must be >= 2.")
    m = columns  # k
    r = strength # t
    if (r > m)
        stop("r must be <= m.")
    l = lambda
    n = l * s^r 
    s = 1:s
    combo = matrix(0, choose(m, r), r)
    index = 1:r
    for (i in 1:nrow(combo)){
        combo[i, ] = index
        index[r] = index[r] + 1
        if (r > 1){
            for (j in r:2){
                if (index[j] > (m+j-r)){
                    index[j-1] = index[j-1] + 1
                    index[j:r] = index[j-1] + 1:(r-j+1) 
                    }
                }
            }
        }
    out = matrix(0, nrow=n, ncol=m)
    for (i in 1:r)
        out[,i] = rep(s, each=l*max(s)^(r-i))

    for (i in (r+1):m)

    for (col in 3:m){
           
        }
    return (out)
    }

oa(3, 4, 2, 1)

