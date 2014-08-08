R.x = function(x, rho){
    n = nrow(x)
    px = ncol(x)
    out = matrix(0, n, n) + diag(n)
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            out[i,j] = 1
            for (k in 1:px){
                out[i,j] = out[i,j] * rho[k]^(4*(x[i,k]-x[j,k])^2)
                }
            out[j,i] = out[i,j]
            }
        }
    return (out)
    }

R.cross = function(x, y, rho){
    n = nrow(x)
    m = nrow(y)
    px = ncol(x)
    out = matrix(0, n, m)
    for (i in 1:n){
        for (j in 1:m){
            out[i,j] = 1
            for (k in 1:px){
                out[i,j] = out[i,j] * rho[k]^(4*(x[i,k]-y[j,k])^2)
                }
            }
        }
    return (out)
    }

diag2 = function(...){
    x = list(...)
    n = 0
    k = 0
    for (i in 1:length(x)){
        n = n + nrow(x[[i]])
        k = k + ncol(x[[i]])
        }
    out = matrix(0, n, k)
    offset.x = 0
    offset.y = 0
    for (i in 1:length(x)){
        out[1:nrow(x[[i]]) + offset.x, 1:ncol(x[[i]]) + offset.y] = x[[i]]
        offset.x = offset.x + nrow(x[[i]])
        offset.y = offset.y + ncol(x[[i]])
        }
    return(out)
    }


matrix.assemble = function(grid, ...){
    # read in by row, NA to denote matrix of 0's
    x = list(...)
    y = length(x)
    grid.m = grid[1] # number of row in matrix grid
    grid.n = grid[2] # number of columns in matrix grid
    if (y != grid.m * grid.n)
        stop("Number of matrices provided does not equal grid area.")
    out.m = 0
    out.n = 0
    check.col = matrix(0, grid.m, grid.n)
    check.row = matrix(0, grid.m, grid.n)
    for (i in 1:y){
        if (all(!is.na(x[[i]]))){
            check.col[floor((i-1) / grid.n)+1, ((i-1) %% grid.n)+1]  = ncol(x[[i]])
            check.row[floor((i-1) / grid.n)+1, ((i-1) %% grid.n)+1]  = nrow(x[[i]])
            }
        }
    if (any(apply(check.row, 1, function(x) length(unique(x[x>0]))) != 1))
        stop("bad sizes")
    if (any(apply(check.col, 2, function(x) length(unique(x[x>0]))) != 1))
        stop("bad sizes")
    out.m = cumsum(c(0, apply(check.row, 1, max)))
    out.n = cumsum(c(0, apply(check.col, 2, max)))
    out = matrix(0, max(out.m), max(out.n))
    for (i in 1:y){
        if (all(!is.na(x[[i]]))){
            out[(out.m[floor((i-1) / grid.n)+1]+1):out.m[floor((i-1) / grid.n)+2],
                (out.n[((i-1) %% grid.n)+1]+1):out.n[((i-1) %% grid.n)+2]] = x[[i]]
            }
        }
    return (out)
    }


# just as fast as hpd()
# hpd.loops = function(x, prob = 0.95, precision = 1000){
#     short = Inf
#     best = 0
#     for (i in seq(0.00, 1-prob, length=precision)){
#         if (diff(quantile(x, c(i, i+prob))) < short){
#             short = diff(quantile(x, c(i, i+prob)))
#             best = i
#             }
#         }
#     return (quantile(x, c(best, best+prob)))
#     }

hpd = function(x, prob = 0.95, precision = 1000){
    range = seq(0, 1-prob, length=precision)
    range = cbind(range, range+prob)
    best = range[which.min(apply(range, 1, function(y)
        diff(quantile(x, y)))),]
    return (quantile(x, best))
    }

trailing = function(x, digits=3)
    formatC(x, digits=digits, format="f")
