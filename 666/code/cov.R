# calculate equation 7.19
# x is a list of matrices each of which as p columns
calc.M = function(x){
    k = length(x)
    # nu = degrees of freedom
    nu = unlist(lapply(x, nrow)) - 1
    p = ncol(x[[1]])
    S = lapply(x, var)
    Spl = matrix(0, p, p)
    for (i in 1:k)
        Spl = Spl + nu[i] * S[[i]]
    Spl = Spl / sum(nu)
    M = 1
    for (i in 1:k)
        M = M * det(S[[i]])^(nu[i]/2)
    M = M / det(Spl)^(sum(nu)/2)
    return (M)
    }

# equation 7.22
calc.c1 = function(x){
    k = length(x)
    nu = unlist(lapply(x, nrow)) - 1
    p = ncol(x[[1]])
    (sum(1/nu) - 1/sum(nu)) * (2*p^2 + 3*p -1) /
        (6*(p+1)*(k-1))
    }

simulate = function(n, k, p){
    x = rep(list(0), k)
    for (i in 1:k)
        x[[i]] = matrix(rnorm(n[i]*p), n[i], p)
    return (x)
    }

# number of variables
p = 5

# number of samples
k = 6

# number of observations for each sample
# must be length k
n = c(18, 46, 30, 25, 23, 32)

x = simulate(n, k, p)

M = calc.M(x)
c1 = calc.c1(x)

# u ~ chi^2, df=0.5*(k-1)*p*(p+1)
# eq 7.23
u = -2 * (1 - c1) * log(M)
pchisq(u, 0.5*(k-1)*p*(p+1), lower.tail = FALSE)
