set.seed(1)
n = 90
R = 100
x = matrix(0, n, R)
x[1,] = rnorm(R)
for (i in 2:n)
    x[i,] = 0.95*x[i-1,] + rnorm(R)
y = x
matplot(x, type = 'l')
nburn = 30000
nmcmc = 20000
chain_init = c(0.5, 0.5)

uni_theta = function(x, u, ord, prior, nburn = 30000, nmcmc = 20000, chain_init = c(0.5, 0.5)){
    require(mwBASE)

    n = NROW(x)
    R = NCOL(x)

    # ``ord'' is for re-ordering the columns when vec'ing
    # (Could result in different interexceedance times)
    if (missing(ord))
        ord = 1:R

    # Vec the matrix (if applicable)
    if (R > 1)
        x = c(x[,ord])

    # Default threshold to the 0.90 quantile
    if (missing(u))
        u = quantile(x, 0.90)

    
    exceed = which(x > u)
    Tu = diff(exceed)
    N = length(exceed)
    m1 = sum(Tu == 1)
    up = mean(x <= u)

    # If multiple realizations, adjust m1 for the possibility of truncating
    # clusters of extremes
    if (R > 1)
        m1 = m1 + (sum(which(x > u) %% n == 0) + sum(which(x > u) %% n == 1)) / 2

    if (missing(prior))
        prior = list("theta_a" = 1/2, "theta_b" = 1,
            "p_a" = up*10, "p_b" = (1-up)*10)

    calc.post = function(x, param){
        theta = param[1]
        p = param[2]
        if (theta <= 0 || theta > 1)
            return (-Inf)
        if (p <= 0 || p >= 1)
            return (-Inf)

        # Likelihood (Ferro, Eq. 3)
        out = m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
            theta*log(p)*sum(x - 1)

        # Priors
        out = out + dbeta(theta, prior$theta_a, prior$theta_b, log = TRUE)
        out = out + dbeta(p, prior$p_a, prior$p_b, log = TRUE)
        return (out)
        }
    mcmc_out = mcmc_sampler(Tu, calc.post, 2, nburn = nburn, nmcmc = nmcmc, chain_init = chain_init)

    theta.hat = mean(mcmc_out$param[,1])


    ### Intervals estimator (bootstrapping)
    C = floor(theta.hat*N)+1
    C = min(C, N-1)
    tmp = sort(Tu, decreasing = TRUE)
    T_C = tmp[C]
    while (!(tmp[C-1] > T_C) && (C > 1)){
        C = C - 1
        T_C = tmp[C]
        }

    # The set of independent intercluster times
    inter.Clust = Tu[Tu > T_C]

    i_j = which(Tu > T_C)
    i_j = c(0, i_j, N)
    ind.seq = rep(list(NULL), C)
    intra.Clust = rep(list(NULL), C)    # The interexceedance times within each cluster
    nice.S = rep(list(NULL), C)     # List of independent clusters, marking when exceedances occur
    nice.C = rep(list(NULL), C)     # The observed value at the exceedance times

    for (k in 2:(C+1)){
        ind.seq[[k-1]] = seq(i_j[k-1]+1, i_j[k]-1)
        if (i_j[k-1]+1 == i_j[k]){
    #       nice.T[[j-1]] = NULL
        } else {
            intra.Clust[[k-1]] = Tu[seq(i_j[k-1]+1, i_j[k]-1)]
            }
        nice.S[[k-1]] = exceed[seq(i_j[k-1]+1, i_j[k])]
        nice.C[[k-1]] = x[nice.S[[k-1]]]
        }

    nice.S[[1]] %% n

    lapply(nice.S, function(x) x %% n)

    which(1:190 %% n == 1)

#   plot(x, bty = 'n', xlim = c(350, 650))
#   abline(h = u, lty = 2, col = 'gray50')
#   for (i in 1:95){
#       points(nice.S[[i]], x[nice.S[[i]]], col = 'red')
#       points(nice.S[[i]][1], x[nice.S[[i]][1]], col = 'green')
#       }



    }


