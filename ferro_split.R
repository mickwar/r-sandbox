library(mwBASE)

set.seed(2)
### Simulated example
R = 10  # number of replicates
S = 9   # number of seasons per replicate
n = 90  # number of observations per season
#n = 810
#n = 400
theta = rnorm(R, 0.1, 0.05)
y = NULL
for (j in 1:R){
    ww = matrix(-1/log(runif(n*S)), n, S)
    x = matrix(0, n, S)
    x[1,] = ww[1,] / theta[j]
    for (i in 2:n)
        x[i,] = apply(rbind((1-theta[j])*x[i-1,], ww[i,]), 2, max)
    y = cbind(y, c(x))
    }



#uu = seq(0.5, 14.5, by = 1.0)
uu = seq(quantile(y, 0.85), quantile(y, 0.97), length = 16)
#uu = c(uu, seq(quantile(wooster, 0.95), quantile(wooster, 0.999), length = 16))
mm = double(length(uu))
vv = matrix(0, length(uu), 2)
mm.p = double(length(uu))
vv.p = matrix(0, length(uu), 2)
a.vec = double(length(uu))
u.vec = double(length(uu))
TC.vec = double(length(uu))

matplot(y, type = 'l')
abline(h = uu, lty = 2)

for (j in 1:length(uu)){
    u = uu[j]
    exceedance = apply(y, 2, function(x) which(x > u))

    Tu = lapply(exceedance, diff)
    N = sapply(exceedance, length)
    m1 = sapply(Tu, function(x) sum(x == 1))

    dat = list("y" = y, "exceed" = exceedance, 
        "Tu" = Tu, "N" = N, "m1" = m1, "up" = mean(y <= u))

    u.vec[j] = dat$up

    calc.post = function(x, param){
        R = length(x$exceed)
        n = NROW(x$y)
        
        # theta_i, p_i for each ``break''
        theta_i = param[1:R]
        p_i = param[R + (1:R)]

        theta = param[2*R + 1]
        p = param[2*R + 2]
        nu = param[2*R + 3]
        tau = param[2*R + 4]

        if (any(theta_i <= 0 | theta_i > 1))
            return (-Inf)
        if (any(p_i <= 0 | p_i >= 1))
            return (-Inf)
        if (theta <= 0 || theta > 1)
            return (-Inf)
        if (p <= 0 || p > 1)
            return (-Inf)
        if (nu <= 0)
            return (-Inf)
        if (tau <= 0)
            return (-Inf)

        m1 = x$m1

        # Likelihood (from Eq. 3)
        out = sum(m1 * log(1 - theta_i*p_i^theta_i) + 
            ifelse(N - 1 - m1 >= 0, N - 1 - m1, 0) * (log(theta_i) + log(1-p_i^theta_i)) + 
            theta_i*log(p_i) * sapply(x$Tu, function(x) sum(x - 1)))

        # Priors
        out = out + sum(dbeta(theta_i, theta*nu, (1-theta)*nu, log = TRUE))
        out = out + sum(dbeta(p_i, p*tau, (1-p)*tau, log = TRUE))
        out = out + dbeta(theta, 1/2, 1, log = TRUE)
        out = out + dbeta(p, dat$up*100, (1-dat$up)*100, log = TRUE)
        out = out + dgamma(nu, 1, 1/10, log = TRUE)
        out = out + dgamma(tau, 1, 1/10, log = TRUE)
#       out = out - log(nu)
#       out = out - log(tau)

#       (m1 * log(1 - theta_i*p_i^theta_i) + 
#           ifelse(N - 1 - m1 >= 0, N - 1 - m1, 0) * (log(theta_i) + log(1-p_i^theta_i)) + 
#           theta_i*log(p_i) * sapply(Tu, function(x) sum(x - 1)))

#       (dbeta(theta_i, theta*nu, (1-theta)*nu, log = TRUE))
#       (dbeta(p_i, p*tau, (1-p)*tau, log = TRUE))
#       dbeta(theta, 1/2, 1, log = TRUE)
#       dbeta(p, 1, 1, log = TRUE)
#       -log(nu)
#       -log(tau)

        return (out)
        }
    out = mcmc_sampler(dat, calc.post, 2*R + 4, nburn = 50000, nmcmc = 30000,
        groups = list(1:R, R + (1:R), 2*R + (1:4)))
#       groups = 0)
    mm[j] = mean(out$param[,2*R+1])
    vv[j,] = quantile(out$param[,2*R+1], c(0.025, 0.975))
    mm.p[j] = mean(out$param[,2*R+2])
    vv.p[j,] = quantile(out$param[,2*R+2], c(0.025, 0.975))
    a.vec[j] = mean(out$accept)

#    plot(0, type='n', xlim = c(0, 1), ylim = c(0, 50), bty = 'n')
#    for (i in 1:R)
#        lines(density(out$param[,i]), col = i, lty = i)
#
#    plot(0, type='n', xlim = c(0, 1), ylim = c(0, 50), bty = 'n')
#    for (i in 1:R)
#        lines(density(out$param[,R+i]), col = i, lty = i)
#
#    plot(density(out$param[,2*R+1]))    # theta
#    plot(density(out$param[,2*R+2]))    # p
#    plot(density(out$param[,2*R+3]))    # nu
#    plot(density(out$param[,2*R+4]))    # tau
#
#    plot(density(out$param[,5]))




#   ### Intervals estimator (bootstrapping)
#   ex.time = which(wooster > u)
#   #   ints.hat[j] = min(1, (2*sum(Tu - 1)^2) / ((N-1)*sum((Tu - 1)*(Tu - 2))))
#   ints.hat[j] = intervals.est(Tu)
#   C = floor(ints.hat[j]*N)+1
#   C = min(C, N-1)
#   tmp = sort(Tu, decreasing = TRUE)
#   T_C = tmp[C]
#   while (!(tmp[C-1] > T_C) && (C > 1)){
#       C = C - 1
#       T_C = tmp[C]
#       }

#   inter.Clust = Tu[Tu > T_C]

#   i_j = which(Tu > T_C)
#   i_j = c(0, i_j, N)
#   ind.seq = rep(list(NULL), C)
#   intra.Clust = rep(list(NULL), C)
#   nice.S = rep(list(NULL), C)
#   nice.C = rep(list(NULL), C)

#   for (k in 2:(C+1)){
#       ind.seq[[k-1]] = seq(i_j[k-1]+1, i_j[k]-1)
#       if (i_j[k-1]+1 == i_j[k]){
#   #       nice.T[[j-1]] = NULL
#       } else {
#           intra.Clust[[k-1]] = Tu[seq(i_j[k-1]+1, i_j[k]-1)]
#           }
#       nice.S[[k-1]] = ex.time[seq(i_j[k-1]+1, i_j[k])]
#       nice.C[[k-1]] = wooster[nice.S[[k-1]]]
#       }


#   ### Bootstrap
#   theta.vec = double(10000)
#   for (i in 1:length(theta.vec)){

#       samp.inter = sample(C-1, replace = TRUE)
#       samp.intra = sample(C, replace = TRUE)

#       tmp = c(inter.Clust[samp.inter], unlist(intra.Clust[samp.intra]))
#       theta.vec[i] = min(1, (2*sum(tmp - 1)^2) / ((N-1)*sum((tmp - 1)*(tmp - 2))))
#       theta.vec[i] = intervals.est(tmp)
#       }

#   ints.mm[j] = mean(theta.vec)
#   ints.vv[j,] = quantile(theta.vec, c(0.025, 0.975))

#   TC.vec[j] = T_C

    ux = seq(min(uu), max(uu), length = 7)
    N.vec = double(length(ux))
    for (i in 1:length(ux))
        N.vec[i] = sum(y > ux[i])

    plot(uu[1:j], mm[1:j], bty = 'n', xlim = range(uu), ylim = c(0, 1), type = 'b', pch = 16)
    lines(uu[1:j], vv[1:j,1], lty = 3)
    lines(uu[1:j], vv[1:j,2], lty = 3)
    points(uu[1:j], mm.p[1:j], col = 'green', pch = 16)
    lines(uu[1:j], vv.p[1:j,1], col = 'green', lty = 3)
    lines(uu[1:j], vv.p[1:j,2], col = 'green', lty = 3)
    points(uu[1:j], u.vec[1:j], col = 'blue', pch = 15)
    axis(3,  at = ux, labels = 1-round(N.vec / (n*R*S), 3))
    axis(3,  at = ux, labels = N.vec, line = 0.75, lty = 0)
#   axis(3,  uu[1:j], labels = TC.vec[1:j], line = 1.5, lty = 0)
#   mtext("T_C", 3, line = 2.6, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
    mtext("# >u", 3, line = 1.9, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
    mtext("1-p", 3, line = 1.2, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
    if (exists("theta"))
        abline(h = mean(theta), lty = 2, col = 'blue')

    }

# system.time({
#     for (i in 1:1000){
#         rnorm(184, 0, 1)
#         }
#     })


