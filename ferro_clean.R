### Simulated example
intervals.est = function(Tu){
    if (length(Tu) == 0)
        return (1)
    if (max(Tu) <= 2){
        out = min(1, 2*(sum(Tu))^2 / ( length(Tu) * sum(Tu^2) ) )
    } else {
        out = min(1, 2*(sum(Tu - 1))^2 / (length(Tu)*sum((Tu-1)*(Tu-2))))
        }
    return (out)
    }
# bayes.est = function(Tu){
#     # theta ~ Beta(a, b)
#     pr.a = 1
#     pr.b = N
# 
#     # q = p^theta ~ Beta(c, d)
#     pr.c = 1
#     pr.d = 1
# 
#     m1 = sum(Tu == 1)
#     N = length(Tu) + 1
# 
#     k = 0:m1
#     v = lchoose(m1, k) + lbeta(k + pr.a + N - 1 - m1 + 1, pr.b) + lbeta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)
# 
#     lchoose(m1, k)
# 
#     lbeta(k + pr.a + N - 1 - m1 + 1, pr.b)
#     lbeta(0 + pr.a + N - 1 - m1 + 1, pr.b)
# 
#     lbeta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)
#     lbeta(0 + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)
# 
# 
# #   v = v + (pr.d + N - 1 - m1)*log(0 + pr.c + sum(Tu - 1))
# #   v = v - lchoose(m1, round(m1/2))
#     v = v - lbeta(0 + pr.a + N - 1 - m1 + 1, pr.b)
#     v = v - lbeta(0 + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)
# 
# num = sum(((-1)^k)*exp(v))
# 
# (((-1)^k)*exp(v))
# 
#     logavec = sort(v[(k %% 2) == 0], decreasing = TRUE)
#     logbvec = sort(v[(k %% 2) == 1], decreasing = TRUE)
# 
#     lsa = logavec[1] + log(1 + sum(exp(logavec[-1] - logavec[1])))
#     lsb = logbvec[1] + log(1 + sum(exp(logbvec[-1] - logbvec[1])))
# 
#     num = lsa + log(1 - exp(lsb - lsa))
# 
#     v = lchoose(m1, k) + lbeta(k + pr.a + N - 1 - m1, pr.b) + lbeta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)
# 
# #   v = v + (pr.d + N - 1 - m1)*log(0 + pr.c + sum(Tu - 1))
# #   v = v - lchoose(m1, round(m1/2))
#     v = v - lbeta(0 + pr.a + N - 1 - m1 + 1, pr.b)
#     v = v - lbeta(0 + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)
# 
# den = sum(((-1)^k)*exp(v))
# 
# num/den
# 
# (((-1)^k)*exp(v))
# 
# 
# ((N - m1 + 1) / (N + 3)) / ((sum(Tu) - 2*N + 3 + m1) / (sum(Tu) - N + 2))
# 
# 
# B = 10000
# rs = rbeta(B, N - m1 + 1, m1 + 2)
# qs = rbeta(B, sum(Tu) -2*N+2+m1+1,N-m1)
# 
# mean(rs/qs)
# 
# plot(density(rs/qs))
# lines(density(out$params), col = 'red')
# 
#     logavec = sort(v[(k %% 2) == 0], decreasing = TRUE)
#     logbvec = sort(v[(k %% 2) == 1], decreasing = TRUE)
# 
#     lsa = logavec[1] + log(1 + sum(exp(logavec[-1] - logavec[1])))
#     lsb = logbvec[1] + log(1 + sum(exp(logbvec[-1] - logbvec[1])))
# 
#     den = lsa + log(1 - exp(lsb - lsa))
# 
#     exp(num - den)
# 
#     sum(choose(m1, k)*((-1)^k)*beta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)*beta(k + pr.a + N - 1 - m1 + 1, pr.b)) / 
#     sum(choose(m1, k)*((-1)^k)*beta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)*beta(k + pr.a + N - 1 - m1, pr.b))
# 
#     sum(choose(m1, k)*((-1)^k)*beta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)*beta(k + pr.a + N - 1 - m1 + 1, pr.b))
#     sum(choose(m1, k)*((-1)^k)*beta(k + pr.c + sum(Tu - 1), pr.d + N - 1 - m1)*beta(k + pr.a + N - 1 - m1, pr.b))
# 
#     }



library(mwBASE)


# theta = 0.5
# ww = -1/log(runif(n))
# ksi = double(n)
# ksi[1] = ww[1] / theta
# for (i in 2:n)
#     ksi[i] = max((1-theta)*ksi[i-1], ww[i])
# 
# plot(log(ksi), xlim = c(1, 200), type='o')


theta.vec = seq(0.05, 0.95, length = 10)
mm = double(length(theta.vec))
vv = matrix(0, length(theta.vec), 2)
ints.mat = matrix(0, 1, length(theta.vec))
for (j in 1:length(theta.vec)){
    set.seed(2)
    theta = theta.vec[j]
    n = 10000
    ww = -1/log(runif(n))
    x = double(n)
    x[1] = ww[1] / theta
    for (i in 2:n)
        x[i] = max((1-theta)*x[i-1], ww[i])
    plot(log(x), xlim = c(0, 100))
    plot((x) )

    u = 250

    Tu = diff(which(x > u))

    N = sum(x > u)
    m1 = sum(Tu == 1)

    calc.post = function(z, param){
        theta = param[1]
        p = param[2]
        if (theta <= 0 || theta > 1)
            return (-Inf)
        if (p <= 0 || p >= 1)
            return (-Inf)

        # Likelihood (from Eq. 3)
        out = m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
            theta*log(p)*sum(Tu - 1)

        # Priors
#       out = out + dbeta(p, 1/2, 1/2, log = TRUE)
        return (out)
        }
    out = mcmc_sampler(z, calc.post, 2, nburn = 20000, chain_init = c(theta, 0.5))
    mm[j] = mean(out$param[,1])
    vv[j,] = quantile(out$param[,1], c(0.025, 0.975))

    ints.mat[1,j] = intervals.est(Tu)

#   ### "Bootstrap" sample for the intervals estimator
#   ### Adds about 10 seconds for each iteration j
#   for (k in 1:NROW(ints.mat)){
#       ww = -1/log(runif(n))
#       x = double(n)
#       x[1] = ww[1] / theta
#       for (i in 2:n)
#           x[i] = max((1-theta)*x[i-1], ww[i])

#       Tu = diff(which(x > u))
#       ints.mat[k,j] = intervals.est(Tu)
#       }
    }


offset = 0.0075
qq = apply(ints.mat, 2, quantile, c(0.025, 0.975))
plot(theta.vec, mm, ylim = range(qq), xlab = expression("True" ~ theta),
    ylab = expression("Posterior" ~ theta))
points(theta.vec, vv[,1], pch="-")
points(theta.vec, vv[,2], pch="-")
segments(x0 = theta.vec, y0 = vv[,1], x1 = theta.vec, y1 = vv[,2])
points(theta.vec, mm, pch = 16)
points(theta.vec+offset, colMeans(ints.mat), col = 'red')
points(theta.vec+offset, qq[1,], pch="-", col = 'red')
points(theta.vec+offset, qq[2,], pch="-", col = 'red')
segments(x0 = theta.vec+offset, y0 = qq[1,], x1 = theta.vec+offset, y1 = qq[2,], col = 'red')
abline(0, 1)



### Wooster data set
#  library(ismev)
#  data(wooster)
#  wooster = -wooster
x = read.csv("~/files/repos/data/wooster_1893_2014.csv")

# Get rid of the extra space
y = gsub(" ", "", x[,6]) # end at 2014
y = gsub(" ", "", x[1:32324,6]) # end at 1987
y = gsub(" ", "", x[1:37803,6]) # end at 2002

y = y[-which(y == ".")]  # remove missing data
y = as.numeric(y)

tail(x, 1)
x[32324,]
x[37803,]

wooster = -y

n = 1000
theta = 0.05
ww = -1/log(runif(n))
y = double(n)
y[1] = ww[1] / theta
for (i in 2:n)
    y[i] = max((1-theta)*y[i-1], ww[i])
wooster = y


plot(tail(wooster, 2000), pch = 16)

n = length(wooster)
#uu = seq(0.5, 14.5, by = 1.0)
uu = seq(quantile(wooster, 0.9), quantile(wooster, 0.99), length = 16)
mm = double(length(uu))
vv = matrix(0, length(uu), 2)
a.vec = double(length(uu))
ints.hat = double(length(uu))
ints.mm = double(length(uu))
ints.vv = matrix(0, length(uu), 2)
for (j in 1:length(uu)){
    u = uu[j]
    Tu = diff(which(wooster > u))
    N = sum(wooster > u)
    m1 = sum(Tu == 1)

    calc.post = function(x, param){
        theta = param[1]
        p = param[2]
        if (theta <= 0 || theta > 1)
            return (-Inf)
        if (p <= 0 || p >= 1)
            return (-Inf)

        # Likelihood (from Eq. 3)
        out = m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
            theta*log(p)*sum(x - 1)
        m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
            theta*log(p)*sum(x - 1)
        (N - 1)*log(1-p^theta) + (N - 1)*log(theta) + theta*log(p)*sum(x - 1)

#   theta*p^(x*theta) - theta*p^(theta*(x+1))

        # Priors
#       out = out + dbeta(theta, 1/2, 1/2, log = TRUE)
        out = out + dbeta(theta, 1/2, 1/20, log = TRUE)
        out = out + dbeta(p, 1/2, 1/2, log = TRUE)
        return (out)
        }
    out = mcmc_sampler(Tu, calc.post, 2, nburn = 50000, nmcmc = 30000, chain_init = c(0.5, 0.5))
    mm[j] = mean(out$param[,1])
    vv[j,] = quantile(out$param[,1], c(0.025, 0.975))
    a.vec[j] = mean(out$accept)


    ### Intervals estimator (bootstrapping)
    ex.time = which(wooster > u)
#   ints.hat[j] = min(1, (2*sum(Tu - 1)^2) / ((N-1)*sum((Tu - 1)*(Tu - 2))))
    ints.hat[j] = intervals.est(Tu)
    C = floor(ints.hat[j]*N)+1
    C = min(C, N-1)
    tmp = sort(Tu, decreasing = TRUE)
    T_C = tmp[C]
    while (!(tmp[C-1] > T_C)){
        C = C - 1
        T_C = tmp[C]
        }

    inter.Clust = Tu[Tu > T_C]

    i_j = which(Tu > T_C)
    i_j = c(0, i_j, N)
    ind.seq = rep(list(NULL), C)
    intra.Clust = rep(list(NULL), C)
    nice.S = rep(list(NULL), C)
    nice.C = rep(list(NULL), C)

    for (k in 2:(C+1)){
        ind.seq[[k-1]] = seq(i_j[k-1]+1, i_j[k]-1)
        if (i_j[k-1]+1 == i_j[k]){
    #       nice.T[[j-1]] = NULL
        } else {
            intra.Clust[[k-1]] = Tu[seq(i_j[k-1]+1, i_j[k]-1)]
            }
        nice.S[[k-1]] = ex.time[seq(i_j[k-1]+1, i_j[k])]
        nice.C[[k-1]] = wooster[nice.S[[k-1]]]
        }


    ### Bootstrap
    theta.vec = double(10000)
    for (i in 1:length(theta.vec)){

        samp.inter = sample(C-1, replace = TRUE)
        samp.intra = sample(C, replace = TRUE)

        tmp = c(inter.Clust[samp.inter], unlist(intra.Clust[samp.intra]))
#       theta.vec[i] = min(1, (2*sum(tmp - 1)^2) / ((N-1)*sum((tmp - 1)*(tmp - 2))))
        theta.vec[i] = intervals.est(tmp)
        }

    ints.mm[j] = mean(theta.vec)
    ints.vv[j,] = quantile(theta.vec, c(0.025, 0.975))

    ux = seq(min(uu), max(uu), length = 7)
    N.vec = double(length(ux))
    for (i in 1:length(ux))
        N.vec[i] = sum(wooster > ux[i])

    plot(uu[1:j], mm[1:j], bty = 'n', xlim = range(uu), ylim = c(0, 1), type = 'b', pch = 16)
    lines(uu[1:j], vv[1:j,1], lty = 3)
    lines(uu[1:j], vv[1:j,2], lty = 3)
    points(uu[1:j], ints.hat[1:j], col = 'green')
    points(uu[1:j], ints.mm[1:j], col = 'red', pch = 16, type='b')
    lines(uu[1:j], ints.vv[1:j,1], lty = 3, col = 'red')
    lines(uu[1:j], ints.vv[1:j,2], lty = 3, col = 'red')
    axis(3,  at = ux, labels = N.vec)
    abline(h = theta, lty = 2, col = 'blue')

    }




range(a.vec)

ints = sapply(uu, function(x) intervals.est(diff(which(wooster > x))))

ux = seq(-10, 20, by = 2)
N.vec = double(length(ux))
for (i in 1:length(ux))
    N.vec[i] = sum(wooster > ux[i])

plot(uu, mm, ylim = c(0, 1), xlab = "Threshold", ylab = "Extremal Index",
    bty = 'n', type='b')
lines(uu, vv[,1], lty=2)
lines(uu, vv[,2], lty=2)
axis(3,  at = ux, labels = N.vec)
points(uu, ints, col = 'red')

# plot(uu, apply(ints.boot, 2, mean, na.rm = TRUE), ylim = c(0, 1), xlab = "Threshold",
#     ylab = "Extremal Index", bty = 'n', type='b')
# lines(uu, apply(ints.boot, 2, quantile, 0.025, na.rm = TRUE), lty=2)
# lines(uu, apply(ints.boot, 2, quantile, 0.975, na.rm = TRUE), lty=2)
# axis(3,  at = ux, labels = N.vec)

