library(mwBASE)
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


### Wooster data set
#  library(ismev)
#  data(wooster)
#  wooster = -wooster
x = read.csv("~/files/repos/data/wooster_1893_2014.csv")

# Get rid of the extra space
y = gsub(" ", "", x[,6]) # end at 2014
y = gsub(" ", "", x[1:32324,6]) # and at 1987
y = gsub(" ", "", x[1:37803,6]) # end at 2002
y = gsub(" ", "", x[1:38168,6]) # end at 2003

y = y[-which(y == ".")]  # remove missing data
y = as.numeric(y)

wooster = -y

# ### Simulated example
# n = 90*9*10
# #n = 810
# #n = 400
# theta = 0.50
# ww = -1/log(runif(n))
# y = double(n)
# y[1] = ww[1] / theta
# for (i in 2:n)
#     y[i] = max((1-theta)*y[i-1], ww[i])
# wooster = y

set.seed(1)
R = 10*9
n = 90
#n = 810
#n = 400
theta = 0.15
ww = matrix(-1/log(runif(n*R)), n, R)
y = matrix(0, n, R)
y[1,] = ww[1,] / theta
for (i in 2:n)
    y[i,] = apply(rbind((1-theta)*y[i-1,], ww[i,]), 2, max)
wooster = y

set.seed(6)
wooster = c(wooster[,sample(90, 90)])
n = length(wooster)
R = 1

#n = length(wooster)
#uu = seq(0.5, 14.5, by = 1.0)
uu = seq(quantile(wooster, 0.85), quantile(wooster, 0.97), length = 16)
#uu = c(uu, seq(quantile(wooster, 0.95), quantile(wooster, 0.999), length = 16))
mm = double(length(uu))
vv = matrix(0, length(uu), 2)
mm.p = double(length(uu))
vv.p = matrix(0, length(uu), 2)
a.vec = double(length(uu))
ints.hat = double(length(uu))
ints.mm = double(length(uu))
ints.vv = matrix(0, length(uu), 2)
TC.vec = double(length(uu))
u.vec = double(length(uu))

plot(wooster, pch = 16, type='l', bty = 'n')
#plot(tail(wooster, 2000), pch = 16, type='l', bty = 'n')
abline(h = uu, lty = 2)

for (j in 1:length(uu)){
    u = uu[j]
    Tu = diff(which(wooster > u))
    N = sum(wooster > u)
    m1 = sum(Tu == 1)
    up = mean(wooster <= u)

    m1 = m1 + (sum(which(wooster > u) %% 90 == 0) + sum(which(wooster > u) %% 90 == 1)) / 2
# which(wooster > u) %% 90
# which(which(wooster > u) %% 90 == 0 | which(wooster > u) %% 90 == 1)
# 
# 91 %% 90

#     u = uu[j]
#     up = mean(wooster <= u)
#     exceedance = apply(wooster, 2, function(x) which(x > u))
# 
#     Tu = unlist(lapply(exceedance, diff))
# #   N = length(Tu) + 1
#     N = sum(sapply(exceedance, length))
# #   NM = max(1, sum(sapply(exceedance, length) == 0))
# #   NM = sum(sapply(exceedance, length) >= 2)
#     NM = sum(sapply(exceedance, length) != 1)
#     m1 = sum(Tu == 1)

    u.vec[j] = up

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

#       m1 * log(1 - theta*p^theta) + (N - 1 - m1)*(log(theta) + log(1-p^theta)) +
#           theta*log(p)*sum(x - 1)
#       (N - 1)*log(1-p^theta) + (N - 1)*log(theta) + theta*log(p)*sum(x - 1)

        # Priors
        out = out + dbeta(theta, 1/2, 1, log = TRUE)
#       out = out + dbeta(theta, 1/2, 1/2+2*R, log = TRUE)
#       out = out + dbeta(theta, 1/(1+R), R/2, log = TRUE)
#       out = out + dbeta(theta, 1/(1+R), 1/2, log = TRUE)
#       out = out + dbeta(theta, 1/2, 1/20, log = TRUE)
 
#       out = out + dbeta(p, 1/2, 1/2, log = TRUE)
#       out = out + dbeta(p, up*2, (1-up)*2, log = TRUE)
#       out = out + dbeta(p, up*100, (1-up)*100, log = TRUE)
        return (out)
        }
    out = mcmc_sampler(Tu, calc.post, 2, nburn = 30000, nmcmc = 20000, chain_init = c(0.5, 0.5))
    mm[j] = mean(out$param[,1])
    vv[j,] = quantile(out$param[,1], c(0.025, 0.975))
    a.vec[j] = mean(out$accept)

    mm.p[j] = mean(out$param[,2])
    vv.p[j,] = quantile(out$param[,2], c(0.025, 0.975))


    ### Intervals estimator (bootstrapping)
    ex.time = which(wooster > u)
#   ints.hat[j] = min(1, (2*sum(Tu - 1)^2) / ((N-1)*sum((Tu - 1)*(Tu - 2))))
    ints.hat[j] = intervals.est(Tu)
    C = floor(ints.hat[j]*N)+1
    C = min(C, N-1)
    tmp = sort(Tu, decreasing = TRUE)
    T_C = tmp[C]
    while (!(tmp[C-1] > T_C) && (C > 1)){
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

    TC.vec[j] = T_C

    ux = seq(min(uu), max(uu), length = 7)
    N.vec = double(length(ux))
    for (i in 1:length(ux))
        N.vec[i] = sum(wooster > ux[i])

    plot(uu[1:j], mm[1:j], bty = 'n', xlim = range(uu), ylim = c(0, 1), type = 'b', pch = 16)
    lines(uu[1:j], vv[1:j,1], lty = 3)
    lines(uu[1:j], vv[1:j,2], lty = 3)

    points(uu[1:j], mm.p[1:j], type = 'b', pch = 16, col = 'green')
    lines(uu[1:j], vv.p[1:j,1], lty = 3, col = 'green')
    lines(uu[1:j], vv.p[1:j,2], lty = 3, col = 'green')
    points(uu[1:j], ints.hat[1:j], col = 'green')
    points(uu[1:j], ints.mm[1:j], col = 'red', pch = 16, type='b')
    lines(uu[1:j], ints.vv[1:j,1], lty = 3, col = 'red')
    lines(uu[1:j], ints.vv[1:j,2], lty = 3, col = 'red')
    axis(3,  at = ux, labels = 1-round(N.vec / (n*R), 3))
    axis(3,  at = ux, labels = N.vec, line = 0.75, lty = 0)
    axis(3,  uu[1:j], labels = TC.vec[1:j], line = 1.5, lty = 0)
    mtext("T_C", 3, line = 2.6, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
    mtext("# >u", 3, line = 1.9, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
    mtext("1-p", 3, line = 1.2, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
    points(uu[1:j], u.vec[1:j], col = 'blue', pch = 15)
    if (exists("theta"))
        abline(h = theta, lty = 2, col = 'blue')

    }


pdf("~/theta_0.85.pdf", width = 9, height = 9)
plot(uu[1:j], mm[1:j], bty = 'n', xlim = range(uu), ylim = c(0, 1), type = 'b', pch = 16)
lines(uu[1:j], vv[1:j,1], lty = 3)
lines(uu[1:j], vv[1:j,2], lty = 3)
points(uu[1:j], mm.p[1:j], type = 'b', pch = 16, col = 'green')
lines(uu[1:j], vv.p[1:j,1], lty = 3, col = 'green')
lines(uu[1:j], vv.p[1:j,2], lty = 3, col = 'green')
points(uu[1:j], ints.hat[1:j], col = 'green')
points(uu[1:j], ints.mm[1:j], col = 'red', pch = 16, type='b')
lines(uu[1:j], ints.vv[1:j,1], lty = 3, col = 'red')
lines(uu[1:j], ints.vv[1:j,2], lty = 3, col = 'red')
axis(3,  at = ux, labels = 1-round(N.vec / (n*R), 3))
axis(3,  at = ux, labels = N.vec, line = 0.75, lty = 0)
axis(3,  uu[1:j], labels = TC.vec[1:j], line = 1.5, lty = 0)
mtext("T_C", 3, line = 2.6, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
mtext("# >u", 3, line = 1.9, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
mtext("1-p", 3, line = 1.2, at = uu[1] - (uu[2] - uu[1]), cex = 0.75)
points(uu[1:j], u.vec[1:j], col = 'blue', pch = 15)
if (exists("theta"))
    abline(h = theta, lty = 2, col = 'blue')
dev.off()

lapply(nice.S, function(x) x %% 90)
