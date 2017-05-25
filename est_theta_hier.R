library(mwEVT)

# ### Wooster data
# x = read.csv("~/files/repos/data/wooster_1893_2014.csv")
# 
# # Get rid of the extra space
# # y = gsub(" ", "", x[,6]) # end at 2014
# # y = gsub(" ", "", x[1:32324,6]) # and at 1987
# # y = gsub(" ", "", x[1:37803,6]) # end at 2002
# y = gsub(" ", "", x[1:38168,6]) # end at 2003
# 
# y = y[-which(y == ".")]  # remove missing data
# y = -as.numeric(y)
# 
# n = length(y)
# 
# uu = seq(0.5, 14.5, by = 1.0)
# K = length(uu)


### Generate dependent random variables with standard Frechet marginals
### For one long sequence
# set.seed(1)
# n = 365*100     # 50 years of data, no leap years
# 
# theta = 0.50
# W = -1/log(runif(n))
# y = double(n)
# y[1] = W[1] / theta
# for (i in 2:n)
#     y[i] = max((1-theta)*y[i-1], W[i])

### For split up sequences
set.seed(1)
R = 10*9
n = 90
#n = 810
#n = 400
theta = 0.85
ww = matrix(-1/log(runif(n*R)), n, R)
y = matrix(0, n, R)
y[1,] = ww[1,] / theta
for (i in 2:n)
    y[i,] = apply(rbind((1-theta)*y[i-1,], ww[i,]), 2, max)
y = c(y)

K = 16
uu = seq(quantile(y, 0.85), quantile(y, 0.97), length = K)

fc.m = double(K)
fc.v = matrix(0, K, 2)

fb.m = double(K)
fb.v = matrix(0, K, 2)

sc.m = double(K)
sc.v = matrix(0, K, 2)

sb.m = double(K)
sb.v = matrix(0, K, 2)

u.vec = double(length(uu))


for (j in 1:K){
    up = mean(y <= uu[j])
    u.vec[j] = up
    
    ferro_cl = theta_uni(y, uu[j], likelihood = "ferro",
        method = "classical")
    fc.m[j] = ferro_cl$theta
    fc.v[j,] = quantile(ferro_cl$bootstrap, c(0.025, 0.975))

    ferro_ba = theta_uni(y, uu[j], likelihood = "ferro",
        method = "bayesian", nburn = 40000, nmcmc = 20000)
    fb.m[j] = mean(ferro_ba$mcmc)
    fb.v[j,] = quantile(ferro_ba$mcmc, c(0.025, 0.975))

    suveges_cl = theta_uni(y, uu[j], likelihood = "suveges",
        method = "classical")
    sc.m[j] = suveges_cl$theta
    sc.v[j,] = quantile(suveges_cl$bootstrap, c(0.025, 0.975))

    suveges_ba = theta_uni(y, uu[j], likelihood = "suveges",
        method = "bayesian", nburn = 40000, nmcmc = 20000)
    sb.m[j] = mean(suveges_ba$mcmc)
    sb.v[j,] = quantile(suveges_ba$mcmc, c(0.025, 0.975))

    plot(0, type = 'n', xlim = range(uu), ylim = c(0, 1), bty = 'n')
    points(uu[1:j], fc.m[1:j], type = 'b', pch = 15, col = 'green')
    points(uu[1:j], fb.m[1:j], type = 'b', pch = 16, col = 'blue')
    points(uu[1:j], sc.m[1:j], type = 'b', pch = 15, col = 'red')
    points(uu[1:j], sb.m[1:j], type = 'b', pch = 16, col = 'orange')

    lines(uu[1:j], fc.v[1:j,1], col = 'lightgreen')
    lines(uu[1:j], fc.v[1:j,2], col = 'lightgreen')
    lines(uu[1:j], fb.v[1:j,1], col = 'lightblue')
    lines(uu[1:j], fb.v[1:j,2], col = 'lightblue')
    lines(uu[1:j], sc.v[1:j,1], col = 'pink')
    lines(uu[1:j], sc.v[1:j,2], col = 'pink')
    lines(uu[1:j], sb.v[1:j,1], col = 'yellow')
    lines(uu[1:j], sb.v[1:j,2], col = 'yellow')
    if (exists("theta"))
        abline(h = theta, lty = 2)
    }

