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
R = 10
n = 500
#n = 810
#n = 400
theta = 0.10
ww = matrix(-1/log(runif(n*R)), n, R)
y = matrix(0, n, R)
y[1,] = ww[1,] / theta
for (i in 2:n)
    y[i,] = apply(rbind((1-theta)*y[i-1,], ww[i,]), 2, max)

K = 16
uu = seq(quantile(y, 0.85), quantile(y, 0.97), length = K)

fer.m = double(K)
fer.v = matrix(0, K, 2)

suv.m = double(K)
suv.v = matrix(0, K, 2)

u.vec = double(length(uu))

for (j in 11:K){
    up = mean(y <= uu[j])
    u.vec[j] = up
    

    ferro = theta_hier(y, uu[j], likelihood = "ferro",
        nburn = 40000, nmcmc = 20000)
    fer.m[j] = mean(ferro$mcmc[,NCOL(ferro$mcmc)])
    fer.v[j,] = quantile(ferro$mcmc[,NCOL(ferro$mcmc)], c(0.025, 0.975))

    suveges = theta_hier(y, uu[j], likelihood = "suveges",
        nburn = 40000, nmcmc = 20000)
    suv.m[j] = mean(suveges$mcmc[,NCOL(suveges$mcmc)])
    suv.v[j,] = quantile(suveges$mcmc[,NCOL(suveges$mcmc)], c(0.025, 0.975))

    plot(0, type = 'n', xlim = range(uu), ylim = c(0, 1), bty = 'n',
        xlab = "Threshold", ylab = expression(theta))
    points(uu[1:j], fer.m[1:j], type = 'b', pch = 16, col = 'blue')
    points(uu[1:j], suv.m[1:j], type = 'b', pch = 16, col = 'orange')

    lines(uu[1:j], fer.v[1:j,1], col = 'lightblue')
    lines(uu[1:j], fer.v[1:j,2], col = 'lightblue')
    lines(uu[1:j], suv.v[1:j,1], col = 'yellow')
    lines(uu[1:j], suv.v[1:j,2], col = 'yellow')
    if (exists("theta"))
        abline(h = theta, lty = 2)
    }

