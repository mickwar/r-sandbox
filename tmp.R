source("./files/repos/mwBASE/R/work_mcmc.R")

set.seed(1)
n = 100

gamma0 = 0.5
gamma1 = 0.75
beta = 3.0

x = runif(n, -2, 4)
y = rgamma(n, exp(gamma0 + gamma1*x), beta)

plot(x, y, pch = 16)

