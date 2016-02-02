# Marin, Pudlo, Robert, Ryder. Approximate Bayesian computational methods

### Algorithm 1 (only for discrete y)
set.seed(1)
n = 15  # VERY slow for basically n > 15 ish (unless the prior is is chosen close to posterior)
y = sort(rpois(n, 2))
plot(table(y))

# prior: theta ~ Gamma(a, b)
a = 8
b = a/2
# sampler is faster if prior is closer to the posterior

N = 400
theta = double(N) # posterior draws

for (i in 1:N){
    cat("\r", i, "/", N)
    z = double(n)
    while(!all(z == y)){
        thetaprime = rgamma(1, a, b)
        z = sort(rpois(n, thetaprime))
        }
    theta[i] = thetaprime
    }

mean(theta)
var(theta)

plot(density(theta))
curve(dgamma(x, a + sum(y), b + n), add = TRUE, col = 'red') # true posterior
curve(dgamma(x, a, b), add = TRUE, col = 'blue') # prior


### Algorithm 2 (from Pritchard et al 1999) (generalized algorithm 1)
set.seed(1)
n = 15

y = rnorm(n, 0, sqrt(1))
plot(density(y))
z = double(n)

eta = function(x) mean(x)
rho = function(etaz, etay) abs(etaz - etay)^1
epsilon = 1e-1

# unknown mean, known variance
a = 0
b = 1

N = 10000
theta = double(N) # posterior draws

for (i in 1:N){
    cat("\r", i, "/", N)
    flag = TRUE
    while(rho(eta(z), eta(y)) > epsilon || flag){
        flag = FALSE
        thetaprime = rnorm(1, a, sqrt(b))
        z = rnorm(n, thetaprime, sqrt(1))
        }
    theta[i] = thetaprime
    }

mean(theta)
var(theta)

plot(density(theta))
curve(dnorm(x, (a/b + sum(y)/b) / (1/b + n), 1/sqrt(1/b + n)), col = 'red', add = TRUE)
curve(dnorm(x, a, sqrt(b)), col = 'blue', add = TRUE)

### Algorithm 3 (MCMC-ABC)
# Using same data used for algorithm 2
set.seed(1)
n = 15

y = rnorm(n, 0, sqrt(1))
plot(density(y))

eta = function(x) c(mean(x), sd(x))
rho = function(etaz, etay) sum(abs(etaz - etay)^2)
epsilon = 1e-2

# unknown mean, known variance
a = 0
b = 1

N = 50000
theta = double(N) # posterior draws
cand.sig = 1/4
accept = double(N)

# Get first draw using algorithm 2
z = matrix(0, N, n)
z[1,] = 10000
while(rho(eta(z[1,]), eta(y)) > epsilon){
    theta[1] = rnorm(1, a, sqrt(b))
    z[1,] = rnorm(n, theta[1], sqrt(1))
    }

# MCMC part
for (i in 2:N){
#   cat("\r", i, "/", N)
    cand.theta = rnorm(1, theta[i-1], cand.sig)
    cand.z = rnorm(n, cand.theta, sqrt(1))
    if ((log(runif(1)) <=
        dnorm(cand.theta, a, sqrt(b), log = TRUE) +
        dnorm(theta[i-1], cand.theta, cand.sig, log = TRUE) -
        dnorm(theta[i-1], a, sqrt(b), log = TRUE) - 
        dnorm(cand.theta, theta[i-1], cand.sig, log = TRUE))
        && (rho(eta(cand.z), eta(y)) <= epsilon)){

        theta[i] = cand.theta
        z[i,] = cand.z
        accept[i] = 1
    } else {
        theta[i] = theta[i-1]
        z[i,] = z[i-1,]
        }
    }

mean(accept)

plot(theta, type='l')

plot(density(theta))
curve(dnorm(x, (a/b + sum(y)/b) / (1/b + n), 1/sqrt(1/b + n)), col = 'red', add = TRUE)
curve(dnorm(x, a, sqrt(b)), col = 'blue', add = TRUE)

