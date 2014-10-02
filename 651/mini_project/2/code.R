dat = as.matrix(read.table("./faculty.dat"))
colnames(dat) <- rownames(dat) <- NULL

y = as.vector(dat)

# Likelihood: Normal, unknown mean and variance
# Prior: mu ~ Normal, sigma^2 ~ I.G, indepdent of each other

# inverse gamma density
igpdf=function(x, a, b)
	(b^a)/gamma(a)*x^(-a-1)*exp(-b/x)

cross = function(...){
    vec = list(...)
    d = length(vec)
    N = double(d+2) + 1
    for (i in 1:d)
        N[i+1] = length(vec[[i]])
    out = matrix(0, prod(N), d)
    for (i in 1:d){
        out[,i] = rep(vec[[i]], times=prod(N[1:i]),
            each=prod(N[(i+2):(d+2)]))
        }
    return(out)
    }

# unnormalize posterior
# m, s2, a, b are hyperparameters
g = function(x, m = 4, s2 = 100, a = 2.5, b = 1.5){
    mu = x[1]
    sig2 = x[2]
    prod(dnorm(y, mu, sqrt(sig2))) * dnorm(mu, m, sqrt(s2)) *
        igpdf(sig2, a, b)
    }

# envelope function
#dW = function(x){
#    mu = x[1]
#    sig2 = x[2]
#    dnorm(mu, mean(y), sqrt(0.04)) * dnorm(sig2, 0.3, sqrt(0.04))
#    }
#rW = function(n)
#    matrix(c(rnorm(n,mean(y),sqrt(0.04)),rgamma(n,5.5,scale=1/16)),n,2)
dW = function(x){
    mu = x[1]
    sig2 = x[2]
    dcauchy(mu-mean(y), scale=1/7) * igpdf(sig2, 5, 0.28*(5+1))
    }
rW = function(n)
    matrix(c(rcauchy(n, scale=1/7)+mean(y),
        1/rgamma(n, 5, rate=0.28*(5+1))),n,2)

mu.vec = seq(4.0, 7.0, length=500)
sig.vec = seq( 0.01, 2.0, length=500)
xy = cross(mu.vec, sig.vec)

z = double(nrow(xy))
w = double(nrow(xy))
for (i in 1:length(z)){
    z[i] = g(xy[i,])
    w[i] = dW(xy[i,])
    }

# estimate of the constant to bring g down
H = G
(G = max(z/w))
g.star = function(x)
    g(x) / G 

#h = double(nrow(xy))
#for (i in 1:length(z))
#    h[i] = g.star(xy[i,])

z.star = z / G
plot3d(cbind(xy, z.star))
points3d(cbind(xy, w), col='red')
plot3d(cbind(xy, w), col='red')
points3d(cbind(xy, z.star))
#points3d(cbind(xy, h), col='blue')

reject = function(M = 100){
    # initialize
    out = matrix(0, M, 5)

    # envelope draws
    out[,1:2] = rW(M)

    # uniform draws
    out[,3] = runif(M)
    
    # calculate ratios
    for (i in 1:M)
        out[i,4] = g.star(out[i,1:2]) / dW(out[i,1:2])

    # accept or not
    for (i in 1:M){
        out[i,5] = ifelse(out[i,3] <= out[i,4], 1, 0)
        }

    return (out)
    }

X = reject(3000000)
# porportion of acceptances (about 0.35 so far)
mean(X[,5])

# count of acceptances
sum(X[,5])

# ratios should not exceed 1, otherwise not a correct envelope
max(X[,4])

# accepted draws
Y = X[X[,5] == 1, 1:2]

# joint posterior
plot(Y[sample(nrow(Y), 10000),], pch=20)

# marginal posterior distributions
plot(density(Y[sample(nrow(Y), 100000),1]))
plot(density(Y[sample(nrow(Y), 100000),2]))

plot(density(y))
curve(dnorm(x, mean(Y[,1]), sqrt(mean(Y[,2]))), add=TRUE, col='blue')
