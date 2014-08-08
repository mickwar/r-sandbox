# Not all computer experiments are costly, some run rather
# quickly. If experiments were run without a proper design,
# there may be design points too close together, which could
# give very similar results and hence are unnecessary.

# It may be desirable to find a subset of the design that
# will still give good analytic properties. This function
# finds a subset of the original design that is similar to a
# maximin design.

# The algorithm is adapted from Morris and Mitchell (1995)
# "Exploratory designs for computational experiments" and
# Joseph and Hung (2008) "Orthogonal-maximin Latin
# hypercube designs"

# note: proposition: within any subset of a design matrix, an upper
#       bound for phi is that phi calculated using the distance
#       matrix of the complete design matrix.
#       (this seems to be the case only for high p (I tested it on
#       p = 15, but p = 1 and p = 2, yielded various ranges))
# i think i can get upper and lower bounds by simulation, the design
# is given and so is the number of observations to be used in the
# subset. start with the overall distance matrix, then take random
# subsets and repeatedly calculate phi (with func.phi())

near.maximin = function(design, m, p = 2, alp1 = 1, alp2 = 1,
    w = 0.5, niter = 10000){

    x = matrix(c(rbeta(50, 1, 3), rbeta(50, 4, 1),
        rbeta(50, 4, 1), rbeta(50, 1, 3)), 100, 2)
    x = rbind(x, matrix(runif(20), 10, 2))
    plot(x, pch=20, xlim=c(0,1), ylim=c(0,1))

    # design: the current design
    # m: the number of design points to use in the subset
    x = as.matrix(design)
    n = nrow(x)
    d = ncol(x)

    # functions
    dbar = (m+1)*d/3
    phi.L = (choose(m, 2)*((ceiling(dbar)-dbar)/(floor(dbar)^p)+
        (dbar-floor(dbar))/(ceiling(dbar)^p)))^(1/p)
    phi.U = (sum(((m-1):1)/(1:(m-1)*d)^p))^(1/p)
    ceiling = function(x)
        trunc(x)+1
    # calculate upper and lower bounds for phi to scale phi to [0,1]
    func.rho.col = function(COR)
        (apply(COR^2, 1, sum)-1)/(d-1)
    func.phi.row = function(DIST)
        apply(1/(DIST^p+diag(Inf, m)), 1, sum)^(1/p)
    func.rho = function(COR)
        (sum(COR^2)-d)/(d*(d-1))
    func.phi = function(DIST)
        (sum(func.phi.row(DIST)^p)/2)^(1/p) / (m-1)
    func.psi = function(COR, DIST){
        rho = func.rho(COR)
        phi = func.phi(DIST)
        w*rho + (1-w)*(phi-phi.L)/(phi.U-phi.L)
        }

    all.dist = as.matrix(dist(x, method="manhattan"))
    all.corr = cor(x)

    # initial sample
    test.index = sample(n, m)
    test.x = x[test.index, ]
    test.dist = all.dist[test.index, test.index]
    test.corr = cor(x[test.index,])

    func.rho(test.corr)
    (func.phi(test.dist) - phi.L)/(phi.U - phi.L)
    func.psi(test.corr, test.dist)

    points(test.x, col='red')
    
    # weights for new samples
    h = double(n)+1

    count = 0
    while (count < niter){
        count = count + 1
        z = calc.phi(all.dist[test.index, test.index])
        move = sample(m, 1, prob=(z^alp1)/sum(z^alp1))
        h[test.index[move]] = h[test.index[move]] + 1
        h[test.index[-move]] = h[test.index[-move]] - 1/(2*m)
        h = ifelse(h < 1, 1, h)
        g = 1/h[-test.index]
        new = sample((1:n)[-test.index], 1, prob=(g^alp2)/sum(g^alp2))
        test.x[move,] = x[new,]
        test.index[move] = new
#       plot(x)
#       points(test.x, pch=20, col='red', lwd=8)
#       Sys.sleep(0.01)
        }
    out = x[order(h)[1:m],]
    plot(out)

    if (print.h)
        print(h)
    return(out)
    }

n = 200
k = 2
x = matrix(runif(n*k), n, k)
m = 10

library(randtoolbox)
x = sobol(n, k)

x = matrix(c(0.1, 0.9, 0.1, 0.8, 0.2, 0.9, 0.2, 0.8,
    0.15, 0.3, 0.7, 0.85, 0.65, 0.35, 0.8, 0.2), 8, 2, byrow=T)
plot(x)
points(near.maximin(x, nrow(x)/2, alp1=10), col='red', cex=2)
points(near.maximin(x, nrow(x)/2, alp2=10), col='blue', cex=4)
points(near.maximin(x, nrow(x)/2, p=10), col='green', cex=6)
points(near.maximin(x, nrow(x)/2, p=10, alp1=10, alp2=10), col='orange', cex=8)

out1 = near.maximin(x, m, p=1, niter=10000, alp1=1, alp2=1)
out2 = near.maximin(x, m, p=1, niter=10000, alp1=1, alp2=20)
out3 = near.maximin(x, m, p=1, niter=10000, alp1=20, alp2=1)
out4 = near.maximin(x, m, p=1, niter=10000, alp1=50, alp2=50)

plot(x, xlim=c(0,1), ylim=c(0,1))
points(jitter(out1, 1), col='red', lwd=2)
points(jitter(out2, 1), col='green', lwd=2)
points(jitter(out3, 1), col='blue', lwd=2)
points(jitter(out4, 1), col='yellow', lwd=2)



### AFRL EXAMPLE
dat = read.csv("~/files/afrl_project/dat/triflate_samples.csv",
    sep=",", header=FALSE)

x = as.vector(as.matrix(dat[1, 23:33]))
ord = order(x)
x = x[ord]

design = matrix(as.numeric(as.matrix(dat[-1, c(1:20)[-18]])), 9930, 19)

# re-scale
for (i in 1:19){
    design[,i] = (design[,i] - min(design[,i])) /
        (max(design[,i]) - min(design[,i]))
    }
# remove duplicate rows
design = unique(design)


n = dim(design)[1]
m = ceiling(n*0.1)

times = matrix(0, 100, 3)
for (i in 1:100)
    times[i,] = system.time(near.maximin(design, m, p=1, niter=i, alp1=20, alp2=20))[1:3]
plot(times[,1])
points(times[,3],col='red')

points(times[,2],col='blue')

y = times[,1]
x = 1:100
mod = lm(y ~ x)
### END AFRL EXAMPLE




make.lh = function(n, d){
    out = matrix(0, n, d)
    for (i in 1:d)
        out[,i] = sample(n)
    return (out)
    }
ceiling = function(x)
    trunc(x)+1
# calculate upper and lower bounds for phi to scale phi to [0,1]
func.rho.col = function(COR)
    (apply(COR^2, 1, sum)-1)/(d-1)
func.phi.row = function(DIST)
    apply(1/(DIST^p+diag(Inf, m)), 1, sum)^(1/p)
func.rho = function(COR)
    (sum(COR^2)-d)/(d*(d-1))
func.phi = function(DIST)
    (sum(func.phi.row(DIST)^p)/2)^(1/p)
func.psi = function(COR, DIST){
    rho = func.rho(COR)
    phi = func.phi(DIST)
    w*rho + (1-w)*(phi-phi.L)/(phi.U-phi.L)
    }
func.psi2 = function(COR, DIST){
    rho = func.rho(COR)
    phi = func.phi(DIST) / (n-1)
    w*rho + (1-w)*(phi-phi.L)/(phi.U-phi.L)
    }

n = 25
d = 2
p = 2.1
ra = matrix(runif(n*d), n, d)
ra.dist = as.matrix(dist(ra, method="manhattan"))

m = 23
niter = 1000
ra.phi = double(niter)
for (i in 1:niter){
    ind = sample(n, m)
    ra.phi[i] = func.phi(ra.dist[ind, ind])
    } 

plot(ra.phi, pch=20); range(ra.phi)
