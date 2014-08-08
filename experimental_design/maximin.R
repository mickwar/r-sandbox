### Distance: Maximin (maximizing the minimum distance, not from a latin hypercube)
maximin=function(default=0,n,dimension=2,domain=c(0,1),
    maxdif=sqrt(n)/100, showplot=TRUE, sobol=TRUE){
    d=dimension
    out=default
    if (!is.matrix(default)){
        # initial scatter
        if (sobol){
            require(randtoolbox)
            out=sobol(n,d)
            }
        if (!sobol)
            out=matrix(runif(n*d),n,d)
        if (d==1)
            out=matrix(out,n,1)
        }
    if (showplot){
        if (d==1)
            plot(out,integer(n),xlim=c(-0.05,1.05),ylim=c(-0.5,0.5))
        if (d==2)
            plot(out,xlim=c(-0.05,1.05),ylim=c(-0.05,1.05))
        if (d > 2)
            print ("Plots not available for dimensions higher than 2.")
        }
    #while ( !all(abs(out - all.mu)<=all.sig) ){
    niter = 0
    while (niter < 10000){
        niter = niter+1
        H = as.matrix(dist(out)) + diag(Inf,n)
        # get pair index for smallest distances
        pair.index = arrayInd(which.min(H), dim(H))
        pair.dist=min(H)
        pair.index=sample(pair.index)
        direction=out[pair.index[1],]-out[pair.index[2],]
        scalar=runif(d,0,maxdif)
        # move one of the points away from the other closest point
        out[pair.index[1],]=out[pair.index[1],]+scalar*direction
        # make sure the points stay in [0,1]
        out[pair.index[1],]=ifelse(out[pair.index[1],]<0, 0, out[pair.index[1],])
        out[pair.index[1],]=ifelse(out[pair.index[1],]>1, 1, out[pair.index[1],])
        if (showplot){
            if (d==1){
                plot(out,integer(n),xlim=c(-0.05,1.05),ylim=c(-0.5,0.5))
                points(out[pair.index,],c(0,0),col='red',pch=20)
                }
            if (d==2){
                plot(out,xlim=c(-0.05,1.05),ylim=c(-0.05,1.05))
                points(out[pair.index,],col='red',pch=20)
                }
            }
        #print(abs(out - all.mu)<all.sig)
        #print(all.sig)
        }
    return(out)
    }

### Maximin from a latin hypercube
lh = function(n, dimensions){
    d = dimensions
    out = matrix(seq(0, 1, length=n), n, d)
    for (i in 1:d)
        out[,i] = sample(out[,i])
    return (out)
    }
exchange = function(latin){
    j = sample(ncol(latin), 1)
    i = sample(nrow(latin), 2)
    latin[i,j] = latin[i[2:1], j]
    return (latin)
    }
maximin = function(n, dimensions, niter){
    d = dimensions
    out = lh(n, d)
    test.lh = out
    min.dist = min(dist(out))
    print(min.dist)
    count = 0
    plot(out, pch=20)
    while (count < niter){
        count = count + 1
        #test.lh = lh(n, d)  # permutation selection
        test.lh = exchange(test.lh)     # exchange
        test.dist = min(dist(test.lh))
        if (test.dist > min.dist){
            min.dist = test.dist
            print(min.dist)
            out = test.lh
            plot(out, pch=20)
            }
        }
    return (out)
    }
out = maximin(16, 2, 100)
out2 = maximin(16, 2)
plot(out2, pch=20); points(out, col='green', pch=20)

### morris and mitchell (1995) algorithm (at least what I can
### glean from joseph and hung (2008))
### should include some "best" design
lh = function(n, dimensions){
    d = dimensions
    out = matrix(seq(0, 1, length=n), n, d)
    for (i in 1:d)
        out[,i] = sample(out[,i])
    return (out)
    }
exchange = function(latin){
    j = sample(ncol(latin), 1)
    i = sample(nrow(latin), 2)
    latin[i,j] = latin[i[2:1], j]
    return (latin)
    }
maximin = function(n, dimensions, t=1, niter=10000){
    func.phi = function(x, p=2)
        (sum(1/dist(x)^p))^(1/p)
    d = dimensions
    out = lh(n, d)
    test.lh = out
    phi.out = func.phi(out, 1)
    print(phi.out)
    count = 0
    plot(out, pch=20)
    while (count < niter){
        count = count + 1
        test.lh = exchange(test.lh)
        phi.test = func.phi(test.lh, 1)
        prob = min(1, exp(-(phi.test-phi.out)/t))
        if (prob > runif(1)){
            phi.out = phi.test
            print(phi.out)
            out = test.lh
            plot(out, pch=20)
            }
        }
    return (out)
    }
out = maximin(25, 2, t=0.1, niter=100000)


### joseph and hung's (2008) algorithm
maximin = function(n, k, w=0.5, t=1, p=2, alpha=1, niter=1000,
    show.psi=FALSE, show.plot=FALSE){
    # n - number of data points
    # k - number of dimensions
    # w - weight parameter 0<w<1 (see func.psi)
    # t - "temperature" parameter t>0 (see calculation of 'prob' in
    #     the while loop, larger t accepts proposal designs more often
    # p - positive integer, 1 or 2 should be good
    # alpha - alpha >= 1, used in choosing row and column in the
    #     exchange, larger alpha inflates probability of choosing
    #     the row/column with highest correlation/smallest distance
    
    # function for random latin hypercube
    lh = function(n, dimensions){
        d = dimensions
        out = matrix(1:n, n, d)
        for (i in 1:d)
            out[,i] = sample(out[,i])
        return (out)
        }
    # FIX. this is what's making psi not between 0 and 1
    dbar = (n+1)*k/3 # only applies when factors take on values 1,..,n
    ceiling = function(x)
        trunc(x)+1
    # calculate upper and lower bounds for phi to scale phi to [0,1]
    phi.L = (choose(n, 2)*((ceiling(dbar)-dbar)/(floor(dbar)^p)+
        (dbar-floor(dbar))/(ceiling(dbar)^p)))^(1/p)
    phi.U = (sum(((n-1):1)/(1:(n-1)*k)^p))^(1/p)
    # set up functions to calculate various components in the algorithm
    func.rho.col = function(COR)
        (apply(COR^2, 1, sum)-1)/(k-1)
    func.phi.row = function(DIST)
        apply(1/(DIST^p+diag(Inf, n)), 1, sum)^(1/p)
    func.rho = function(COR)
        (sum(COR^2)-k)/(k*(k-1))
    func.phi = function(DIST)
        (sum(func.phi.row(DIST)^p)/2)^(1/p)
    func.psi = function(COR, DIST){
        rho = func.rho(COR)
        phi = func.phi(DIST)
        w*rho + (1-w)*(phi-phi.L)/(phi.U-phi.L)
        }

    # begin with random latin hypercube
    out.lh = lh(n, k)
    out.corr = cor(out.lh)
    out.dist = as.matrix(dist(out.lh))
    out.psi = func.psi(out.corr, out.dist)
    if (show.psi)
        print(out.psi)

    # initialize variables for proposal design
    test.lh = out.lh
    test.corr = out.corr
    test.dist = out.dist

    if (show.plot && k == 2)
        plot(out.lh, pch=20)

    count = 0
    while (count < niter){
        count = count + 1

        y = func.rho.col(test.corr)
        z = func.phi.row(test.dist)

        # choosing column
        if (all(y==0)){ # in the case of no correlation
            l.star = sample(k, 1)
        } else {
            l.star = sample(k, 1, prob=(y^alpha)/sum(y^alpha))
            }
        # choosing row
        if (all(z==0)){ # probably won't occur
            i.star = sample(n, 1)
        } else {
            i.star = sample(n, 1, prob=(z^alpha)/sum(z^alpha))
            }
        # random row with which to swap
        i.prime = sample((1:n)[-i.star], 1)

        # make the swap
        test.lh[c(i.star, i.prime), l.star] = test.lh[c(i.prime, i.star), l.star]

        # calculate anew correlation and distance variables
        test.corr = cor(test.lh)
        test.dist = as.matrix(dist(test.lh))

        test.psi = func.psi(test.corr, test.dist)

        prob = min(1, exp(-(test.psi-out.psi)/t))

        if (prob > runif(1)){
            # update out.lh and other variables
            out.lh = test.lh
            out.corr = test.corr
            out.dist = test.dist
            out.psi = test.psi
            if (show.psi)
                print(out.psi)
            if (show.plot && k == 2)
                plot(out.lh, pch=20)
            }
        }
    return (out.lh)
    }
out = maximin(50, 2, w=0.5, t=0.1, p=15, alpha=10, niter=1000,
    show.psi=TRUE, show.plot=TRUE)



##################################
maxi.subset = function(design, lh, ask=FALSE, bonus=FALSE){
    design = as.matrix(design)
    lh = as.matrix(lh)
    n = nrow(design)
    k = ncol(design)
    if (k != ncol(lh))
        stop("Latin hypercube must have same dimensions as design.")
    m = nrow(lh)
#    radius = min(dist(lh))/2
    dist2 = function(x, y){ # the lower left matrix of
                            # dist(rbind(x, y))
        x = as.matrix(x) # m x k
        y = as.matrix(y) # n x k
        out = matrix(0, nrow(y), nrow(x)) # n x m
        for (j in 1:nrow(x)){
            for (i in 1:nrow(y)){
                A = t(y) - x[j,]
                out[i,j] = sqrt(t(A[,i])%*%A[,i])
                }
            }
        return (out)
        }
    if (bonus && k == 2){
        plot(design, xlim=c(0,1), ylim=c(0,1))
        points(lh, col='red',pch=20)
        }
#    circle = function(h,k,r){    # for 2-d
#        out = matrix(0, 200, 2)
#        out[1:100,1] = seq(-r, r, length=100)
#        out[101:200,1] = seq(r, -r, length=100)
#        out[,2] = sqrt(r^2 - out[,1]^2)
#        out[1:100,2] = -out[1:100,2]
#        out[,1] = out[,1]+h
#        out[,2] = out[,2]+k
#        return(out)
#        }
    new = matrix(NA, m, k)
    dists = rbind(as.matrix(dist(lh))+diag(Inf, m), dist2(lh, design))
    total.index = NULL
    for (i in 1:m){
        if (ask)
            readline("Hit ENTER for next iteration:")
#        dists = dist(lh[i,], design)
#        index = which(dists<radius)
#        if (length(index) > 0)
#            index = which.min(dists) # get closest point regardless
        index = which.min(dists[,i])
        dists[i,] = Inf
        if (index > m){      # closet point is in the design set
            dists[index,] = Inf # remove design point from being chosen twice
            new[i,] = design[index-m,]
            total.index = c(total.index, index-m)
            }
#        total.dists = c(total.dists, min(dists, na.rm=TRUE))
#        design[index,] = NA
#        if (bonus && k == 2){
#            white = design[index,]
#            points(white, col='white',lwd=2)
#            points(lh, col='red',pch=20)
#            points(new, col='green',pch=20)
#            circ = circle(lh[i,1], lh[i,2], radius)
#            points(circ, col='green',type='l')
#            }
        }
    return(list("points"=new[!is.na(new[,1]),], "coverage"=length(total.index)/m,
        "index"=total.index))#, "dists"=total.dists))
    }

### AFRL EXAMPLE
dat = read.csv("~/files/afrl_project/dat/triflate_samples.csv",
    sep=",", header=FALSE)

x = as.vector(as.matrix(dat[1, 23:33]))
ord = order(x)
x = x[ord]

inputs = matrix(as.numeric(as.matrix(cbind(dat[-1, 1:20], 0))), 9930, 21)

pairs(inputs[,11:13])

dat = dat[-1, c((23:33)[ord], (35:45)[ord], (47:57)[ord], (59:69)[ord],
	(71:81)[ord])]

colors = c("green", "red", "blue", "orange", "black")
labelz = c("CF3SO3-", "CF3SO2-", "F-", "CF3-", "SO3-")

m = dim(dat)[1]
px = length(x)
pt = dim(inputs)[2]

design = matrix(rep(inputs, px, each=px), m*px, pt)
design[,pt] = rep(x, times=m)
test.set = maximin(1000, pt, w=0.5, t=0.1, alpha=10,
    niter=10, show.psi=TRUE)
test.set = (test.set-1)/999 # scale back to [0,1]
out = maxi.subset(design, test.set, ask=FALSE, bonus=FALSE)
### END AFRL EXAMPLE


design.set = matrix(rbeta(n*k, 1.3, 0.6), n, k)
#design.set = matrix(runif(n*k), n, k)
test.set = maximin(m,k,w=0.5,t=0.1,alpha=1,niter=1000)
#test.set = seq(0, 1, length=m)
#pairs(test.set)
out = maxi.subset(design.set, test.set, ask=FALSE, bonus=FALSE)
plot(design.set, xlim=c(0,1), ylim=c(0,1))
points(test.set, col='red')
points(out$points, col='green', pch=20)


out = maxi.subset(design.set, test.set, ask=TRUE, bonus=TRUE)
# getting low coverage with high k, but is that a problem?
# not working with high dimensions, maybe use a rectangles
# instead of balls?

pairs(rbind(test.set, out$points), col=rep(c("red","green"),
    each=nrow(test.set)),pch=rep(c(1,20),each=nrow(test.set)),
    xlim=c(0,1), ylim=c(0,1))

n = 100
m = 10
nreps = 10
maxk = 5
coverage = matrix(0, nreps, maxk)
for (k in 1:maxk){
    for (i in 1:nreps){
        design.set = matrix(runif(n*k), n, k)
        if (k > 1){
            test.set = maximin(m,k,w=0.5,t=0.1,alpha=1,niter=1000)
        } else {
            test.set = seq(0, 1, length=m)
            }
        out = maxi.subset(design.set, test.set, ask=FALSE, bonus=FALSE)
        coverage[i, k] = out$coverage
        }
    }
apply(coverage, 2, mean)
apply(coverage, 2, sd)
apply(coverage, 2, range)

plot(design.set)
points(test.set, col='red')
points(out$points, col='green', pch=20)


pairs(test.set)

out = maxi.subset(design.set, test.set, ask=FALSE, bonus=FALSE)



# library(rgl)
# plot3d(design.set)
# points3d(out$points, col='green',size=5)
# points3d(test.set,col='red')

pairs(design.set)
pairs(out$points, col='red')

# metrics
nreps = 10
loops = 100
n = floor(seq(100, 1000, length=nreps))
sums = double(nreps)
for (i in 1:nreps){
    temp = double(loops)
    for (j in 1:loops){
        design.set = matrix(runif(n[i]*k), n[i], k)
        out = maxi.subset(design.set, test.set, ask=FALSE, bonus=FALSE)
        temp[j] = sum(out$dists)
        }
    sums[i] = mean(temp)
    }
plot(n, sums, type='l')


# do with real data now
