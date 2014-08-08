# need to add some "best" design instead of always switching
# with that probability
# look at morris and mitchel for more parts on the algorithm

### joseph and hung's (2008) algorithm
maximin = function(n, k, w=0.5, t=1, p=2, alpha=1, niter=1000,
    show.psi=FALSE, show.plot=FALSE, pb=FALSE){
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
    if (pb)
        source("~/files/R/pb_linux.R")
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
        if (pb)
            pb.linux(count, niter)

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
