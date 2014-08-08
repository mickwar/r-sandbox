### joseph and hung (2008) algorithm, a modification
### to the simulated annealing algorithm proposed
### by morris and mitchell (1995), both papers
### should be referred to in coding the algorithm
omlhd = function(n, k, w=0.5, t=1, fac.t = 0.95, p=15, alpha=5, imax=100){
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
    # lemma 1
    dbar = (n+1)*k/3

    # ceiling function as defined by joseph and hung (2008) page 176
    # (doesn't work for negative numbers, but none are used in the
    # algorithm)
    ceiling = function(x)
        trunc(x)+1

    # calculate upper and lower bounds for phi to scale phi to [0,1]
    phi.L = (choose(n, 2)*((ceiling(dbar)-dbar)/(floor(dbar)^p)+
        (dbar-floor(dbar))/(ceiling(dbar)^p)))^(1/p)
    phi.U = (sum(((n-1):1)/(1:(n-1)*k)^p))^(1/p)

    # set up functions to calculate various components in the algorithm
    # equation (1)
    func.rho = function(COR)
        (sum(COR^2)-k)/(k*(k-1))
    # phi_p on page 176
    func.phi = function(DIST)
        (sum(func.phi.row(DIST)^p)/2)^(1/p)
    # psi_p on page 176
    func.psi = function(COR, DIST){
        rho = func.rho(COR)
        phi = func.phi(DIST)
        w*rho + (1-w)*(phi-phi.L)/(phi.U-phi.L)
        }
    # equation (3)
    func.rho.col = function(COR)
        (apply(COR^2, 2, sum)-1)/(k-1)
    # equation (4)
    func.phi.row = function(DIST)
        apply(1/(DIST^p+diag(Inf, n)), 1, sum)^(1/p)

    # begin with random latin hypercube
    out.lh = lh(n, k)
    out.corr = cor(out.lh)
    out.dist = as.matrix(dist(out.lh, method="manhattan"))
    out.psi = func.psi(out.corr, out.dist)

    # initialize variables for proposal design
    try.lh = out.lh
    try.corr = out.corr
    try.dist = out.dist

    # initialize best design
    best.lh = out.lh
    best.psi = out.psi

    tflag = TRUE
    at = 1
    while (tflag){
        tflag = FALSE
        i = 1
        print(c(best.psi, t, at))
        while (i < imax){
            y = func.rho.col(try.corr)
            z = func.phi.row(try.dist)

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
            try.lh[c(i.star, i.prime), l.star] = try.lh[c(i.prime, i.star), l.star]

            # calculate new correlation and distance variables
            try.corr = cor(try.lh)
            try.dist = as.matrix(dist(try.lh, method="manhattan"))

            try.psi = func.psi(try.corr, try.dist)

            if (log(runif(1)) < -(try.psi - out.psi)/t){
                # update out.lh and other variables
                out.lh = try.lh
                out.corr = try.corr
                out.dist = try.dist
                out.psi = try.psi
                if (!tflag)
                    at = i 
                tflag = TRUE
                }
            if (try.psi < best.psi){
                best.psi = try.psi
                best.lh = try.lh
                i = 1
            } else {
                i = i + 1
                }
            }
        if (tflag)
            t = t * fac.t
        }
    return (best.lh)
    }
