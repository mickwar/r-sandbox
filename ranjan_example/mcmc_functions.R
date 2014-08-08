row.col = function(x, rows){
    out = double(ncol(x))
    for (i in 1:ncol(x))
        out[i] = x[rows[i], i]
    return (out)
    }
covmat = function(i = NA){
    LAMW = row.col(lam.w, apply(lam.w, 2, function(x) max(which(!is.na(x)))))
    RHOW = row.col(rho.w, apply(rho.w, 2, function(x) max(which(!is.na(x)))))
    LAMV = row.col(lam.v, apply(lam.v, 2, function(x) max(which(!is.na(x)))))
    RHOV = row.col(rho.v, apply(rho.v, 2, function(x) max(which(!is.na(x)))))
    LAMETA = row.col(lam.eta, apply(lam.eta, 2, function(x) max(which(!is.na(x)))))
    LAMY = row.col(lam.y, apply(lam.y, 2, function(x) max(which(!is.na(x)))))
    THETA = row.col(theta, apply(theta, 2, function(x) max(which(!is.na(x)))))
    if (!is.na(i)){
        LAMW = lam.w[i,]
        RHOW = rho.w[i,]
        LAMV = lam.v[i,]
        RHOV = rho.v[i,]
        LAMETA = lam.eta[i,]
        LAMY = lam.y[i,]
        THETA = theta[i,]
        }
    out = matrix(0, n*(p.delta+p.eta)+m*p.eta, n*(p.delta+p.eta)+m*p.eta)
    # would need to change if Fgroups>1 or p.eta>1, check notes in 2.2.4 higdon (2008)
    sig.v = kronecker(diag(1/LAMV, p.delta), R.x(field.X, RHOV))
    sig.u = 1/LAMW*R.x(cbind(field.X, THETA), RHOW)
    sig.uw = 1/LAMW*R.cross(cbind(field.X, THETA), sim.X, RHOW)
    sig.w = 1/LAMW*R.x(sim.X, RHOW)
    out = matrix.assemble(c(3,3), sig.v, NA, NA, NA, sig.u, sig.uw, NA, t(sig.uw), sig.w)
    out = out + diag2(solve(LAMY*t(B)%*%W%*%B + diag(ridge,n*(p.delta+p.eta))), solve(LAMETA*t(K)%*%K))
    return (out)
    }

#joint.out = rbind(as.matrix(field.yStd), t(Ksi)) # for regular
# using result 1

log.Like = function(G){
    CH = chol(G)
    -sum(log(diag(CH)))-1/2*t(z.vec) %*% chol2inv(CH) %*% z.vec
#   -1/2*determinant(G, logarithm=TRUE)$modulus[1]-1/2*t(z.vec)%*%solve(G)%*%z.vec
#    -1/2*determinant(G, logarithm=TRUE)$modulus[1]-1/2*t(joint.out)%*%solve(G)%*%joint.out
    }
log.Beta = function(x, params){
    a = params[1]
    b = params[2]
    (a-1)*log(x) + (b-1)*log(1-x)
    }
log.Gamma = function(x, params){
    a = params[1]
    b = params[2]
    (a-1)*log(x) - b*x
    }

compute.posterior = function(G){
    LAMW = row.col(lam.w, apply(lam.w, 2, function(x) max(which(!is.na(x)))))
    RHOW = row.col(rho.w, apply(rho.w, 2, function(x) max(which(!is.na(x)))))
    LAMV = row.col(lam.v, apply(lam.v, 2, function(x) max(which(!is.na(x)))))
    RHOV = row.col(rho.v, apply(rho.v, 2, function(x) max(which(!is.na(x)))))
    LAMETA = row.col(lam.eta, apply(lam.eta, 2, function(x) max(which(!is.na(x)))))
    LAMY = row.col(lam.y, apply(lam.y, 2, function(x) max(which(!is.na(x)))))
    THETA = row.col(theta, apply(theta, 2, function(x) max(which(!is.na(x)))))

    out = log.Like(G)
    # lam.w
    for (i in 1:p.eta)
        out = out + log.Gamma(LAMW[i], params[["lamw"]])

    # rho.w
    for (i in 1:(px + pt))
        out = out + log.Beta(RHOW[i], params[["rhow"]])

    # lam.v
    for (i in 1:Fgroups)
        out = out + log.Gamma(LAMV[i], params[["lamv"]])
   
    # rho.v
    for (i in 1:px)
        out = out + log.Beta(RHOV[i], params[["rhov"]])

    # lam.eta
    out = out + log.Gamma(LAMETA, params[["lameta"]])

    # lam.y
    out = out + log.Gamma(LAMY, params[["lamy"]])

    return (out)
    }

sig.calc = function(x, a, b, c)
    (sqrt(2*pi)*((b-a)*exp(-c*x)+a))^(-1)
