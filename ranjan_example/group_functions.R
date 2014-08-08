new.seq = function(x, back){
    out = seq(x-back+1, x)
    return (out[out > 0])
    }
last.na = function(x) # not really, more like last non na
    x[max(which(!is.na(x)))]
covmat = function(i = NA){ # computer sigma_z_hat
    LAMW = last.na(draws[["lamW"]])
    LAMV = last.na(draws[["lamV"]])
    LAMY = last.na(draws[["lamY"]])
    LAMETA = last.na(draws[["lamEta"]])
    RHOW1 = last.na(draws[["rhoW1"]])
    RHOW2 = last.na(draws[["rhoW2"]])
    RHOW3 = last.na(draws[["rhoW3"]])
    RHOV1 = last.na(draws[["rhoV1"]])
    RHOV2 = last.na(draws[["rhoV2"]])
    THETA = last.na(draws[["theta"]])
    if (!is.na(i)){
        LAMW = draws[["lamW"]][i]
        LAMV = draws[["lamV"]][i]
        LAMY = draws[["lamY"]][i]
        LAMETA = draws[["lamEta"]][i]
        RHOW1 = draws[["rhoW1"]][i]
        RHOW2 = draws[["rhoW2"]][i]
        RHOW3 = draws[["rhoW3"]][i]
        RHOV1 = draws[["rhoV1"]][i]
        RHOV2 = draws[["rhoV2"]][i]
        THETA = draws[["theta"]][i]
        }
    RHOW = c(RHOW1, RHOW2, RHOW3)
    RHOV = c(RHOV1, RHOV2)
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
    }
#log.Like = function(G)
#    -1/2*determinant(G, logarithm=TRUE)$modulus[1]-1/2*t(z.vec) %*% solve(G) %*% z.vec
#    -1/2*determinant(G, logarithm=TRUE)$modulus[1]-1/2*t(joint.out)%*%solve(G)%*%joint.out
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
    out = log.Like(G)
    out = out + log.Gamma(last.na(draws[["lamW"]]), params[["lamW"]])
    out = out + log.Gamma(last.na(draws[["lamV"]]), params[["lamV"]])
    out = out + log.Gamma(last.na(draws[["lamY"]]), params[["lamY"]])
    out = out + log.Gamma(last.na(draws[["lamEta"]]), params[["lamEta"]])
    out = out + log.Beta(last.na(draws[["rhoW1"]]), params[["rhoW1"]])
    out = out + log.Beta(last.na(draws[["rhoW2"]]), params[["rhoW2"]])
    out = out + log.Beta(last.na(draws[["rhoW3"]]), params[["rhoW3"]])
    out = out + log.Beta(last.na(draws[["rhoV1"]]), params[["rhoV1"]])
    out = out + log.Beta(last.na(draws[["rhoV2"]]), params[["rhoV2"]])
    return (out)
    }
