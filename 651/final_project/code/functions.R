# logit function
logit = function(x)
    log(x/(1-x))
# covariance for simulator input (eq. 1)
cov.simul = function(params, x, t){
    lam.eta = params[1]
    rho.x = params[2:(1+px)]
    rho.t = params[(2+px):(2+pt)]


    }
