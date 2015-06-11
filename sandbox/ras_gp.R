x = c(-4.6, -4.2, -4.0, -3.8, -3.7,
      -3.1, -2.0, -1.9, -1.6, -0.5,
      -0.1,  0.1,  0.4,  0.5,  0.8,
       1.2,  1.5,  1.7,  2.5,  3.1)

y = c( 5.1,  4.8,  3.9,  3.5,  3.9,
       2.0,  1.3,  1.7,  1.9, -2.3,
      -3.5, -2.8, -2.4, -2.2, -0.4,
       1.0,  1.1,  1.3,  1.0,  2.1)

plot(x, y, pch=20, xlim=c(-5, 5), ylim=c(-4, 8))
abline(v=(-5:5), h=seq(-4, 8, by=2), col="gray50", lty=2)

m = function(x, theta_m){
    a = theta_m[1]
    b = theta_m[2]
    c = theta_m[3]
    return (a*x^2 + b*x + c)
    }
k = function(x, xprime, theta_k){
    # if xprime = x, leave xprime missing (i.e. a covariance matrix
    # is created for all possible pairs in x)
    n = NROW(x)
    
    sig_y = theta_k[1]
    sig_n = theta_k[2]
    l = theta_k[3]

    if (missing(xprime)){
        d = as.matrix(dist(x, method = "manhattan"))
        out = sig_y * exp(-(d^2) / (2*l^2)) + diag(sig_n, n)
    } else {
        d = as.matrix(dist(c(x, xprime), method = "manhattan"))
        d = d[1:n, (n+1):(n+m)]
        out = sig_y * exp(-(d^2) / (2*l^2)) 
        }
    return (out)
    }

# log posterior
calc.post = function(params){
    mu = m(x, theta_m = params[1:3])
    sigma = k(x, theta_k = params[4:6])
    # likelihood
    out = -0.5 * determinant(sigma)$modulus[1] +
        -0.5 * t(y-mu) %*% (solve(sigma) %*% (y-mu))
    # priors
    # none for now
    return (out)
    }


