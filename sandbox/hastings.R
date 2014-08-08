# section 3 of hastings (1970) describes a matrix E(theta)

# m is positive integer
# 1 <= i != j <= m integers
# theta is in [0, 2*pi]

set.seed(1)
E = function(m, theta, i, j){
    out = diag(m)
    out[i, i] = cos(theta)
    out[i, j] = sin(theta)
    out[j, i] = -sin(theta)
    out[j, j] = cos(theta)
    return (out)
    }

make.H = function(m){
    out = matrix(0, m, m)
    for (j in 1:m){
        out[1,j] = 1/sqrt(m)
        out[2,j] = 1/sqrt(m)*cos((j-1)/pi)
        for (r in 3:(0.5*m+1))
            out[r,j] = sqrt(2/m)*cos((j-1)*(r-2)*2*pi/m)
        for (s in (0.5*m+2):m)
            out[s,j] = sqrt(2/m)*sin((j-1)*(s-0.5*m-1)*2*pi/m)
        }
    return (out)
    }

# starting point
m = 50
#H0 = svd(matrix(rnorm(6^2), 6, 6))$v
H0 = make.H(m)


# example 4
L = 25
niter = 1000
big.J = double(L)
for (l in 1:L){
    H = rep(list(matrix(0, m, m)), niter)
    i = sample(1:m, 1)
    j = sample((1:m)[-i], 1)
    theta = runif(1, 0, 2*pi)
    H[[1]] = E(m, theta, i, j) %*% H0
    for (t in 2:niter){
        i = sample(1:m, 1)
        j = sample((1:m)[-i], 1)
        theta = runif(1, 0, 2*pi)
        H[[t]] = E(m, theta, i, j) %*% H[[t-1]]
        }

    f = function(H)
        sum(diag(H)^2)
    J = double(niter)
    for (t in 1:niter)
        J[t] = f(H[[t]])
    big.J[l] = mean(J)
    # J = matrix(0, m, m)
    # for (t in 1:niter)
    #     J = J + f(H[[t]])
    # J = J / niter
}

mean(big.J)
sd(big.J)
plot(density(big.J))

# for (i in 1:5000){
#     plot(density(J[i:niter]))
#     plot(J[i:niter], type='l')
#     Sys.sleep(0.1)
#     }
