# draws from the unit circle (Gibbs)

n = 10000
p = 2
r = 1
x = matrix(0, n, p)

for (i in 2:n){
    x[i, 1] = runif(1, -sqrt(r^2 - x[i-1,2]^2), sqrt(r^2 - x[i-1,2]^2))
    x[i, 2] = runif(1, -sqrt(r^2 - x[i,1]^2),   sqrt(r^2 - x[i,1]^2))
    }

# marginally, these are Wigner semicircle distributions (are they?)
# conditionally, they are uniform

# plot thedraws
plot(x, pch=20, cex=0.5, xlim=c(-2,2), ylim=c(-2,2))

# squeeze into an ellips
A = matrix(c(2,0,0,0.5), 2, 2)
y = x %*% A
points(y, pch=20, cex=0.5, xlim=c(-2,2), ylim=c(-2,2), col='red')

# rotate the ellipse
R = matrix(c(cos(pi/4), sin(pi/4), sin(pi/4), -cos(pi/4)), 2, 2)
z = y %*% R
points(z, pch=20, cex=0.5, xlim=c(-2,2), ylim=c(-2,2), col='blue')

# do them both in one step (plotted with a slight shift)
H = matrix(c(2*cos(pi/4), 0.5*sin(pi/4), 2*sin(pi/4), -0.5*cos(pi/4)), 2, 2)
w = x %*% H # i.e., H = A %*% R
points(w+0.01, pch=20, cex=0.5, xlim=c(-2,2), ylim=c(-2,2), col='green')

# How to generalize to higher dimension ellipses?
# Issues with the Gibbs sampling? How long should burn in be?
# Drawing directly from the ellipse?

cov(x)

cov(y)
t(A) %*% cov(x) %*% A

cov(z)
t(R) %*% cov(y) %*% R
t(A %*% R) %*% cov(x) %*% A %*% R




### Draw from n-sphere
rsphere = function(n, d = 2, r = 1, burn){
    d = round(ifelse(d <= 2, 2, d))
    if (missing(burn))
        burn = d^2
    out = matrix(0, d, n+burn)
    for (i in 2:(n+burn)){
        for (j in 1:d){
            t = seq((i-1)*d+j-d+1, (i-1)*d+j-1, by = 1)
            out[(i-1)*d+j] = runif(1, -sqrt(r^2 - sum(out[t]^2)), sqrt(r^2 - sum(out[t]^2)))
            }
        }
    return (tail(t(out), n))
    }

x = rsphere(2000, 15, 1, burn = 100)

library(rgl)

A = matrix(c(1.5,0,0.7,0,0.5,-0.3,0.7,-0.3,1), 3, 3)
y = x %*% A

plot3d(x)
points3d(y, col = 'red')

plot(x[,1:2], pch = 20)

pairs(x, pch = 20, cex = 0.5)

dists = apply(x, 1, function(x) sqrt(sum(x^2)))
mean(dists)

l = function(lambda){
    A = matrix(c(13,12,2,12,13,-2,2,-2,8), 3, 3)
    det(A-lambda*diag(3))
    }

l(0)

U = matrix(c(1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2)), 2, 2)
V = matrix(c(1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(18), -1/sqrt(18), 4/sqrt(18), 2/3, -2/3, -1/3), 3, 3)
S = matrix(c(5, 0, 0, 3, 0, 0), 2, 3)

solve(V)
t(V)

U %*% S %*% t(V)

A = matrix(c(1,1,2,1,2,1,2,1,1),3,3)

b = eigen(A)

b$vec %*% diag(b$val) %*% t(b$vec)
svd(A)
