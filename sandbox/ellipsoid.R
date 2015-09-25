# draws from the unit circle (Gibbs)

n = 10000
p = 2
r = 1
x = matrix(0, n, p)

for (i in 2:n){
    x[i, 1] = runif(1, -sqrt(r^2 - x[i-1,2]^2), sqrt(r^2 - x[i-1,2]^2))
    x[i, 2] = runif(1, -sqrt(r^2 - x[i,1]^2),   sqrt(r^2 - x[i,1]^2))
    }

# marginally, these are Wigner semicircle distributions
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
