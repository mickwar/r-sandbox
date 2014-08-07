source("./density.R")

# example
xx = rnorm(10000)
dx = density(xx)

plot(dx)
color.den(dx, -3, -2, "blue")
color.den(dx, -1.5, 0, "green")
color.den(dx, 1, 4, "red")

# multimodal example (using hpd)
n = 50000
xx = double(n)
for (i in 1:n){
    r = runif(1)
    if (r < 0.25)
        xx[i] = rnorm(1, -14, 2)
    if (0.25 <= r & r < 0.6)
        xx[i] = rnorm(1, 0, 4)
    if (r >= 0.6)
        xx[i] = rnorm(1, 14, 2)
    }
dx = density(xx)
int = hpd.mult(xx)

plot(dx)
for (i in 1:(length(int)/2))
    color.den(dx, int[2*i-1], int[2*i], "cyan")
