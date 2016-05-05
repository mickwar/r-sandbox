### Example of slice sampling for the standard normal
### Requires the inverse pdf (easy for the normal)
n = 10000
y = double(1000)
for (i in 2:n){
    h = runif(1, 0, dnorm(y[i-1]))
    y[i] = runif(1, -sqrt(-2*log(h)-log(2)-log(pi)), sqrt(-2*log(h)-log(2)-log(pi)))
    }
curve(dnorm(x), xlim = c(-4, 4), lwd = 2)
lines(density(y), col = 'red', lwd = 2)
lines(density(rnorm(n)), col = 'blue', lwd = 2)
