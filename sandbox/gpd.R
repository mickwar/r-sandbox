### Generalized Pareto distribution
dgpd = function(x, mu, sigma, ksi){
    y = 1/sigma * (1 + ksi * (x - mu) / sigma)^(-1/ksi - 1)
    z = 1*(x >= mu)
    if (ksi < 0)
        z = z * (x <= mu - sigma/ksi)
    return (y*z)
    }
    

xx = seq(0, 5, length = 1000)
plot(xx, dgpd(xx, 0, 1, 1), type='l', xlim = c(0, 5), ylim = c(0, 1), col = 'darkblue', lwd=2)
lines(xx, dgpd(xx, 0, 1, 5), col = 'green', lwd=2)
lines(xx, dgpd(xx, 0, 1, 20), col = 'red', lwd=2)
lines(xx, dgpd(xx, 0, 2, 1), col = 'purple', lwd=2)
lines(xx, dgpd(xx, 0, 2, 5), col = 'orange', lwd=2)
lines(xx, dgpd(xx, 0, 2, 20), col = 'lightblue', lwd=2)

