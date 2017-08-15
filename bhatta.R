bhatta.dist = function(x, y, ...){
    # Bhattacharyya distance for two continuous random variables x and y
    left = max(min(x), min(y))
    right = min(max(x), max(y))

    if (left > right)
        return (Inf)

    dx = density(x, from = left, to = right, ...)
    dy = density(y, from = left, to = right, ...)

    xx = dx$x
    p = dx$y
    q = dy$y

    BC = sum(mean(diff(xx)) * sqrt(p*q))
    BC = ifelse(BC < 0, 0, BC)
    BC = ifelse(BC >= 1, 1, BC)

    DB = -log(BC)
    
    return (DB)
    }


mu1 = 0
sig1 = 5
mu2 = 0
sig2 = 1

BC.true = sqrt(2*sig1*sig2/(sig1^2+sig2^2))*exp(-(sig1^2+sig2^2)/(4*sig1^2*sig2^2)*( (mu1^2*sig2^2+mu2^2*sig1^2)/(sig1^2+sig2^2) - ( (mu1*sig2^2+mu2*sig1^2)/(sig1^2+sig2^2) )^2))
DB.true = -log(BC.true)

m = 4096
n = 100000

x = rnorm(n, mu1, sig1)
y = rnorm(n, mu2, sig2)

 x = rt(n, sig1)
#y = rt(n, sig2)

#x = (-log(runif(n)))^(-1/4)
#y = 1+(-log(runif(n)))^(-1/8)

dx = density(x, from = min(x, y), to = max(x, y), n = m)
dy = density(y, from = min(x, y), to = max(x, y), n = m)

xx = dx$x
p = dx$y
q = dy$y

par(mfrow = c(2,1))
plot(xx, p, type = 'l', lwd = 3)
lines(xx, q, col = 'blue', lwd = 3)

BC1 =      sum(mean(diff(xx)) * sqrt(p*q))
DB1 =  -log(sum(mean(diff(xx)) * sqrt(p*q)))

dx = density(x, from = max(min(x), min(y)), to = min(max(x), max(y)), n = m)
dy = density(y, from = max(min(x), min(y)), to = min(max(x), max(y)), n = m)

xx = dx$x
p = dx$y
q = dy$y

lines(xx, p, col = 'red')
lines(xx, q, col = 'green')

plot(xx, sqrt(p*q), type = 'l', col = 'orange', lwd= 3)

DB2 = bhatta.dist(x, y, n = m)
BC2 = exp(-DB2)


### Monte Carlo
BC3 = mean(sqrt(dnorm(y, mu1, sig1) / dnorm(y, mu2, sig2)))
DB3 = -log(BC3)

BC4 = mean(sqrt(dnorm(x, mu2, sig2) / dnorm(x, mu1, sig1)))
DB4 = -log(BC4)

BC5 = mean(sqrt(dnorm(y, mean(x), sd(x)) / dnorm(y, mean(y), sd(y))))
DB5 = -log(BC5)

BC6 = mean(sqrt(dnorm(x, mean(y), sd(y)) / dnorm(x, mean(x), sd(x))))
DB6 = -log(BC6)

abs(cbind("BC"=c(BC1, BC2, BC3, BC4, BC5, BC6)-BC.true, "DB"=c(DB1, DB2, DB3, DB4, DB5, DB6)-DB.true))
cbind("BC"=c(BC1, BC2, BC3, BC4, BC5, BC6, BC.true), "DB"=c(DB1, DB2, DB3, DB4, DB5, DB6, DB.true))
