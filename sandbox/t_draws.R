# draws from a t-dstribution
# from robert and casella "monte carlo statistical methods", p.131

make.sig = function(n){
    out = matrix(0, n, n)
    xseq = 1:n
    for (i in 1:n){
        out[i,] = 1/xseq
        xseq[1:i] = xseq[1:i] + 1
        }
    return (out)
    }

S = make.sig(8)
A = solve(S)
(A = ifelse(abs(A) < 1, 0, A))

ft = function(x, nu, theta, sigma)
    (gamma((nu+1)/2)/gamma(nu/2))/sqrt(sigma^2*nu*pi)*
    (1+((x-theta)^2)/(nu*sigma^2))^(-(nu+1)/2)
xx = seq(mu - 5*sqrt(sig2), mu + 5*sqrt(sig2), length=1000)

mu = 0
sig2 = 3
nu = 7

n = 10000
Yinv = rgamma(n, nu/2, nu/2)
Xt = rnorm(n, mu, sqrt(sig2/Yinv))

plot(density(Xt), xlim=c(mu-5*sqrt(sig2),mu+5*sqrt(sig2)))
lines(xx, ft(xx, nu, mu, sqrt(sig2)), col='red')
