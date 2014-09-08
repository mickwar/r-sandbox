### 1

f = function(y, sig)
    0.5 * dnorm(y, 66, sig) + 0.5 * dnorm(y, 70, sig)

sig = 3
xx = seq(66 - 3 * sig, 70 + 3 * sig, length=100)
pdf("fig1a.pdf")
plot(xx, f(xx, sig), type='l', ylab="f(y)", xlab="y", main="1(a)")
dev.off()

dnorm(66, 66, sig) * 0.5 / f(66, sig)

pos.1 = function(y, sig)
    0.5 * dnorm(y, 66, sig) / f(y, sig)

pdf("fig1c.pdf")
ss = seq(0.05, 7, length=100)
plot(ss, pos.1(64, ss), type='l', ylim=c(0,1), xlab=expression(sigma),
    ylab=expression(paste("p(",theta,"=1)")), main = "1(c)")
lines(ss, pos.1(66, ss), type='l', col='red')
lines(ss, pos.1(68, ss), type='l', col='blue')
lines(ss, pos.1(70, ss), type='l', col='darkgreen')
lines(ss, pos.1(72, ss), type='l', col='darkorange')
legend(6, 1, col=c("black", "red", "blue", "darkgreen", "darkorange"),
    legend=c(64, 66, 68, 70, 72), lty=1, title="Height")
dev.off()

cross = function(...){
    vec = list(...)
    d = length(vec)
    N = double(d+2) + 1
    for (i in 1:d)
        N[i+1] = length(vec[[i]])
    out = matrix(0, prod(N), d)
    for (i in 1:d){
        out[,i] = rep(vec[[i]], times=prod(N[1:i]),
            each=prod(N[(i+2):(d+2)]))
        }
    return(out)
    }

X = cross(xx, ss)

Z = double(nrow(X))
for (i in 1:nrow(X))
j   Z[i] = pos.1(X[i,1], X[i,2])

library(rgl)
plot3d(cbind(X, Z))

### 2
# bernoulli trials, follow binomial distribution, assuming
# the null is true: p=0.5
pbinom(12, 17, 0.5, lower.tail = FALSE) +
    pbinom(17-13, 17, 0.5, lower.tail = TRUE)

pbinom(28, 44, 0.5, lower.tail = FALSE) +
    pbinom(15, 44, 0.5, lower.tail = TRUE)

X1 = 0:17
X2 = 0:27

X = cross(X1, X2)

f = function(x){
    x1 = x[1]
    x2 = x[2]
    dbinom(x1, 17, 0.5) * dbinom(x2, 27, 0.5)
    }

joint = apply(X, 1, f)

Y = cbind(X, joint)

calc = 0
for (i in 1:nrow(Y)){
    if (Y[i,1] + Y[i,2] <= 15)
        calc = calc + Y[i,3]
    if (Y[i,1] + Y[i,2] >= 29)
        calc = calc + Y[i,3]
    }

p = pbinom(12, 17, 0.5, lower.tail = FALSE) +
    pbinom(17-13, 17, 0.5, lower.tail = TRUE)

p +    calc * (1-p)


f = function(x1)
    (1 - pbinom(28-x1, 27, 0.5) + pbinom(15-x1, 27, 0.5)) * dbinom(x1, 17, 0.5)
p + sum(f(5:12))

### 3
simulate = function(){
    # arrivals (shouldn't be more than 100 patients to arrive
    # before closing time)
    x = rexp(100, 1/10)
    x = cumsum(x[which(cumsum(x) <= 7*60)])

    # length of doctor visit for each patient
    visit = runif(length(x), 5, 20)
    wait = double(length(x))
    n.doc = 3
    j = 1:n.doc

    # compute wait time
    # no matter what, the first n.doc patients won't have to wait
    for (i in (n.doc+1):length(x)){
        temp = x[j] + visit[j] + wait[j]
        j[which.min(temp)] = i
        wait[i] = max(min(temp) - x[i], 0)
        }

    npatients = length(x)
    nwait = sum(wait != 0)
    meanwait = 0
    if (nwait > 0)
        meanwait = mean(wait[wait != 0])
#   meanwait = mean(wait)
    # the clinic stays open at least until 4 p.m.
    close = 420
    for (i in 1:length(x))
        close = max(close, x[i] + visit[i] + wait[i])
    return (c(npatients, nwait, meanwait, close))
    }

set.seed(1)
simulate()

nrep = 1000
out = matrix(0, nrep, 4)

set.seed(1)
for (i in 1:nrep)
    out[i,] = simulate()

nice = apply(out, 2, quantile, c(0.1, 0.5, 0.9))
colnames(nice) = c("patients", "waited", "mean_wait", "closed")
nice

plot(out[,2], out[,3], pch=20)
