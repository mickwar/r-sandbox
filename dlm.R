
##### Model M1
library(Matrix)
### Part 1


times = 1:(365*6)
x = sin(2*pi / 365 * times) + 0.5*sin(2*pi / (365/2)*times)
y = x + rnorm(length(times))
plot(times, y, type='l')

y = read.csv("./googletrendsUCSC.csv")[,2]
times = read.csv("./googletrendsUCSC.csv")[,1]

substrRight = function(x, n)
    substr(x, nchar(x)-n+1, nchar(x))

times = unique(c(as.character(times),
    sort(paste0(2017:2024,"-",rep(substrRight(paste0("0", 1:12), 2), each = length(2017:2024))))))

x.ind = round(seq(1, length(y), by = 12))
#as.numeric(times)

Gj = function(j, p)
    matrix(c(cos(2*pi*j/p), -sin(2*pi*j/p), sin(2*pi*j/p), cos(2*pi*j/p)), 2, 2)

F.vec = matrix(c(c(1,0), rep(c(1,0), 12/2 - 1), 1), ncol = 1)
G.mat = as.matrix(bdiag(matrix(c(1,0,1,1),2,2),
    Gj(1, 365), Gj(2, 365), Gj(3, 365), Gj(4, 365), Gj(5, 365), -1))

F.vec = matrix(c(1, 0), ncol = 1)
G.mat = matrix(c(1,0,1,1), 2, 2)

#F.vec = matrix(c(c(1,0), rep(c(1,0), 3), 1), ncol = 1)
#G.mat = as.matrix(bdiag(matrix(c(1,0,1,1),2,2),
#    Gj(1, 12), Gj(2, 12), Gj(5, 12), -1))

# F.vec = matrix(1, 1, 1)
# G.mat = matrix(1, 1, 1)

dat = list("y" = y, "F" = F.vec, "G" = G.mat)

### Part 2
get.vals = function(dat, delta, m0, C0, n0, d0){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    n = length(y)
    p = length(m0)

    Rt = array(0, c(n+1, p, p))
    qt = double(n+1)
    At = matrix(0, n+1, p)
    ft = double(n+1)
    et = double(n+1)
    mt = matrix(0, n+1, p)
    nt = double(n+1)
    dt = double(n+1)
    st = double(n+1)
    Ct = array(0, c(n+1, p, p))

    mt[1,] = m0
    Ct[1,,] = C0
    nt[1] = n0
    dt[1] = d0
    st[1] = dt[1] / nt[1]

    for (i in 2:(n+1)){
#       Rt[i,,] = G.mat %*% Ct[i-1,,] %*% t(G.mat) + (1-delta)/delta*Ct[i-1,,]
        Rt[i,,] = 1/delta * G.mat %*% Ct[i-1,,] %*% t(G.mat)
        qt[i] = t(F.vec) %*% Rt[i,,] %*% F.vec + st[i-1]
        At[i,] = (Rt[i,,] %*% F.vec)/qt[i]   
        ft[i] = t(F.vec) %*% G.mat %*% mt[i-1,]
        et[i] = y[i-1] - ft[i]
        mt[i,] = G.mat %*% mt[i-1,] + At[i,] * et[i]
        nt[i] = nt[i-1] + 1
        dt[i] = dt[i-1] + st[i-1]*(et[i]^2) / qt[i]
        st[i] = st[i-1] + st[i-1] / nt[i]*( (et[i]^2) / qt[i] - 1)
        Ct[i,,] = (st[i] / st[i-1])*(Rt[i,,] - At[i,] %*% t(At[i,]) * qt[i])
        }

    
    return (list("Rt"=Rt[-1,,], "qt"=qt[-1], "At"=At[-1,],
        "ft"=ft[-1], "et"=et[-1], "mt"=mt[-1,], "nt"=nt[-1],
        "dt"=dt[-1], "st"=st[-1], "Ct"=Ct[-1,,], "delta"=delta,
        "m0"=m0, "C0"=C0, "n0"=n0, "d0"=d0))
    }
obs.pred.dens = function(dat, params){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    df = c(params$n0, params$nt[-length(params$nt)])
    mu = params$ft
    sig2 = params$qt
    out = lgamma((df+1)/2) - lgamma(df/2) - 1/2*log(df*pi*sig2) -
        (df+1)/2*log(1+(y - mu)^2 / (df*sig2))
    return (sum(out))
    }
smooth = function(dat, params){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    G.inv = solve(G.mat)
    n = length(y)
    p = length(F.vec)
    delta = params$delta
    out.at = matrix(0, n, p)
    out.Rt = array(0, c(n, p, p))
    out.at[n,] = params$mt[n,]
    out.Rt[n,,] = params$Ct[n,,]
    out.ft = double(n)
    out.qt = double(n)
    out.ft[n] = t(F.vec) %*% out.at[n,]
    out.ft[n] = t(F.vec) %*% out.at[n,]
    for (i in (n-1):1){
        out.at[i,] = (1-delta)*params$mt[i,] + delta*(G.inv %*% out.at[i+1,])
        out.Rt[i,,] = (1-delta)*params$Ct[i,,] + delta^2*(G.inv %*% out.Rt[i+1,,] %*% t(G.inv))
        out.ft[i] = t(F.vec) %*% out.at[i,]
        out.qt[i] = t(F.vec) %*% out.Rt[i,,] %*% F.vec + params$st[n]
        }
    return (list("at.s"=out.at, "Rt.s"=out.Rt, "ft.p"=out.ft, "qt.p"=out.qt))
    }
forecast = function(dat, params, h = 12){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    n = length(y)
    p = length(F.vec)
    delta = params$delta

#   at.p = rep(params$mt[n], h)
#   Rt.p = (1/delta^(1:h))*params$Ct[n]
#   ft.p = rep(params$mt[n], h)
#   qt.p = Rt.p + params$St[n]
#   return (list("at.p"=at.p, "Rt.p"=Rt.p,
#       "ft.p"=ft.p, "Qt.p"=Qt.p))

    at.p = matrix(0, h+1, p)
    Rt.p = array(0, c(h+1, p, p))
    ft.p = double(h+1)
    qt.p = double(h+1)

    at.p[1,] = params$mt[n,]
    Rt.p[1,,] = params$Ct[n,,]
    for (i in 2:(h+1)){
        at.p[i,] = G.mat %*% at.p[i-1,]
        Rt.p[i,,] = G.mat %*% (Rt.p[i-1,,] + (1-delta)/delta*params$Ct[n,,]) %*% t(G.mat)
        ft.p[i] = t(F.vec) %*% at.p[i,]
        qt.p[i] = t(F.vec) %*% Rt.p[i,,] %*% F.vec + params$st[n]
        }
    return (list("at.p"=at.p[-1,], "Rt.p"=Rt.p[-1,,],
        "ft.p"=ft.p[-1], "qt.p"=qt.p[-1]))
    }

# Priors
n = length(dat$y)
p = 2
m0 = double(p)
m0[1] = 0
C0 = 1*diag(p)
n0 = 1
d0 = 10

#p = 1
#m0 = 60
#C0 = 200

# Optimal discount factor
del = seq(0.8, 1, by = 0.01)
opd = double(length(del))
for (i in 1:length(del)){
    out = get.vals(dat, delta = del[i], m0 = m0, C0 = C0, n0 = n0, d0 = d0)
    opd[i] = obs.pred.dens(dat, out)
    }
(d.max = del[which.max(opd)])

pdf("./figs/m1_discount.pdf", height = 8, width = 8)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(del, opd, type='l', ylab = "log observed predictive density", bty = 'n',
    xlab = "Discount Factor" ~ delta, cex.lab = 1.4, cex.main = 2)
abline(v = del[which.max(opd)], lty = 2, col = 'gray', lwd = 3)
points(d.max, max(opd), pch = 16)
text(d.max+0.02, max(opd) + 1.0, round(max(opd), 2), cex = 2.0)
dev.off()

out = get.vals(dat, delta = d.max, m0 = m0, C0 = C0, n0 = n0, d0 = d0)
spar = smooth(dat, out)         # Smoothing
h = 365*5
ppar = forecast(dat, out, h)    # Forecasting

ppar$ft.p = spar$ft.p
ppar$qt.p = spar$qt.p



## 1 (a)
## Marginal filtering distributions for polynomial trend and each harmonic
pdf("./figs/m1_components.pdf", width = 12, height = 16)
par(mfrow = c(3,2), mar = c(5.1, 4.1, 4.1, 2.1))
plot(y, bty='n', type='o', ylim = c(0, 120), col = 'gray50',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = "Trend")
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
center = out$mt[,1]
bound = sqrt(out$Ct[,1,1])*qt(0.975, out$nt)
lines(center, col = 'blue', lwd = 3)
lines(center + bound, col = 'lightblue', lwd = 2, lty = 2)
lines(center - bound, col = 'lightblue', lwd = 2, lty = 2)
for (i in 1:5){
    plot(0, bty='n', type='n', xlim = c(1,n), ylim = c(-20, 20),
        xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
        axes = FALSE, main = paste("Harmonic", i))
    axis(2)
    axis(1, at = x.ind, labels = times[x.ind])
    center = out$mt[,2*i+1]
    bound = sqrt(out$Ct[,2*i+1,2*i+1])*qt(0.975, out$nt)
    lines(center, col = 'blue', lwd = 3)
    lines(center + bound, col = 'lightblue', lwd = 2, lty = 2)
    lines(center - bound, col = 'lightblue', lwd = 2, lty = 2)
    text(n-5, 18, round(mean(center + bound < 0), 2))
    text(n-5, -18, round(mean(center - bound > 0), 2))
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

## One-step ahead
pdf("./figs/m1_onestep.pdf", width = 12, height = 8)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(y, bty='n', type='o', col = 'gray50',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ y[t] ~ "|" ~ D[t-1] )
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(ppar$ft.p, col = 'red', lwd = 3)
lines(ppar$ft.p + sqrt(ppar$qt.p)*qt(0.975, out$nt-1), col = 'pink', lwd = 1, lty = 2)
lines(ppar$ft.p - sqrt(ppar$qt.p)*qt(0.975, out$nt-1), col = 'pink', lwd = 1, lty = 2)
lines(times, x, col = 'blue')
lines(out$ft, col = 'green')
dev.off()

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(exp(y), bty='n', type='o', col = 'gray50',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ y[t] ~ "|" ~ D[t-1] )
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
lines(exp(ppar$ft.p + ppar$qt.p/2), col = 'red', lwd = 3)
lines(qlnorm(0.025, ppar$ft.p, sqrt(ppar$qt.p)), col = 'pink')
lines(qlnorm(0.975, ppar$ft.p, sqrt(ppar$qt.p)), col = 'pink')
(exp(ppar$qt.p/2)-1)*exp(2*ppar$ft.p+ppar$qt.p), col = 'pink')
lines(ppar$ft.p + sqrt(ppar$qt.p)*qt(0.975, out$nt-1), col = 'pink', lwd = 1, lty = 2)
lines(ppar$ft.p - sqrt(ppar$qt.p)*qt(0.975, out$nt-1), col = 'pink', lwd = 1, lty = 2)
lines(times, exp(x), col = 'blue')
lines(exp(out$ft + out$qt/2), col = 'green')

## Marginal smoothing distributions for polynomial trend and each harmonic
pdf("./figs/m1_smooth_components.pdf", width = 12, height = 16)
par(mfrow = c(3,2), mar = c(5.1, 4.1, 4.1, 2.1))
plot(y, bty='n', type='o', col = 'gray50',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = "Trend")
axis(2)
axis(1, at = x.ind, labels = times[x.ind])
center = spar$at.s[,1]
bound = sqrt(out$st[n] / out$st * spar$Rt.s[,1,1]) * qt(0.975, out$nt)
lines(center, col = 'green', lwd = 3)
lines(center + bound, col = 'lightgreen', lwd = 2, lty = 2)
lines(center - bound, col = 'lightgreen', lwd = 2, lty = 2)
for (i in 1:5){
    plot(0, bty='n', type='n', xlim = c(1,n),
        xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
        axes = FALSE, main = paste("Harmonic", i))
    axis(2)
    axis(1, at = x.ind, labels = times[x.ind])
    center = spar$at.s[,2*i+1]
    bound = sqrt(out$st[n] / out$st * spar$Rt.s[,2*i+1,2*i+1]) * qt(0.975, out$nt)
    lines(center, col = 'green', lwd = 3)
    lines(center + bound, col = 'lightgreen', lwd = 2, lty = 2)
    lines(center - bound, col = 'lightgreen', lwd = 2, lty = 2)
    text(n-5, 1, round(mean(center + bound < 0), 2))
    text(n-5, -1, round(mean(center - bound > 0), 2))
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

## 5-year forecast
pdf("./figs/m1_forecast.pdf", width = 12, height = 8)
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(1:(n+h), c(dat$y, rep(NA, h)), bty='n', type='o', ylim = c(0, 120), col = 'gray50',
    xlab = "Time", ylab = "Search interest", cex.lab = 1.4, cex.main = 2,
    axes = FALSE, main = ~ y[T+h] ~ "|" ~ D[T], xlim = c(n-36, n + h))
axis(2)
axis(1, at = seq(1, n+h*12, by = 12), labels = times[seq(1, n+h*12, by = 12)])
abline(v = n, lty=2, col = 'gray')
lines(out$ft, col = 'red', lwd = 4)
lines(out$ft + sqrt(out$qt)*qt(0.975, out$nt-1), col = 'pink', lwd = 1, lty = 2)
lines(out$ft - sqrt(out$qt)*qt(0.975, out$nt-1), col = 'pink', lwd = 1, lty = 2)
lines((n+1):(n+h), ppar$ft.p, col = 'firebrick1', lwd = 3)
lines((n+1):(n+h), ppar$ft.p + sqrt(ppar$qt.p)*qt(0.975, out$nt[n]),
    col = 'pink2', lwd = 1, lty = 2)
lines((n+1):(n+h), ppar$ft.p - sqrt(ppar$qt.p)*qt(0.975, out$nt[n]),
    col = 'pink2', lwd = 1, lty = 2)
dev.off()

