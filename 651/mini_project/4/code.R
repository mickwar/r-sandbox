# functions
hpd.uni = function(x, prob = 0.95, precision = 1000){
    range = seq(0, 1-prob, length=precision)
    range = cbind(range, range+prob)
    best = range[which.min(apply(range, 1, function(y)
        diff(quantile(x, y)))),]
    return (quantile(x, best))
    }
color.den = function(dens, from, to, col.inside = 1, col.border = NULL){
    if (is.null(col.border))
        col.border = col.inside
    index = which(dens$x > from & dens$x < to)
    polygon(c(dens$x[index][1], dens$x[index],
        dens$x[index][length(index)]), c(0, dens$y[index], 0),
        col=col.inside, border=col.border)
    }
col.mult = function(col1 = 0x000000, col2 = "black"){
    if (is.character(col1))
        val1 = t(col2rgb(col1) / 255)
    if (is.numeric(col1))
        val1 = t(int2rgb(col1) / 255)
    if (is.character(col2))
        val2 = t(col2rgb(col2) / 255)
    if (is.numeric(col2))
        val2 = t(int2rgb(col2) / 255)
    rgb(val1 * val2)
    }
int2rgb = function(x){
    hex = as.character(as.hexmode(x))
    hex = paste0("#", paste0(rep("0", 6-nchar(hex)), collapse=""), hex)
    col2rgb(hex)
    }
hpd.plot = function(dens, hpd, col1, col2 = NULL, multiply = TRUE, border = "black", ...){
    if (is.null(col2))
        col2 = "gray50"
    if (multiply)
        col2 = col.mult(col1, col2)
    plot(dens, type='n', ...)
    polygon(dens, col=col1)
    color.den(dens, hpd[1], hpd[2], col2)
    lines(dens, col = border)
    }

y = scan("~/files/R/651/data/faculty.dat")

k = length(y)
n = k
ybar = mean(y)

igpdf=function(x, a, b)
    (b^a)/gamma(a)*x^(-a-1)*exp(-b/x)

# The model:
# Y_i ~ N(theta_i, sigma^2)
# theta_i ~ N(mu, tau^2)
# sigma^2 ~ IG(a, b)
# mu ~ N(m, s^2)
# tau^2 ~ IG(c, d)

# priors
# for mu
m = 6       # prior mean
s2 = 0.25^2  # prior variance on the prior mean

# for tau^2 ("between" faculty variance)
c = 3
d = 2
d/(c+1) # mode
d/(c-1) # mean
d^2/((c-1)^2*(c-2)) # var
plot(density(sqrt(1/rgamma(10000, c, scale = d))), xlim=c(0, 2))

# for sigma^2 ("within" faculty variance)
a = 4
b = 1.5
b/(a+1) # mode
b/(a-1) # mean
b^2/((a-1)^2*(a-2)) # var
plot(density(sqrt(1/rgamma(10000, a, scale = b))), xlim=c(0, 2))

nburn = 100
nmcmc = 100000
p.theta = matrix(0, nburn + nmcmc, k)
p.mu = double(nburn + nmcmc)
p.sig2 = double(nburn + nmcmc)
p.tau2 = double(nburn + nmcmc)

# init
p.theta[1,] = m
p.mu[1] = m
p.sig2[1] = 1/(b * (a+1))   # mode for IG
p.tau2[1] = 1/(d * (c+1))

for (i in 2:(nburn+nmcmc)){
    # update thetas
    p.theta[i,] = rnorm(n, (y*p.tau2[i-1] + p.mu[i-1]*p.sig2[i-1]) /
        (p.tau2[i-1] + p.sig2[i-1]), sqrt(p.tau2[i-1]*p.sig2[i-1] /
        (p.tau2[i-1] + p.sig2[i-1])))

    # update mu
    p.mu[i] = rnorm(1, (mean(p.theta[i,])*k*s2 + m*p.tau2[i-1]) /
        (k*s2+p.tau2[i-1]), sqrt(s2*p.tau2[i-1] / (k*s2+p.tau2[i-1])))

    # update tau^2
    p.tau2[i] = 1/rgamma(1, shape = c + k/2, scale = 1/(1/d + 0.5 *
        sum((p.theta[i,] - p.mu[i])^2)))

    # update sigma^2
    p.sig2[i] = 1/rgamma(1, shape = a + n/2, scale = 1/(1/b + 0.5 *
        sum((y - p.theta[i,])^2)))
    }

p.theta = p.theta[(nburn+1):(nburn+nmcmc),]
p.mu = p.mu[(nburn+1):(nburn+nmcmc)]
p.tau2 = p.tau2[(nburn+1):(nburn+nmcmc)]
p.sig2 = p.sig2[(nburn+1):(nburn+nmcmc)]

# display the parameters
pdf("figs/post_theta.pdf")
par(mar = c(5.1, 5.1, 4.1, 2.1))
boxplot(p.theta, pch=20, cex=0.5, main="", ylab=expression(theta[i]), cex.lab=2,
    xlab = "i", axes = FALSE, col='dodgerblue')
box()
axis(1, at = 1:23, labels = FALSE)
axis(2)
mtext(side = 1, line = 1, at = 1:23, text = as.character(1:23))
dev.off()

dens.mu = density(p.mu)
hpd.mu = hpd.uni(p.mu)

dens.tau = density(p.tau2)
hpd.tau = hpd.uni(p.tau2)

dens.sig = density(p.sig2)
hpd.sig = hpd.uni(p.sig2)

pdf("figs/post_other.pdf", width=12, height=4)
par(mfrow=c(1,3))
hpd.plot(dens.mu, hpd.mu, "dodgerblue", main="", ylab="", xlab=expression(mu), cex.lab=2)
hpd.plot(dens.tau, hpd.tau, "dodgerblue", main="", ylab="", xlab=expression(tau^2), cex.lab=2)
hpd.plot(dens.sig, hpd.sig, "dodgerblue", main="", ylab="", xlab=expression(sigma^2), cex.lab=2)
dev.off()

rr = cor(cbind(p.theta, "mu"=p.mu, "sig^2"=p.sig2, "tau^2"=p.tau2))
library(fields)

pdf("figs/post_var.pdf")
image.plot(rr)
box()
dev.off()

library(corrplot)
pdf("figs/post_var.pdf")
corrplot(rr, method="color")
dev.off()

# posterior mean
c(apply(p.theta, 2, mean), mean(p.mu), mean(p.tau2), mean(p.sig2))

# posterior variance
apply(cbind(p.theta, p.mu, p.tau2, p.sig2), 2, var)



# posterior predictive
preds = rnorm(nmcmc, rnorm(nmcmc, p.mu, sqrt(p.tau2)), sqrt(p.sig2))

pdf("figs/pred.pdf")
hist(y, col='gray', freq=FALSE, breaks=6, ylim=c(0, 0.8), xlim=c(3.8,7.2),
    main="", ylab="", xlab="Faculty Rating", cex.lab=2)
points(density(y), col="black", type='l', lwd=3)
points(density(preds), col="chartreuse4", type='l', lwd=5)
dev.off()


mean(preds)
var(preds)
quantile(preds)
mean(preds > 5)



mean(preds > 7)

# "bootstrap"
boot = double(100)
phil.boot = double(100)
for (i in 1:length(boot)){
    preds = rnorm(nmcmc, rnorm(nmcmc, p.mu, sqrt(p.tau2)), sqrt(p.sig2))
    phil.pred = rnorm(n*nmcmc, c(p.theta), sqrt(rep(p.sig2, times = n)))
    boot[i] = mean(preds > 5)
    phil.boot[i] = mean(phil.pred > 5)
    }
mean(boot)
mean(phil.boot)
plot(density(boot), ylim=c(0, 2500))
points(density(phil.boot), type='l', col='red')

plot(density(preds))
points(density(phil.pred), type='l', col='red')

#plot(density(preds), col='red', type='l')
