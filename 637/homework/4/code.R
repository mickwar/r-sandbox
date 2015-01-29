dat = sleep
logit = function(x)
    log(x / (1-x))
logistic = function(x)
    exp(x) / (1 + exp(x))

### part 1
# proportion of increases found
m = 10
y1 = sum(dat[dat[,2] == 1, 1] > 0)
y2 = sum(dat[dat[,2] == 2, 1] > 0)
phat1 = mean(dat[dat[,2] == 1, 1] > 0) # y1 / m
phat2 = mean(dat[dat[,2] == 2, 1] > 0) # y2 / m

phat1 + c(-1, 1)*qnorm(0.975) * sqrt(0.5 * (1-0.5) / m)
phat2 + c(-1, 1)*qnorm(0.975) * sqrt(0.5 * (1-0.5) / m)

phat1 + c(-1, 1)*qnorm(0.975) * sqrt(phat1 * (1-phat1) / m)
phat2 + c(-1, 1)*qnorm(0.975) * sqrt(phat2 * (1-phat2) / m)

psquigs1 = (phat1 * 10 + 2)/(m + 4)
psquigs2 = (phat2 * 10 + 2)/(m + 4)
psquigs1 + c(-1, 1)*qnorm(0.975) * sqrt(psquigs1 * (1-psquigs1) / (m+4))
psquigs2 + c(-1, 1)*qnorm(0.975) * sqrt(psquigs2 * (1-psquigs2) / (m+4))

eta1 = log(y1 / (m - y1)) + c(-1, 1)*qnorm(0.975)*sqrt(1/y1 + 1/(m-y1))
eta2 = log(y2 / (m - y2)) + c(-1, 1)*qnorm(0.975)*sqrt(1/y2 + 1/(m-y2))
logistic(eta1)
logistic(eta2)

### part 3
smooth = function(x, y, m, d){
    # m = number of points to do the smoothing,
    # should be greater than length(x) to look nice
    if (missing(m))
        m = max(100, 8*length(x))
    if (missing(d))
        d = diff(range(x)) / length(x)

    # gaussian kernel, delta is the bandwidth parameter
    kern = function(x, y)
        exp(-1/(2*d^2)*(x-y)^2)

    # go outside the range of the data by some
    w = diff(range(x))*0.25

    outx = seq(min(x) - w, max(x) + w, length=m)
    outy = double(m)
    # loop to compute the weighted value at each new x
    for (i in 1:length(outy))
        outy[i] = sum(y*kern(x, outx[i]))/
            sum(kern(x, outx[i]))
    return (list("x"=outx, "y"=outy))
    }

li = c(8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 34, 38)
nc = c(2, 2, 3, 3, 3, 1, 3, 2, 1, 1, 1, 1, 1, 3)
nr = c(0, 0, 0, 0, 0, 1, 2, 1, 0, 1, 1, 0, 1, 2)

ss = smooth(li, nr/nc, d = 4)

pdf("links.pdf", height = 12, width = 6)
# logit link
par(mfrow=c(3,1), mar = c(4.1, 4.2, 3.1, 2.1))
mod = glm(cbind(nr, nc-nr) ~ li, family = binomial(link = "logit"))
mod.pred = predict(mod, se.fit = TRUE)
plot(li, logistic(mod.pred$fit), pch = 20, lwd = 5, ylim = c(0, 1),
    main = "Logit Link", xlab = "LI", ylab = quote(widehat("p")))
polygon(c(li, li[length(li):1]), c(logistic(mod.pred$fit - 1.96*mod.pred$se.fit), logistic(mod.pred$fit + 1.96*mod.pred$se.fit)[length(li):1]), col = 'lightgreen', border = "white")
points(li, logistic(mod.pred$fit), pch = 20, lwd = 8, ylim = c(0, 1), col = rgb(0, 0.5, 0))
lines(ss, lty=2, lwd=3)
summary(mod)

# probit link
mod = glm(cbind(nr, nc-nr) ~ li, family = binomial(link = "probit"))
mod.pred = predict(mod, se.fit = TRUE)
plot(li, logistic(mod.pred$fit), pch = 20, lwd = 5, ylim = c(0, 1),
    main = "Probit Link", xlab = "LI", ylab = quote(widehat("p")))
polygon(c(li, li[length(li):1]), c(logistic(mod.pred$fit - 1.96*mod.pred$se.fit), logistic(mod.pred$fit + 1.96*mod.pred$se.fit)[length(li):1]), col = 'lightgreen', border = "white")
points(li, logistic(mod.pred$fit), pch = 20, lwd = 8, ylim = c(0, 1), col = rgb(0, 0.5, 0))
lines(ss, lty=2, lwd=3)
summary(mod)

# cloglog link
mod = glm(cbind(nr, nc-nr) ~ li, family = binomial(link = "cloglog"))
mod.pred = predict(mod, se.fit = TRUE)
plot(li, logistic(mod.pred$fit), pch = 20, lwd = 5, ylim = c(0, 1),
    main = "C-log-log Link", xlab = "LI", ylab = quote(widehat("p")))
polygon(c(li, li[length(li):1]), c(logistic(mod.pred$fit - 1.96*mod.pred$se.fit), logistic(mod.pred$fit + 1.96*mod.pred$se.fit)[length(li):1]), col = 'lightgreen', border = "white")
points(li, logistic(mod.pred$fit), pch = 20, lwd = 8, ylim = c(0, 1), col = rgb(0, 0.5, 0))
lines(ss, lty=2, lwd=3)
summary(mod)
dev.off()

# go with logit link
mod = glm(cbind(nr, nc-nr) ~ li, family = binomial(link = "logit"))
summary(mod)
exp(coef(mod))

mod.pred = predict(mod, se.fit = TRUE)
logistic(mod.pred$fit)

# LI[10] is closest to p = 0.5
li[10]

# OR
coef(mod)[1] / coef(mod)[2] * (-1)

### CI for LI = 8
eta = coef(mod)[1] + coef(mod)[2] * li
cbind(logistic(eta - qnorm(0.975) * mod.pred$se.fit),
    logistic(eta + qnorm(0.975) * mod.pred$se.fit))

### s.e.
# sqrt(1/(m*phat*(1-phat))) ## not the same that R gives

# bootstrap
ps = double(10000)
for (i in 1:length(ps)){
    train = sample(length(nr), replace = TRUE)
    mod = tryCatch(glm(cbind(nr, nc-nr) ~ li, family = binomial(link = "logit"), subset = train), warning = function(x) 0)
    if (is(mod)[1] == "glm"){
        ps[i] = coef(mod)[1] + coef(mod)[2] * 8
    } else {
        ps[i] = 0
        }
    }
ps = ps[ps != 0]
var(ps)
plot(density(ps))
