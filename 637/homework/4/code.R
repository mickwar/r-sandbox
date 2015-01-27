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

### part 2
