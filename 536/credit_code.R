dat = read.csv("../data/credit.csv")
dat = dat[,-1]
dat = dat[-c(84,86,257,262,324),] # outliers/leverage

Y = dat$Balance
X = dat[-ncol(dat)]

k = ncol(X)
n = nrow(X)

# deal with collinearity via orthogonalization of a variable against another
ortho = X[,2] - X[,3] %*% solve(t(X[,3]) %*% X[,3]) %*% t(X[,3]) %*% X[,2] # limit
X[,2] = ortho
names(X)[2] = "Limit*"

(out = step(lm(Y~.+.^2,data=X), k=log(n), direction="both"))
#summary(out)
#plot(out)

k = 13
A = qt(0.975,n-k)*summary(out)[[4]][,2]
A = cbind(A, A)*rep(c(-1,1),each=13)
summary(out)[[4]][,1] + A
out$coef

orig.X = X

X = X[,c(1,2,3,4,5,8)]
X[,6] = ifelse(X[,6]=="Yes",1,0)
X = cbind(X, X[,1]*X[,2], X[,1]*X[,4], X[,1]*X[,6], X[,2]*X[,3],
    X[,3]*X[,4], X[,3]*X[,6])
X = as.matrix(X)

# subset
val.size = floor(3*n/4)
preds = n-val.size
P = ncol(X)+1
reps = 10000
coverage = double(reps)
for (i in 1:reps){
    val.samp = sort(sample(n, val.size))
    predict = (1:n)[-val.samp]
    mod = lm(formula = Y[val.samp] ~ Income + `Limit*` + Rating + Cards + Age + Student + 
        Income:`Limit*` + Income:Cards + Income:Student + `Limit*`:Rating + 
        Rating:Cards + Rating:Student, data = orig.X[val.samp,])
    s2 = (summary(mod)$sigma)^2
    ranges = qt(0.975, val.size-P)*rep(c(-1,1),each=preds)*matrix(
        sqrt(diag(s2*(1+X[predict,] %*% solve(t(X[val.samp,]) %*%
        X[val.samp,])%*%t(X[predict,])))), nrow=preds, ncol=2)
    ests = cbind(1,X[predict,])%*%mod$coef
    ests = cbind(ests, ests)
    intervals = ests + ranges
    intervals = ifelse(intervals<0, 0, intervals)
    within = apply(cbind(Y[predict] >= intervals[,1], Y[predict] <= intervals[,2]), 1, all)
    coverage[i] = mean(within)
    }
mean(coverage)
var(coverage)

pdf("./coverage.pdf")
hist(coverage, col='gray', xlab="Coverage",main="Histogram of Coverage",
    ylab="Density", freq=F)
dev.off()
