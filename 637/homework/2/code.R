dat = read.table("~/files/data/637/aids.txt", header = TRUE)
y = as.numeric(t(matrix(unlist(dat[,-1]), 5, 4)))
X = cbind(1, log(1:length(y)))


### part a
pdf("figs/data_plot.pdf", width = 9, height = 9)
plot(y, pch=20, xaxt = "n", ylab = "Number of Cases of AIDS",
    xlab = "Time period", lwd=5, ylim=c(0,max(y)+10))
abline(v=c(4.5, 8.5, 12.5, 16.5), lty=2)
axis(1, at = 1:20, labels = paste0(rep("Q", 20), 1:4))
for (i in 4:8)
    text(4*(i-3)-2+0.5, max(y)+10, paste0("198", i), cex = 2.5)
dev.off()

### part b
pdf("figs/log_data.pdf", width = 9, height = 9)
plot(log(1:20), log(y), pch=20, xaxt = "n", ylab = "log Number of Cases of AIDS",
    xlab = "log Time period")
abline(v=log(c(4.5, 8.5, 12.5, 16.5)), lty=2)
axis(1, at = log(1:20), labels = paste0(rep("Q", 20), 1:4))
dev.off()

### part c
beta = c(0, 1)  # starting point
eps = 1e-6      # stopping rule
diff = eps + 1
iter = 0

while (diff > eps){
    iter = iter + 1
    oldbeta = beta
    xbeta = X %*% oldbeta
    W = diag(as.numeric(exp(xbeta)), length(y))
    z = xbeta + y * exp(-xbeta) - 1
    beta = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
    diff = sqrt(sum((beta - oldbeta)^2))
    if (TRUE){
        print(iter)
        print(beta)
        cat("\n")
        }
    }

### part d
summary(glm(y ~ 0 + X, family = poisson))

### part e
# covariance matrix (using the latest W matrix)
solve(t(X) %*% W %*% X)

# standard errors
(se = sqrt(diag(solve(t(X) %*% W %*% X))))

# z-statistics for beta[2] (slope parameter, for x_i = log(i))
beta[2] / se[2]

# confidence intervals
cbind(beta - qnorm(0.975) * se, beta + qnorm(0.975) * se)
