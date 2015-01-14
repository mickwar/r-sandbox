dat = read.table("~/files/data/637/aids.txt", header = TRUE)
y = as.numeric(t(matrix(unlist(dat[,-1]), 5, 4)))
X = cbind(1, log(1:length(y)))


### part a
plot(y, pch=20, xaxt = "n", ylab = "Number of Cases of AIDS",
    xlab = "Time period")
abline(v=c(4.5, 8.5, 12.5, 16.5), lty=2)
axis(1, at = 1:20, labels = paste0(rep("Q", 20), 1:4))
for (i in 4:8)
    text(4*(i-3)-2, max(x), paste0("198", i), cex = 2.5)

### part b
plot(log(1:20), log(y), pch=20, xaxt = "n", ylab = "Number of Cases of AIDS",
    xlab = "Time period")
abline(v=log(c(4.5, 8.5, 12.5, 16.5)), lty=2)
axis(1, at = log(1:20), labels = paste0(rep("Q", 20), 1:4))

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
    z = xbeta + (y - exp(xbeta)) * exp(-xbeta)
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
