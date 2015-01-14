dat = read.table("~/files/data/637/aids.txt", header = TRUE)
y = as.numeric(t(matrix(unlist(dat[,-1]), 5, 4)))


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
:w
