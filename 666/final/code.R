library(ISLR)
library(tree)
library(randomForest)

### Beginning of chapter 8 in the 536 book
dat = Hitters # Hitters is from ISLR package
dat = dat[!is.na(dat$Salary),]
names(dat)


# figure 8.1
tree.hitters = tree(log(Salary) ~ Hits + Years, data = dat, mincut = 45)
plot(tree.hitters)
text(tree.hitters)

# figure 8.2
plot(dat$Years, dat$Hits, pch=20, xlab = "Years", ylab = "Hits",
    cex.lab = 1.5, col = "darkorange")
text(2, 117.5, expression(R[1]), cex = 2)
text(11.8, 52.5, expression(R[2]), cex = 2)
text(11.8, 170.5, expression(R[3]), cex = 2)
abline(v=4.5, lwd=2, col="darkgreen")
lines(c(4.5, 50), c(117.5, 117.5), lwd=2, col="darkgreen")
