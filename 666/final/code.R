library(ISLR)
library(tree)
library(randomForest)

### Beginning of chapter 8 in the 536 book
dat = Hitters # Hitters is from ISLR package
dat = dat[!is.na(dat$Salary),]
names(dat)


pdf("./figs/ex_tree.pdf", width = 18, height = 9)
par(mfrow=c(1,2), mar=c(5.1,4.6,4.1,1.6))
# figure 8.1
tree.hitters = tree(log(Salary) ~ Hits + Years, data = dat)
tree.hitters = prune.tree(tree.hitters, best = 3)
plot(tree.hitters, lwd=2)
text(tree.hitters, cex=3)

# figure 8.2
plot(dat$Years, dat$Hits, pch=20, xlab = "Years", ylab = "Hits",
    cex.lab = 2.0, col = "darkorange")
text(2, 117.5, expression(R[1]), cex = 3)
text(11.8, 52.5, expression(R[2]), cex = 3)
text(11.8, 170.5, expression(R[3]), cex = 3)
abline(v=4.5, lwd=3, col="darkgreen")
lines(c(4.5, 50), c(117.5, 117.5), lwd=3, col="darkgreen")
dev.off()

