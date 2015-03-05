# pred_v was generated from matlab code
v = as.matrix(read.table("./pred_v.txt"))

# get observations
a = v[1,]

g = double(ncol(v))
for (i in 1:ncol(v)){
    g[i] = names(which(quantile(v[,i], seq(0, 1, length=nrow(v))) == a[i]))[1]
    }
g = as.numeric(substring(g, 1, 4)) / 100
g = unique(g)

ks.test(g, punif)
pdf("./figs/ks.pdf")
plot(density(g), main = "Density of the quantiles", xlab = "", lwd = 3, ylab = "")
abline(h = 0)
dev.off()
