### functions
get.pc.scores = function(x, k = 1, p = NULL){
    # x is the data
    # k is the _index_ of eigenvectors to use
    # if p is given, set k to be the number of eigenvalues which
    # cumulative proportion greater than p (overrides k)
    S = var(x)
    eigS = eigen(S)
    lam = eigS$values
    if (!is.null(p)){
        k = 1:which.max(cumsum(lam) / sum(lam))
        }
    prop = sum(lam[k]) / sum(lam)
    scores = as.matrix(x) %*% eigS$vectors[,k]
    return (list("p"=prop, "scores"=scores))
    }
mw.tree = function(x, k, scores, method = "complete", dist = "euclidean"){
    par(mfrow=c(2,1))
    clust.out = hclust(dist(x, dist), method)
    plot(clust.out)
    cutree.out = cutree(clust.out, k)
    table.out = table(cutree.out)
    center = matrix(0, length(table.out), 2)
    for (i in 1:length(table.out))
        center[i,] = apply(as.matrix(scores[cutree.out == i, 1:2]), 2, mean)
    plot(scores[,1], scores[,2], col = cutree.out,
        pch = as.character(cutree.out), cex = 1.5)
    points(center, pch=as.character(1:length(table.out)), cex = 5, col="gray50")
    par(mfrow=c(1,1))
    return (list("cluster"=clust.out, "cutree"=cutree.out, "counts"=table.out))
    }
classify = function(data, mw.tree, method, k = 5, seed = 1){
    # method is one of: linear, quadratic, quadnormal
    set.seed(seed)
    x = data
    y = mw.tree
    n = y$counts
    g = length(n)

    # get index groups for each variable for each sample k
    index.k = matrix(0, g, k+1)
    for (i in 1:g)
        for (j in 2:(k+1))
            index.k[i,j] = index.k[i,j-1] + round((n[i] - index.k[i,j-1])/(k+2-j))
    
    # randomize the samples
    groups = NULL
    for (i in 1:g)
        groups[[i]] = sample(which(y$cutree == i))

    error = double(k)
    for (i in 1:k){
        cat("\rFold:",i,"/",k)   
        # retrieve indices for test and training cases
        test.index = NULL
        train.index = NULL
        for (j in 1:g){
            test.index = c(test.index, groups[[j]][(index.k[j,i]+1):index.k[j,i+1]])
            train.index = c(train.index, groups[[j]][(1:n[j])[-((index.k[j,i]+1):index.k[j,i+1])]])
            }

        # compute xbars and S_i on training set
        if (method == "quadratic" || method == "quadnormal"){
            xbar = matrix(0, g, ncol(x))
            S = NULL
            for (j in 1:g){
                xbar[j,] = apply(x[train.index[y$cutree[train.index] == j],], 2, mean)
                S[[j]] = var(x[train.index[y$cutree[train.index] == j],])
                }
            }
        if (method == "linear")
            Spl = var(x[train.index,])
        # calculate distance and put into class
        pred.class = double(length(test.index))
        d = double(g)
        for (j in 1:length(test.index)){
            for (l in 1:g){
                if (method == "linear")
                    d[l] = as.double(x[test.index[j],] - xbar[l,]) %*% solve(Spl) %*%
                        t(x[test.index[j],] - xbar[l,])
                if (method == "quadratic")
                    d[l] = as.double(x[test.index[j],] - xbar[l,]) %*% solve(S[[l]]) %*%
                        t(x[test.index[j],] - xbar[l,])
                if (method == "quadnormal")
                    d[l] = as.double(x[test.index[j],] - xbar[l,]) %*% solve(S[[l]]) %*%
            pred.class[j] = which.min(d)
            }
        error[i] = 1-mean(pred.class == y$cutree[test.index])
        }
    cat("\n")
    mean(error)
    }
#classify(x, s4.3); classify(x, s4.4); classify(x, s4.5); classify(x, s4.6); classify(x, s4.7)

kmeans = function(data, mw.tree, kfold, knear, seed = 1){
    set.seed(seed)
    k = kfold
    x = data
    nr = nrow(x)
    y = mw.tree
    n = y$counts
    g = length(n)

    # get index groups for each variable for each sample k
    index.k = matrix(0, g, k+1)
    for (i in 1:g)
        for (j in 2:(k+1))
            index.k[i,j] = index.k[i,j-1] + round((n[i] - index.k[i,j-1])/(k+2-j))
    
    # randomize the samples
    groups = NULL
    for (i in 1:g)
        groups[[i]] = sample(which(y$cutree == i))

    error = double(k)
    for (i in 1:k){
        cat("\rFold:",i,"/",k)   
        # retrieve indices for test and training cases
        test.index = NULL
        train.index = NULL
        for (j in 1:g){
            test.index = c(test.index, groups[[j]][(index.k[j,i]+1):index.k[j,i+1]])
            train.index = c(train.index, groups[[j]][(1:n[j])[-((index.k[j,i]+1):index.k[j,i+1])]])
            }
        m = length(test.index)

        # calculate distance matrix
        d = as.matrix(dist(rbind(x[test.index,], x[train.index,])))
        # only consider a partition of the matrix
        d = d[1:m, (m+1):nr]

        pred.class = double(m)
        for (j in 1:m){
            pred.class[j] = as.numeric(names(which.max(table(y$cutree[train.index[order(d[j,])[1:knear]]]))))
            }
        error[i] = 1-mean(pred.class == y$cutree[test.index])
        }
    cat("\n")
    mean(error)
    }


### read in the data
dat = read.table("~/files/R/666/data/collins.txt", header = TRUE)

### remove unnecessary rows
dat = dat[,-which(apply(dat, 2, function(x) length(unique(x))) == nrow(dat))]

### Step 1
# use only first 18 variables
x = dat[,1:18]

# get principal component scores
scores = get.pc.scores(x, k = 1:3)
scores$p                # proportion of variance
scores = scores$scores  # just retain the scores for each observation

# crap
#s1.3 = mw.tree(x, 3, scores, "single"); s1.3$counts
#s1.4 = mw.tree(x, 4, scores, "single"); s1.4$counts
#s1.5 = mw.tree(x, 5, scores, "single"); s1.5$counts
#s1.6 = mw.tree(x, 6, scores, "single"); s1.6$counts
#s1.7 = mw.tree(x, 7, scores, "single"); s1.7$counts

#s2.3 = mw.tree(x, 3, scores, "complete"); s2.3$counts
#s2.4 = mw.tree(x, 4, scores, "complete"); s2.4$counts
#s2.5 = mw.tree(x, 5, scores, "complete"); s2.5$counts
#s2.6 = mw.tree(x, 6, scores, "complete"); s2.6$counts
#s2.7 = mw.tree(x, 7, scores, "complete"); s2.7$counts

#s3.3 = mw.tree(x, 3, scores, "average"); s3.3$counts
#s3.4 = mw.tree(x, 4, scores, "average"); s3.4$counts
#s3.5 = mw.tree(x, 5, scores, "average"); s3.5$counts
#s3.6 = mw.tree(x, 6, scores, "average"); s3.6$counts
#s3.7 = mw.tree(x, 7, scores, "average"); s3.7$counts

s4.3 = mw.tree(x, 3, scores, "ward.D"); s4.3$counts
s4.4 = mw.tree(x, 4, scores, "ward.D"); s4.4$counts
s4.5 = mw.tree(x, 5, scores, "ward.D"); s4.5$counts
s4.6 = mw.tree(x, 6, scores, "ward.D"); s4.6$counts
s4.7 = mw.tree(x, 7, scores, "ward.D"); s4.7$counts

s5.3 = mw.tree(x, 3, scores, "ward.D2"); s5.3$counts
s5.4 = mw.tree(x, 4, scores, "ward.D2"); s5.4$counts
s5.5 = mw.tree(x, 5, scores, "ward.D2"); s5.5$counts
s5.6 = mw.tree(x, 6, scores, "ward.D2"); s5.6$counts
s5.7 = mw.tree(x, 7, scores, "ward.D2"); s5.7$counts

#s6.3 = mw.tree(x, 3, scores, "mcquitty"); s6.3$counts
#s6.4 = mw.tree(x, 4, scores, "mcquitty"); s6.4$counts
#s6.5 = mw.tree(x, 5, scores, "mcquitty"); s6.5$counts
#s6.6 = mw.tree(x, 6, scores, "mcquitty"); s6.6$counts
#s6.7 = mw.tree(x, 7, scores, "mcquitty"); s6.7$counts
#
#s7.3 = mw.tree(x, 3, scores, "median"); s7.3$counts
#s7.4 = mw.tree(x, 4, scores, "median"); s7.4$counts
#s7.5 = mw.tree(x, 5, scores, "median"); s7.5$counts
#s7.6 = mw.tree(x, 6, scores, "median"); s7.6$counts
#s7.7 = mw.tree(x, 7, scores, "median"); s7.7$counts
#
#s8.3 = mw.tree(x, 3, scores, "centroid"); s8.3$counts
#s8.4 = mw.tree(x, 4, scores, "centroid"); s8.4$counts
#s8.5 = mw.tree(x, 5, scores, "centroid"); s8.5$counts
#s8.6 = mw.tree(x, 6, scores, "centroid"); s8.6$counts
#s8.7 = mw.tree(x, 7, scores, "centroid"); s8.7$counts

classify(x, s4.3); classify(x, s5.3)
classify(x, s4.4); classify(x, s5.4)
classify(x, s4.5); classify(x, s5.5)
classify(x, s4.6); classify(x, s5.6)
classify(x, s4.7); classify(x, s5.7)

classify(x, s5.7, "linear", k = 2)
classify(x, s5.7, "linear", k = 3)
classify(x, s5.7, "linear", k = 5)
classify(x, s5.7, "linear", k = 10)

at = 10:30
error.rates = double(length(at))
for (i in 1:length(at))
    error.rates[i] = kmeans(x, s5.7, kfold = 5, knear = at[i])
plot(at, error.rates, type='l')

library(rgl)
plot3d(scores, col=s4.7$cutree, size=5)
