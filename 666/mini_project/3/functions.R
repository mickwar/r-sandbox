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
    # proportion of variance explained by this set of eigenvectors
    prop = sum(lam[k]) / sum(lam)
    scores = as.matrix(x) %*% eigS$vectors[,k]
    return (list("p"=prop, "scores"=scores))
    }

# create and cut a tree, returns the hclust object, the vector
# of cluster numbers when cutting, and a count for the groups
# also plots the tree and the first two principal components
# marked by group
mw.tree = function(x, k, scores, method = "ward.D", dist = "euclidean"){
    # x is data
    # k is number of groups to obtain when cutting the tree
    # method is the linkage used in creating the tree (see hclust())
    # dist is for the distance metric to be used (see dist())
    par(mfrow=c(1,2), mar=c(5.1, 5.5, 4.1, 2.1))
#   par(mfrow=c(1,2))
    clust.out = hclust(dist(x, dist), method)
    plot(clust.out, sub="", cex.main=1.5, labels=FALSE, cex.lab = 1.5,
        xlab = "Writing Samples")
    cutree.out = cutree(clust.out, k)
    table.out = table(cutree.out)
    center = matrix(0, length(table.out), 2)
    for (i in 1:length(table.out))
        center[i,] = apply(as.matrix(scores[cutree.out == i, 1:2]), 2, mean)
    plot(scores[,1], scores[,2], col = cutree.out, xlab="PC1", ylab="PC2", cex.lab=1.5,
        pch = as.character(cutree.out), cex = 1.5)
    points(center, pch=as.character(1:length(table.out)), cex = 5, col="dodgerblue")
    par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
#   par(mfrow=c(1,1))
    return (list("cluster"=clust.out, "cutree"=cutree.out, "counts"=table.out))
    }

# get cross-validated error rates using linear/quadratic
# discriminant (analysis?)
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
                        t(x[test.index[j],] - xbar[l,]) + determinant(S[[l]])$modulus[1]
                }
            pred.class[j] = which.min(d)
            }
        error[i] = 1-mean(pred.class == y$cutree[test.index])
        }
    cat("\n")
    mean(error)
    }

# get a cross-validated error rate using k-nearest neighbors
knearest = function(data, mw.tree, kfold, knear, seed = 1){
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

        # get predicted classes for the test set
        pred.class = double(m)
        for (j in 1:m){
            pred.class[j] = as.numeric(names(which.max(table(
                y$cutree[train.index[order(d[j,])[1:knear]]]) /
                n[as.numeric(names(table(y$cutree[train.index[order(d[j,])[1:knear]]])))])))
            }
        error[i] = 1-mean(pred.class == y$cutree[test.index])
        }
    cat("\n")
    mean(error)
    }

# for knearest, calculate error.rates and a range of k values
get.error = function(data, mw.tree, k, kfold){
    if (missing(k))
        k = 1:10
    error.rates = double(length(k))
    for (i in 1:length(k))
        error.rates[i] = knearest(x, mw.tree, kfold = kfold, knear = k[i])
    plot(k, error.rates, type='l')
    return(error.rates)
    }

# check cluster validity
validate = function(data, k, scores, method = "ward.D", dist = "euclidean", seed = 1){
    set.seed(seed)
    x = data
    temp.index = sample(1:nrow(x))
    m = floor(nrow(x)/2)
    A.ind = temp.index[1:m]
    B.ind = temp.index[-(1:m)]
    A = x[A.ind,]
    B = x[B.ind,]

    A.clust = hclust(dist(A, dist), method)
    A.cut = cutree(A.clust, k)
    A.centroids = matrix(0, k, ncol(x))
    for (i in 1:k)
        A.centroids[i,] = apply(A[A.cut == i,], 2, mean)

    B1.clust = hclust(dist(B, dist), method)
    B1.cut = cutree(B1.clust, k)

    B2.cut = double(length(B.ind))
    bdist = as.matrix(dist(rbind(B, A.centroids), dist))[1:m,1:k]
    for (i in 1:length(B.ind))
        B2.cut[i] = which.min(bdist[i,])

    temp.g = order(table(B1.cut))
    temp.cut = B1.cut
    for (i in 1:k)
        temp.cut[B1.cut == temp.g[i]] = i
    B1.cut = temp.cut
    temp.g = order(table(B2.cut))
    temp.cut = B2.cut
    for (i in 1:k)
        temp.cut[B2.cut == temp.g[i]] = i
    B2.cut = temp.cut
    

    par(mfrow=c(2,1), mar=c(2.1,4.1,4.1,2.1))
    plot(scores[B.ind,1], scores[B.ind,2], pch=as.character(B1.cut),
        col=B1.cut, cex=1.5)
    plot(scores[B.ind,1], scores[B.ind,2], pch=as.character(B2.cut),
        col=B2.cut, cex=1.5)
    par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
    
    # returns the error rate, but don't put too much stock in this, the
    # cluster labels won't necessarily agree; look at the positions of
    # the clusters
    return(1 - mean(B1.cut == B2.cut))
    }
