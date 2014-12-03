### Multivariate random forests
### (as explained by Segal and Xiao (2011))
### (for regression)


library(mvpart)
y = spider[,1:12]
x = spider[,13:18]

mrf = function(y, x, mtry, ntree, subset, minnodesize, maxnodesize, V){
    p = ncol(x)

    if (missing(mtry))
        mtry = round(sqrt(p))
    if (missing(subset))
        subset = 1:n

    y = as.matrix(y[subset,])
    n = nrow(y)
    q = ncol(y)

    if (missing(V))
        V = diag(q)

    # initialize nodes, each row indicates a split
    # column 1: node indicator
    # column 2: which variable was split
    # column 3: where the variable was split
    # column 4: parent node
    # column 5: number of observations inside the node
    # column 6: 1 for left, 2 for right
    nodes = matrix(0, 1, 6)
    nodes[1,1] = 1  # first parent node
    nodes[1,5] = n

    get.terminal.nodes = function(nodes){
        if (nrow(nodes) > 1)
            return ((1:nrow(nodes))[-unique(nodes[,4])])
        return (1)
        }

    # get the index of observations in node "which"
    get.node.obs = function(which.node){
        if (which.node == 1)
            return (subset)
        keep.bool = rep(TRUE, nodes[1, 5])
        current = which.node
        while (current > 0){
            if (nodes[current, 6] == 1)
                keep.bool = keep.bool & x[,nodes[current, 2]] <= nodes[current, 3]
            if (nodes[current, 6] == 2)
                keep.bool = keep.bool & x[,nodes[current, 2]] > nodes[current, 3]
            # current is now the parent node
            current = nodes[current, 4]
            }
        return (which(keep.bool))
        }

    # assumes numeric for now (binary will also work)
    # assumes euclidaen distance
    split.node = function(var, which.node){
        t.obs = get.node.obs(which.node)
        temp = matrix(0, length(var), 5)
        temp[,1] = var
        for (v in var){
            # get obs inside

            yv = y[t.obs,]
            xv = x[t.obs, v]
            mu = apply(y[t.obs,], 2, mean)
            V = diag(q)
            SSt = 0
            for (i in t.obs)
                SSt = SSt + t(y[i,] - mu) %*% solve(V) %*% (y[i,] - mu)
        
            # where to try the splits
            s = sort(unique(xv))
            s = diff(s)/2 + s[-length(s)]
            phi = double(length(s))

            for (j in 1:length(s)){
                tL = which(xv <= s[j])
                tR = which(xv > s[j])

                SSL = 0
                muL = apply(matrix(yv[tL,], nrow = length(tL)), 2, mean)
                for (i in tL)
                    SSL = SSL + t(yv[i,] - muL) %*% solve(V) %*% (yv[i,] - muL)

                SSR = 0
                muR = apply(matrix(yv[tR,], nrow = length(tR)), 2, mean)
                for (i in tR)
                    SSR = SSR + t(yv[i,] - muR) %*% solve(V) %*% (yv[i,] - muR)

                phi[j] = SSt - SSL - SSR
                }
            temp[temp[,1] == v, 2] = s[which.max(phi)]                      # the cutoff
            temp[temp[,1] == v, 3] = length(which(xv <= s[which.max(phi)])) # nobs to the left
            temp[temp[,1] == v, 4] = length(which(xv > s[which.max(phi)]))  # nobs to the right
            temp[temp[,1] == v, 5] = max(phi)                               # optimzation functional
            if (is.na(temp[temp[,1] == v, 2]))
                temp[temp[,1] == v, 2:5] = -Inf
            }
        return (temp[which.max(temp[,5]),])
        }

    # todo: need a stopping rule for when to stop growing the tree
    #       min node size?

    # grow the tree
    for (i in 1:10){
        n.ind = get.terminal.nodes(nodes)
        temp = matrix(0, length(n.ind), 6)
        for (j in 1:length(n.ind)){
            temp[j,] = c(split.node(1:p, n.ind[j]), n.ind[j])
            }
        newsplit = temp[which.max(temp[,5]),]
        nodes = rbind(nodes, c(nrow(nodes)+1, newsplit[1], newsplit[2], newsplit[6], newsplit[3], 1))
        nodes = rbind(nodes, c(nrow(nodes)+1, newsplit[1], newsplit[2], newsplit[6], newsplit[4], 2))
        }

    # check index should contain every observation in subset
#   j = get.terminal.nodes(nodes)
#   index = NULL
#   for (i in 1:length(j))
#       index = c(index, get.node.obs(j[i]))
    

    }
