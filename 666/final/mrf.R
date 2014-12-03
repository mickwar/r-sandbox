### Multivariate random forests


library(mvpart)
y = spider[,1:12]
x = spider[,13:18]

mrf = function(y, x, mtry, ntree, subset){
    p = ncol(x)

    if (missing(mtry))
        mtry = round(sqrt(p))
    if (missing(subset))
        subset = 1:n

    y = as.matrix(y[subset,])
    n = nrow(y)
    q = ncol(y)

    # initialize nodes, each row indicates a split
    # column 1: node indicator
    # column 2: which variable was split
    # column 3: where the variable was split
    # column 4: parent node
    # column 5: number of observations inside the number
    nodes = matrix(0, 1, 5)
    nodes[1,1] = 1  # first parent node
    nodes[1,5] = n

    get.terminal.nodes = function(nodes){
        if (nrow(nodes) > 1)
            return ((1:nrow(nodes))[-unique(nodes[,4])])
        return (1)
        }

    # get the index of observations in node "which"
    get.node.obs = function(nodes, which){
        if (which == 1)
            return (subset)
        keep = NULL
        current = which
        while (current > 0){
            keep = c(keep, 1)
            }
        }

    # assumes numeric for now (binary will also work)
    # assumes euclidaen distance
    split.node = function(var, which.node){
        temp = matrix(0, length(var), 5)
        temp[,1] = var
        for (v in var){
            # get obs inside
            t.obs = get.node.obs(nodes, which.node)

            yv = y[t.obs,]
            xv = x[t.obs, v]
            mu = apply(y[t.obs,], 2, mean)
            V = diag(q)
            SSt = 0
            for (i in t.obs)
                SSt = SSt + t(yv[i,] - mu) %*% solve(V) %*% (yv[i,] - mu)
        
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
            temp[temp[,1] == v, 2] = s[which.max(phi)]
            temp[temp[,1] == v, 3] = length(which(xv <= s[which.max(phi)]))
            temp[temp[,1] == v, 4] = length(which(xv > s[which.max(phi)]))
            temp[temp[,1] == v, 5] = max(phi)
            }
        temp[which.max(temp[,5]),]
        }
    

    }
