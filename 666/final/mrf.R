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

    # assumes numeric for now (binary will also work)
    # assumes euclidaen distance
    split = function(var, which.nodes){
        # vector of the terminal node indicators
        t.vec = which.nodes
        for (t in t.vec){
            mu = apply(y, 2, mean)
            V = diag(q)
            SSt = 0
            for (i in 1:n)
                SSt = SSt + t(y[i,] - mu) %*% solve(V) %*% (y[i,] - mu)
        
        # where to try the splits
        s = sort(unique(var))
        s = diff(s)/2 + s[-length(s)]


        SStL.1 = 0
        muL.1 = apply(y[tL,], 2, mean)
        for (i in tL)
            SStL.1 = SStL.1 + t(y[i,] - muL.1) %*% solve(V) %*% (y[i,] - muL.1)

        SStR.1 = 0
        muR.1 = apply(y[tR,], 2, mean)
        for (i in tR)
            SStR.1 = SStR.1 + t(y[i,] - muR.1) %*% solve(V) %*% (y[i,] - muR.1)

(       phi1.1 = SSt.1 - SStL.1 - SStR.1)

        for (i in 1:length(s)){
            tL = which(var <= s[i])
            tR = which(var > s[i])
            SSt = 
            }
        }

    }
