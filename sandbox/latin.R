solver = function(LS){
    old.LS = LS
    n = nrow(LS)
    sub.sq = function(index){
        if (length(index) == 0)
            return (list("sq"=LS, "indx"=NULL))
        out = matrix(0, length(index), 2)
        for (j in 1:length(index)){
            col = (index[j] - 1) %/% n + 1
            row = index[j] - n*(col-1)
            out[j,] = c(row, col)
            }
        return (list("sq"=as.matrix(LS[-out[,1], -out[,2]]), "indx"=out))
        }

    fill.sq = function(sub, char){
        m = nrow(sub$sq)

        # char already solved
        if (m == 0)
            return (LS)

        for (j in 1:m){
            # by column
            if (sum(sub$sq[,j] == 0) == 1 && all(sub$sq[,j] != char))
                sub$sq[which(sub$sq[,j] == 0), j] = char
            }
        for (j in 1:m){
            # by row
            if (sum(sub$sq[j,] == 0) == 1 && all(sub$sq[j,] != char))
                sub$sq[j, which(sub$sq[j,] == 0)] = char
            }
        if (is.null(sub$indx))
            return (sub$sq)
        LS[-sub$indx[,1], -sub$indx[,2]] = sub$sq
        return (LS)
        }

    while (sum(LS == 0) > 0){
        for (i in 1:n){
            # get the subsquare, the subset of rows and columns with no i
            sub = sub.sq(which(LS == i))
            LS = fill.sq(sub, i)
            }
        # checking whether the square is solvable
        if (all(old.LS == LS))
            stop("Latin square is not solvable or is not a latin square.")
        old.LS = LS
        }
    return (LS)
    }

# LS = matrix(0, 5, 5)
# LS[1,1] = 1
# LS[2,2] = 1
# LS[4,3] = 2
# LS[3,4] = 2
# LS[4,4] = 3
# LS[1,2] = 4
# LS[2,4] = 5
# LS

gen.full.square = function(n){
    if (n > 9){
        n = 9
        warning("Computation is too much for that n.")
        }
    char = 1:n
    out = matrix(char[1], n, n)
    out[,1] = sample(char)
    for (i in 2:(n-1))
        while (any(out[,i] == out[,1:(i-1)]))
            out[,i] = sample(char)


    # the last column only has one possibility.
    for (j in 1:n)


    return (out)
    }

random.remove = function(LS, ntry){
    n = nrow(LS)
    if (missing(ntry))
        ntry = n^2
    while (ntry > 0){
        rm.row = sample(n, 1)
        rm.col = sample(n, 1)
        rm.val = LS[rm.row, rm.col]
        LS[rm.row, rm.col] = 0
        if (class(try(solver(LS), silent = TRUE)) == "try-error"){
            LS[rm.row, rm.col] = rm.val
            ntry = ntry - 1
            }
        }
    return (LS)
    }

LS = gen.full.square(8, as.char = FALSE)
(puzzle = random.remove(LS)); sum(puzzle == 0) / nrow(LS)^2
solver(puzzle)

