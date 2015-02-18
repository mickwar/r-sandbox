solver = function(LS){
    old.LS = LS
    n = nrow(LS)

    # check if any value shows up multiple times in one row/column
    check = function(LS){
        n = nrow(LS)
        for (i in 1:n){
            # by row
            if (any(table(LS[i, LS[i,] != 0]) > 1))
                return ("error")
            if (any(table(LS[LS[,i] != 0, i]) > 1))
                return ("error")
            }
        return ("good")
        }

    # get the subsquare (the part of the LS that does not have index in
    # any row or column)
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

    # solve the subsquare
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

    # looping mechanism
    while (sum(LS == 0) > 0){
        for (i in 1:n){
            # get the subsquare, the subset of rows and columns with no i
            sub = sub.sq(which(LS == i))
            LS = fill.sq(sub, i)
            check(LS)
            if (check(LS) == "error"){
#               stop("Not a valid latin square.")
                return (0)
                }
            }
        # checking whether the square is solvable
        if (all(old.LS == LS)){
#           stop("Latin square is not solvable or is not a latin square.")
            return (LS)
            }
        old.LS = LS
        }
    return (LS)
    }

LS = matrix(0, 5, 5)
LS[1,1] = 1
#LS[2,2] = 1
LS[4,3] = 2
LS[3,4] = 2
LS[4,4] = 3
LS[1,2] = 4
LS[5,2] = 5
LS
(LS = solver(LS))
LS[4,1] = 5
(LS = solver(LS))
LS[2,3] = 5
solver(LS)

gen.full.square = function(n){
    if (n > 16){
        n = 16
        warning("Computation is too much for that n.")
        }
    char = 1:n
    out = matrix(char[1], n, n)
    out[,1] = sample(char)
    for (i in 2:(n-1))
        while (any(out[,i] == out[,1:(i-1)]))
            out[,i] = sample(char)

    # the last column only has one possibility.
    out[,n] = 0
    out = solver(out)

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

random.add = function(n){
    check = function(LS){
        n = nrow(LS)
        for (i in 1:n){
            # by row
            if (any(table(LS[i, LS[i,] != 0]) > 1))
                return ("error")
            if (any(table(LS[LS[,i] != 0, i]) > 1))
                return ("error")
            }
        return ("good")
        }
    get.inds = function(i){
        col = (i - 1) %/% n + 1
        row = i - n*(col-1)
        unique(c((0:(n-1))*n + row, ((col-1)*n + 1):(col * n)))
        }
    LS = matrix(0, n, n)
    space = 1:(n^2)
    draw.space = rep(list(c(1:n)), n^2)
    while (length(space) > 0){
        add.ind = sample(space, 1)
        add.val = sample(draw.space[[add.ind]], 1)
        LS[add.ind] = add.val
        if (check(LS) == "good"){
            if (all(solver(LS) == 0)){
                LS[add.ind] = 0
                draw.space[[add.ind]] = draw.space[[add.ind]][draw.space[[add.ind]] != add.val]
            } else {
                old.LS = LS
                LS = solver(LS)
                # additional squares may be added, these need to be accounted for
                # as well when reducing the draw space
                new.ind = which(c(old.LS) != c(LS))
                new.val = LS[new.ind]
                if (any(c(LS) == 0) && length(new.ind) > 0){
                    for (k in 1:length(new.ind))
                        for (j in get.inds(new.ind[k]))
                            draw.space[[j]] = draw.space[[j]][draw.space[[j]] != new.val[k]]
                    }

                # further reduce draw space when that value appears in a row or column
                for (j in get.inds(add.ind))
                    draw.space[[j]] = draw.space[[j]][draw.space[[j]] != add.val]
                }
        } else {
            LS[add.ind] = 0
            draw.space[[add.ind]] = draw.space[[add.ind]][draw.space[[add.ind]] != add.val]
            }
        space = which(c(LS) == 0)
        }
    return (LS)
    }

### example
LS = gen.full.square(8, as.char = FALSE)
(puzzle = random.remove(LS)); sum(puzzle == 0) / nrow(LS)^2
solver(puzzle)

system.time(random.add(10))
random.add(9)
gen.full.square(9)


### simulation
time.gen = matrix(0, 10, 6)
for (i in 1:ncol(time.gen))
    for (j in 1:nrow(time.gen))
        time.gen[j,i] = system.time(gen.full.square(i + 3))[3]

time.add = matrix(0, 10, 6)
for (i in 1:ncol(time.add))
    for (j in 1:nrow(time.add))
        time.add[j,i] = system.time(random.add(i + 3))[3]

mean.gen = apply(time.gen, 2, mean)
mean.add = apply(time.add, 2, mean)

plot(4:9, mean.gen)
points(4:9, mean.add)
