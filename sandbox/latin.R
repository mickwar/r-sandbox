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

random.sequential = function(n, res.max = 100){
    out = matrix(0, n, n)
    total.iter = 0
    reset = FALSE
    reset.counter = 0
    row = 1
    col = 1
    out[1,] = sample(n, n)
    draw = (1:n)[-out[1,1]]
    out[2:n,1] = sample(draw, length(draw))
    while (row < n && col < n){
        total.iter = total.iter + 1
        # restart
        if (reset){
            out = matrix(0, n, n)
            row = 1
            col = 1
            out[1,] = sample(n, n)
            draw = (1:n)[-out[1,1]]
            out[2:n,1] = sample(draw, length(draw))
            reset = FALSE
            }
        # randomize column
        if ((row + col) %% 2){
            col = col + 1
            draw = (1:n)[-out[1:row, col]]
            out[(row+1):n, col] = sample(draw, length(draw))
            reset.counter = 0
            while (reset.counter < res.max &&
                any(out[(row+1):n, col] == out[(row+1):n, 1:(col-1)])){
                reset.counter = reset.counter + 1
                out[(row+1):n, col] = sample(draw, length(draw))
                }
        # randomize row
        } else {
            row = row + 1
            draw = (1:n)[-out[row, 1:col]]
            out[row, (col+1):n] = sample(draw, length(draw))
            reset.counter = 0
            while (reset.counter < res.max && 
                any(out[row, (col+1):n] == t(out[1:(row-1), (col+1):n]))){
                reset.counter = reset.counter + 1
                out[row, (col+1):n] = sample(draw, length(draw))
                }
            }
        if (reset.counter == res.max)
            reset = TRUE
        }
    print(total.iter)
    return (out)
    }
random.sequential(15, res.max = 500)

random.add = function(n, print = FALSE){
    get.inds = function(i){
        col = (i - 1) %/% n + 1
        row = i - n*(col-1)
        # first n-1 are in the row, last n-1 are in the column
        out = c((0:(n-1))*n + row, ((col-1)*n + 1):(col * n))
        return(out[-which(out == i)])
        }
    LS = matrix(0, n, n)
    space = 1:(n^2)
    draw.space = matrix(rep(1:n, n^2), n^2, n, byrow = TRUE)
    while (length(space) > 0){

        if (print)
            print("next iteration")

        lens = 1
        if (length(space) > 1)
            lens = apply(draw.space[space,], 1, function(x) sum(x != 0))
        add.ind = space[which(lens == min(lens))]

        if (print){
            print("add.ind1")
            print(add.ind)
            }

        if (length(add.ind) > 1)
            add.ind = sample(add.ind, 1)
        temp.draw = draw.space[add.ind,]

        if (print){
            print("temp.draw")
            print(temp.draw)
            }

        add.val = temp.draw[temp.draw != 0]
        if (length(add.val) > 1)
            add.val = sample(temp.draw[temp.draw != 0], 1)
        LS[add.ind] = add.val

        if (print){
            print("add.ind2")
            print(add.ind)
            print("add.val")
            print(add.val)
            }
        LS = solver(LS)
        a.temp = table(c(space, which(c(LS) == 0)))
        for (j in as.numeric(names(a.temp[a.temp == 1])))
            draw.space[get.inds(j), LS[j]] = 0
        space = which(c(LS) == 0)

        if (print){
            print("LS")
            print(LS)
#           print("a.temp")
#           print(a.temp)
            }

        }
    return (LS)
    }

### example
LS = gen.full.square(8, as.char = FALSE)
(puzzle = random.remove(LS)); sum(puzzle == 0) / nrow(LS)^2
solver(puzzle)

system.time(random.add(10))
random.add(8, TRUE)

gen.full.square(5)


### simulation
reps = 10
ns = seq(3, 20)
time.gen = matrix(0, reps, length(ns))
for (i in 1:length(ns))
    for (j in 1:reps)
        time.gen[j,i] = system.time(gen.full.square(ns[i]))[3]

time.add = matrix(0, reps, length(ns))
for (i in 1:length(ns))
    for (j in 1:reps)
        time.add[j,i] = system.time(random.add(ns[i]))[3]

mean.gen = apply(time.gen, 2, median)
min.gen = apply(time.gen, 2, min)
max.gen = apply(time.gen, 2, max)

mean.add = apply(time.add, 2, median)
min.add = apply(time.add, 2, min)
max.add = apply(time.add, 2, max)

plot(ns, mean.gen, type='l',# ylim=c(0, max(c(max.gen, max.add))),
    col = 'blue', lwd=3, ylim=c(0,200))
points(ns, min.gen, type='l', lty=2, col = 'blue')
points(ns, max.gen, type='l', lty=2, col = 'blue')
points(ns, mean.add, type='l', lwd=3, col = 'red')
points(ns, min.add, type='l', lty=2, col = 'red')
points(ns, max.add, type='l', lty=2, col = 'red')

