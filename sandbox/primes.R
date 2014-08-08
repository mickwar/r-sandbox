# function to get all prime numbers up to n
get.primes1 = function(start = 2, n = 3){
    out = NULL
    if (start <= 2){
        out = 2
        start = 3
        }
    for (i in start:n){
        if (sum(i %% 1:ceiling(sqrt(i)) == 0) == 1)
            out = c(out, i)
        }
    return (out)
    }

get.primes2 = function(n = 3){
    out = 2
    for (i in 3:n){
        if (all(i %% out != 0))
            out = c(out, i)
        }
    return (out)
    }

get.primes3 = function(n = 3){
    out = 2
    to = 1
    index = 2
    for (i in 3:n){
        if (any(sqrt(i) == out)){
            to = to + 1
            index = c(index, out[to])
            }
        if (all(i %% index != 0))
            out = c(out, i)
        }
    return (out)
    }

get.primes1(2, 100)
get.primes2(100)
get.primes3(100)

system.time(get.primes1(2, 10000))
system.time(get.primes2(10000))
system.time(get.primes3(10000))
