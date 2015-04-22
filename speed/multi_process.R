trailing = function(x, digits = 4)
    formatC(x, digits=digits, format="f")

nice.time = function(seconds){
    # floor() or round() would work as well
    seconds = ceiling(seconds)
    days = seconds %/% (60*60*24)
    seconds = seconds %% (60*60*24)
    hours = seconds %/% (60*60)
    seconds = seconds %% (60*60)
    minutes = seconds %/% (60)
    seconds = seconds %% (60)
    return (paste0(days,"d ",hours,"h ",minutes,"m ",seconds, "s"))
    }

time_elapse = function(iter, nmax, every = 100){
    if (iter == 1){
        begin_time <<- as.numeric(Sys.time())
    } else {
        curr_time <<- as.numeric(Sys.time()) - begin_time
        spi = curr_time / iter
        if (floor(iter/every) == iter/every){
            cat(paste0("\rTime remaining: ", nice.time(spi*(nmax-iter))," "))
            }
        }
    }


nproc = 3
library(doMC)
registerDoMC(nproc)


set.seed(1)
X = matrix(rnorm(1000^2), 1000, 1000)
mat = t(X) %*% X

system.time(solve(mat))

k = 1000

### without printing the time
actual_start1 = as.numeric(Sys.time())
x = foreach(i = 1:k) %dopar% {
    solve(mat)
    }
actual_end1 = as.numeric(Sys.time())
nice.time(actual_end1 - actual_start1) # total time elapsed


### just catting the iteration number each time
actual_start2 = as.numeric(Sys.time())
x = foreach(i = 1:k) %dopar% {
    solve(mat)
    cat(i, "\n")
    }
actual_end2 = as.numeric(Sys.time())
nice.time(actual_end2 - actual_start2) # total time elapsed



actual_start3 = as.numeric(Sys.time())
x = foreach(i = 1:k) %dopar% {
    solve(mat)
    time_elapse(i, k, every = 1)
    }
actual_end3 = as.numeric(Sys.time())
nice.time(actual_end3 - actual_start3) # total time elapsed



actual_start4 = as.numeric(Sys.time())
for (i in 1:k){
    solve(mat)
    time_elapse(i, k, every = 1)
    }
actual_end4 = as.numeric(Sys.time())
nice.time(actual_end4 - actual_start4) # total time elapsed

