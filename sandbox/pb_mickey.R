### This code provides functions for implementing
### a simple progress bar in the console

# function to force a numeric to go out to
# "digits" decimal places
trailing = function(x, digits=2)
    formatC(x, digits=digits, format="f")

# begin the progress bar
pb.Init = function()
    cat("\r    |",rep(" ", 50),"| 0.00%", sep="")

# assumes iteration takes on values 1:max
# i.e. adjustments need to be made for other sequences,
# especially for non-sequential (e.g. i in c(1,3,4,7))
pb.Update = function(iteration, max){
    fill = floor(50*iteration/max)
    cat("\r    |",rep("=", fill),rep(" ", 50-fill),"| ",
        trailing(100*iteration/max, 2),"%", sep="")
    if (iteration == max)
        cat("\n")
    }

# example
# start = 10
# niter = 70
# for (i in start:niter){
#     if (i == start){ # beginning
#         pb.Init()
#         Sys.sleep(1.00)
#         }
#     # do functions
#     Sys.sleep(0.10)
#     pb.Update(i-start+1, niter-start+1)
#     }
