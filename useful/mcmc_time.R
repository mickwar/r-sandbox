### Functions for calculating elapsed time of executed code and
### estimated remaining time until completion of the code.
### Made specifically with MCMC chains in mind that would be
### run in the background for several days. A .txt file is created
### with mcmc_time() that displays the information as well as
### outputting the draws periodically.

# get current list of objects in the console
TEMP_LIST_O = ls()

# trailing()
# Truncate trailing digits
#
# Params: x      - numeric, to be truncated
#         digits - numeric, the number of decimal places that
#                  x will be trunacted to
trailing = function(x, digits = 4)
    formatC(x, digits=digits, format="f")

# nice.time()
# Returns a string in nice time format given the number of seconds
#
# Params: seconds - numeric, converted to days/hours/minutes/seconds
nice.time = function(seconds){
    # floor() or round() would work as well
    seconds = ceiling(seconds)
    days = seconds %/% (60*60*24)
    seconds = seconds %% (60*60*24)
    hours = seconds %/% (60*60)
    seconds = seconds %% (60*60)
    minutes = seconds %/% (60)
    seconds = seconds %% (60)
    out = ""
    if (days > 0)
        out = paste0(out, days, "d ", hours, "h ", minutes, "m ", seconds, "s")
    if (days == 0 && hours > 0)
        out = paste0(out, hours, "h ", minutes, "m ", seconds, "s")
    if (days == 0 && hours == 0 && minutes > 0)
        out = paste0(out, minutes, "m ", seconds, "s")
    if (days == 0 && hours == 0 && minutes == 0)
        out = paste0(out, seconds, "s")
    return (out)
    }

# mcmc_time()
# Estimate remaining time for an MCMC and output posterior draws
# and other objects periodically.
#
# Params: do     - logical, should the function be run? 
#         iter   - numeric, value of 1 initializes the process, any
#                  other value increments the process. M-H chains
#                  often began at 2 (iter = 1 are the starting values),
#                  and this function assumes that.
#         every  - numeric, greater than 1. How often should the results
#                  be written to external text files. It should be the case
#                  that (nmcmc + nburn) %% every == 0. For computationally
#                  intensive loops, it's a good idea to output the results
#                  periodically. If the loop takes only a few minutes, then set
#                  every = nburn + nmcmc (output only at the end)
#         params - (nburn + nmcmc) by (number of parameters) matrix of
#                  posterior draws
#         accept - similar to params, but of the acceptances (binary matrix)
#           sigs - candidate sigmas for proposal distributions
#          nburn - number of burnin draws
#          nmcmc - number of mcmc iterations (does not include burnin,
#                  i.e. total iterations = nburn + nmcmc
#
# This function uses variables in global scope and so should be reserved.
# These are:
#     mcmc_begin_time
#     mcmc_curr_time
#
# Requires function nice.time()
# (note: The periodic outputting of results should be worked on. Not every
# sets up their MCMC the same way as I do).
mcmc_time = function(do = TRUE, iter, every = 100, params, accept,
    sigs, nburn, nmcmc, dir = ".", prefix = "mcmc_"){
    if (!do)
        return (0)
    if (iter == 1){
        mcmc_begin_time <<- as.numeric(Sys.time())
    } else {
        mcmc_curr_time <<- as.numeric(Sys.time()) - mcmc_begin_time
        spi = mcmc_curr_time / (iter - 1) # seconds per iteration
        write.table(paste0("Iteration: ", iter, "\n",
            "  Burn-in: ", nburn, "\n",
            "    NMCMC: ", nmcmc, "\n\n",
            "  Elapsed: ", nice.time(mcmc_curr_time), "\n",
            "Remaining: ", nice.time(spi*(nmcmc+nburn-iter)), "\n\n",
            " Complete: ", trailing(100*iter/(nburn+nmcmc)),"%"),
            paste0(dir, "/", prefix, "status.txt"),
            col.names=F, row.names=F, quote=F)
        if (floor(iter/every) == iter/every){
            write.table(signif(params[(1+every*(ceiling(iter/every)
                -1)):iter,]), paste0(dir, "/", prefix, "params.txt"),
                row.names=F, col.names=F, quote=F, append=T)
            write.table(signif(accept[(1+every*(ceiling(iter/every)
                -1)):iter,]), paste0(dir, "/", prefix, "accept.txt"),
                row.names=F, col.names=F, quote=F, append=T)
            # print the latest parameter values and candidate sigmas
            # until burn-in is complete, so the values may be used
            # as starting point for future runs
            if (iter <= nburn + 1){
                write.table(params[iter,],
                    paste0(dir, "/", prefix, "init_params.txt"),
                    row.names=F, col.names=F, quote=F)
                write.table(sigs,
                    paste0(dir, "/", prefix, "cand_sigmas.txt"),
                    row.names=F, col.names=F, quote=F)
                }
            }
        }
    }


###
# compare TEMP_LIST with ls() which now has additional functions
TEMP_LIST_N = ls()
TEMP_LIST_N = TEMP_LIST_N[-which(TEMP_LIST_N == "TEMP_LIST_O")] # remove TEMP_LIST_O
TEMP_LIST_N = TEMP_LIST_N[!(TEMP_LIST_N %in% TEMP_LIST_O)] # get only added functions

# print the sourced functions
if (length(TEMP_LIST_N) > 0){
    cat("Sourced functions:\n    ")
    cat(TEMP_LIST_N, sep="\n    ")
    }

rm(TEMP_LIST_O, TEMP_LIST_N) # remove temporary variables
