### Density function and obtaining a random draws from the Skew Normal distribution
### Use the 'sn' package for more Skew Normal functions

# get current list of objects in the console
TEMP_LIST_O = ls()

### dskewnorm()
# Calculate the density function of a Skew Normal random variable
#
# Params: x     - numeric, real (support)
#         ksi   - numeric, real (location)
#         omega - numeric, positive (scale)
#         alpha - numeric, real (shape)
#         log   - logical, should the log be calculated? (note that
#                 this takes slightly longer to compute)
dskewnorm = function(x, ksi, omega, alpha, log = FALSE){
    if (log){
        log(2) - log(omega) + dnorm((x-ksi)/omega, 0, 1, log = TRUE) + 
            pnorm(alpha * (x - ksi)/omega, 0, 1, log.p = TRUE)
    } else {
        2/omega * dnorm((x-ksi)/omega, 0, 1) * pnorm(alpha*(x-ksi)/omega, 0, 1)
        }
    }

### rskewnorm()
# 
#
# Params: n     - numeric, positive integer
#         ksi   - numeric, real (location)
#         omega - numeric, positive (scale)
#         alpha - numeric, real (shape)
rskewnorm = function(n, ksi, omega, alpha){
    delta = alpha / sqrt(1 + alpha^2)

    # draw 2*n standard normals
    u_0 = rnorm(n, 0, 1)
    v = rnorm(n, 0, 1)

    # impose skewness
    u_1 = delta*u_0 + sqrt(1-delta^2)*v

    z = ifelse(u_0 >= 0, u_1, -u_1)

    # transform location and scale
    z = z*omega + ksi

    return (z)
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
