# balances

dat = read.table("./electric_bill.txt", header=TRUE, row.names=1)

# place NA for people who don't need to pay for that month
# calculate everyone's portion of the bill, with those not
# making the actual payment (not column 2) rounded to the
# nearest dollar
update = function(dat){
    total_payee = 0
    for (i in 1:nrow(dat)){
        if (any(dat[i,] == 0, na.rm = TRUE)){
            ind = which(!is.na(dat[1:(i-1),2]))
            if (i == 1)
                ind = logical(0)
            should = sum(dat[ind, 1] / apply(dat, 1, function(x)
                sum(!is.na(x)) - 1)[ind])
            if (any(is.na(dat[i,]))){
                j = (1:ncol(dat))[-c(1,2,which(is.na(dat[i,])))]
            } else {
                j = 3:ncol(dat)
                }
            if (total_payee >= should){
                dat[i,j] = ceiling(dat[i,1] / (length(j) + 1))
            } else {
                dat[i,j] = floor(dat[i,1] / (length(j) + 1))
                }
            dat[i,2] = dat[i,1] - sum(dat[i,3:ncol(dat)], na.rm = TRUE)
            }
        total_payee = total_payee + dat[i,2]
        }
    return (dat)
    }

# this needs work, it's to implement specific payments made
# by others, if they over and underpaid in a month for instance
# not complete, don't use it
update2 = function(dat){
    total_payee = 0
    for (i in 1:nrow(dat)){
        k = which(dat[i,] == 0)
        if (length(k) > 0){
            ind = which(!is.na(dat[1:(i-1),2]))
            if (i == 1)
                ind = logical(0)
            should = sum(dat[ind, 1] / apply(dat, 1, function(x)
                sum(!is.na(x)) - 1)[ind])
            k = k[-which(k == 2)]
            if (length(k) > 0){
                if (total_payee >= should){
                    dat[i,j] = ceiling(dat[i,1] / (length(j) + 1))
                } else {
                    dat[i,j] = floor(dat[i,1] / (length(j) + 1))
                    }
                }
            dat[i,2] = dat[i,1] - sum(dat[i,3:ncol(dat)], na.rm = TRUE)
            }
        total_payee = total_payee + dat[i,2]
        }
    return (dat)
    }

(newdat= update2(dat))
(newdat= update(dat))

# save to the text file

# check to see we're all making the right payments the value
# from apply() should be close to each of the sums
 apply(newdat, 2, sum, na.rm=TRUE)
 ind = which(!is.na(newdat[,2]))
 sum(dat[ind, 1]/apply(newdat, 1, function(x) sum(!is.na(x)) - 1)[ind])
 ind = which(!is.na(newdat[,3]))
 sum(dat[ind, 1]/apply(newdat, 1, function(x) sum(!is.na(x)) - 1)[ind])
#ind = which(!is.na(newdat[,4]))
#sum(dat[ind, 1]/apply(newdat, 1, function(x) sum(!is.na(x)) - 1)[ind])
#ind = which(!is.na(newdat[,5]))
#sum(dat[ind, 1]/apply(newdat, 1, function(x) sum(!is.na(x)) - 1)[ind])
