### How random is this string of characters?

### In this code, I'm going to attempt to look at the
### randomness of a string of characters (a-z and 0-9)
### and try to determine if they are sufficiently random.


### So far (12 Jan 2013) the function only takes into account
### the actual characters used.

### To do:  Test randomness by the order of the characters.
###        Look at the grouping of characters in the string.
###        Fix "need finite 'ylim' values" plot error.
###            (This occurs when there is only one unique
###             character in the string being tested.)
###        Take several samples of various size "n" to get
###            data for passing proportions.

random.string=function(n,replace=TRUE){
    order="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890!@#$%^&*()`~-_=+[];,'./{}|:<>?"
    numbers=sample(nchar(order),n,replace=replace)
    string=NULL
    for (i in 1:n){
        string=paste(string,sep="",substring(order,numbers[i],numbers[i]))
        }
    return(string)
    }


random.r.string=function(n,replace=TRUE){
    order="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
    numbers=sample(nchar(order),n,replace=replace)
    string=NULL
    for (i in 1:n){
        string=paste(string,sep="",substring(order,numbers[i],numbers[i]))
        }
    return(string)
    }

random.test=function(pw){
    if (nchar(pw)<2)
        stop("Use a string of at least two characters.")
    order=random.string(92,replace=FALSE)
    L=nchar(pw)
    value=NULL
    for (i in 1:L){
        failcount=0
        for (j in 1:nchar(order)){
            if (substring(pw,i,i) != substring(order,j,j))
                failcount=failcount+1
            else (value[i]=j)
            }
        if (failcount==nchar(order))
            stop("Invalid character string.  Use the alphabet and numbers 0-9 only.")
        }
    exp=mean(seq(1,nchar(order),1))
    stdev=sd(seq(1,nchar(order),1))
    realmu=mean(value)
    realsd=sd(value)
    par(mfrow=c(1,2))
    curve(dnorm(x,exp,stdev),from=min(exp-3*stdev,realmu-3*realsd),
        to=max(exp+3*stdev,realmu+3*realsd),col='black',lwd=2,
        main="Expected (Black), vs. Real (Red)",
        ylim=c(0,max(dnorm(exp,exp,stdev),dnorm(realmu,realmu,realsd))))
    curve(dnorm(x,realmu,realsd),add=T,col='red',lwd=2)
    stdev=stdev/sqrt(L)
    curve(dnorm(x,exp,stdev),from=0,to=nchar(order)+1,lwd=2,
        main="Character Distribution Test")
    for (i in c(1,-1,2,-2,3,-3)){
        lines(c(exp+i*stdev,exp+i*stdev),c(0,dnorm(exp+i*stdev,exp,stdev)),
            col=(abs(i)-1+(3/2)*(abs(i)-2)*(abs(i)-3)))
        }
    abline(v=realmu,col='blue',lwd=2)
    perc=pnorm(exp-abs(exp-realmu),exp,stdev)
    range=pnorm(exp-stdev,exp,stdev)
    if (perc>range)
        ctest="PASSED"
    else (ctest="FAILED")
    print("Character Distribution Test")
    print(paste("     Strength:",perc))
    print(paste("     Desired strength should be >",range))
    print(paste("Result:",ctest))
    }
