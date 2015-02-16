library(RCurl)

dat = read.table(text = getURL("https://raw.githubusercontent.com/luiarthur/Fall2014/master/Stat637/5/sore.txt"),
    header = TRUE)

dat = dat[,-1]

dat = dat[order(dat[,3]),]
dat = dat[order(dat[,2]),]
dat = dat[order(dat[,1]),]
