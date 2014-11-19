### read in the data
dat = read.table("~/files/R/666/data/collins.txt", header = TRUE)
source("functions.R")

### remove unnecessary rows
dat = dat[,-which(apply(dat, 2, function(x) length(unique(x))) == nrow(dat))]

### Step 1
# use only first 18 variables
x = dat[,1:18]

# get principal component scores
scores = get.pc.scores(x, k = 1:3)
scores$p                # proportion of variance
scores = scores$scores  # just retain the scores for each observation

# crap
#s1.3 = mw.tree(x, 3, scores, "single"); s1.3$counts
#s1.4 = mw.tree(x, 4, scores, "single"); s1.4$counts
#s1.5 = mw.tree(x, 5, scores, "single"); s1.5$counts
#s1.6 = mw.tree(x, 6, scores, "single"); s1.6$counts
#s1.7 = mw.tree(x, 7, scores, "single"); s1.7$counts

#s2.3 = mw.tree(x, 3, scores, "complete"); s2.3$counts
#s2.4 = mw.tree(x, 4, scores, "complete"); s2.4$counts
#s2.5 = mw.tree(x, 5, scores, "complete"); s2.5$counts
#s2.6 = mw.tree(x, 6, scores, "complete"); s2.6$counts
#s2.7 = mw.tree(x, 7, scores, "complete"); s2.7$counts

#s3.3 = mw.tree(x, 3, scores, "average"); s3.3$counts
#s3.4 = mw.tree(x, 4, scores, "average"); s3.4$counts
#s3.5 = mw.tree(x, 5, scores, "average"); s3.5$counts
#s3.6 = mw.tree(x, 6, scores, "average"); s3.6$counts
#s3.7 = mw.tree(x, 7, scores, "average"); s3.7$counts

s4.3 = mw.tree(x, 3, scores, "ward.D"); s4.3$counts
s4.4 = mw.tree(x, 4, scores, "ward.D"); s4.4$counts
s4.5 = mw.tree(x, 5, scores, "ward.D"); s4.5$counts
s4.6 = mw.tree(x, 6, scores, "ward.D"); s4.6$counts
s4.7 = mw.tree(x, 7, scores, "ward.D"); s4.7$counts

s5.3 = mw.tree(x, 3, scores, "ward.D2"); s5.3$counts
s5.4 = mw.tree(x, 4, scores, "ward.D2"); s5.4$counts
s5.5 = mw.tree(x, 5, scores, "ward.D2"); s5.5$counts
s5.6 = mw.tree(x, 6, scores, "ward.D2"); s5.6$counts
s5.7 = mw.tree(x, 7, scores, "ward.D2"); s5.7$counts

#s6.3 = mw.tree(x, 3, scores, "mcquitty"); s6.3$counts
#s6.4 = mw.tree(x, 4, scores, "mcquitty"); s6.4$counts
#s6.5 = mw.tree(x, 5, scores, "mcquitty"); s6.5$counts
#s6.6 = mw.tree(x, 6, scores, "mcquitty"); s6.6$counts
#s6.7 = mw.tree(x, 7, scores, "mcquitty"); s6.7$counts
#
#s7.3 = mw.tree(x, 3, scores, "median"); s7.3$counts
#s7.4 = mw.tree(x, 4, scores, "median"); s7.4$counts
#s7.5 = mw.tree(x, 5, scores, "median"); s7.5$counts
#s7.6 = mw.tree(x, 6, scores, "median"); s7.6$counts
#s7.7 = mw.tree(x, 7, scores, "median"); s7.7$counts
#
#s8.3 = mw.tree(x, 3, scores, "centroid"); s8.3$counts
#s8.4 = mw.tree(x, 4, scores, "centroid"); s8.4$counts
#s8.5 = mw.tree(x, 5, scores, "centroid"); s8.5$counts
#s8.6 = mw.tree(x, 6, scores, "centroid"); s8.6$counts
#s8.7 = mw.tree(x, 7, scores, "centroid"); s8.7$counts

classify(x, s5.7, "linear", k = 2)
classify(x, s5.7, "linear", k = 3)
classify(x, s5.7, "linear", k = 5)
classify(x, s5.7, "linear", k = 10)

at = 5:30
error.rates = double(length(at))
for (i in 1:length(at))
    error.rates[i] = knearest(x, s5.6, kfold = 10, knear = at[i])
plot(at, error.rates, type='l')

k6.k = at[which.min(error.rates)]
k6.e = min(error.rates)

c(k3.k, k3.e)
c(k4.k, k4.e)
c(k5.k, k5.e)
c(k6.k, k6.e)
c(k7.k, k7.e)

library(rgl)
plot3d(scores, col=s4.7$cutree, size=5)
