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

# construct the 5 best clusterings
s4.3 = mw.tree(x, 3, scores, "ward.D"); s4.3$counts
s4.4 = mw.tree(x, 4, scores, "ward.D"); s4.4$counts
s4.5 = mw.tree(x, 5, scores, "ward.D"); s4.5$counts
s4.6 = mw.tree(x, 6, scores, "ward.D"); s4.6$counts
s4.7 = mw.tree(x, 7, scores, "ward.D"); s4.7$counts


errors = double(5)
errors[1] = knearest(x, s4.3, kfold = 10, floor(sqrt(min(s4.3$counts))))
errors[2] = knearest(x, s4.4, kfold = 10, floor(sqrt(min(s4.4$counts))))
errors[3] = knearest(x, s4.5, kfold = 10, floor(sqrt(min(s4.5$counts))))
errors[4] = knearest(x, s4.6, kfold = 10, floor(sqrt(min(s4.6$counts))))
errors[5] = knearest(x, s4.7, kfold = 10, floor(sqrt(min(s4.7$counts))))

errors
