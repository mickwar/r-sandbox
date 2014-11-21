### read in the data
dat = read.table("~/files/R/666/data/collins.txt", header = TRUE)
source("functions.R")

### remove unnecessary rows
dat = dat[,-which(apply(dat, 2, function(x) length(unique(x))) == nrow(dat))]

### Step 1
# use only first 18 variables
x = scale(dat[,1:18])

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

pdf("figs/dendrogram.pdf", width = 18, height = 9)
s4.5 = mw.tree(x, 5, scores, "ward.D")
dev.off()
1
# don't worry about the error rates so much, but ask, "are the clusters
# roughly in the same location in both plots?"
validate(x, 3, scores, seed = 11)
validate(x, 4, scores, seed = 3)
validate(x, 5, scores, seed = 3)
validate(x, 6, scores, seed = 53)
validate(x, 7, scores, seed = 13)

# use k-nearest neighbors, with k-fold cross-validation, to estimate
# error rates for newly acquired data
(errors = c(knearest(x,s4.3, kfold = 10, floor(sqrt(min(s4.3$counts)))),
    knearest(x,s4.4, kfold = 10, floor(sqrt(min(s4.4$counts)))),
    knearest(x,s4.5, kfold = 10, floor(sqrt(min(s4.5$counts)))),
    knearest(x,s4.6, kfold = 10, floor(sqrt(min(s4.6$counts)))),
    knearest(x,s4.7, kfold = 10, floor(sqrt(min(s4.7$counts))))))
round(errors, 3)

# using super genre to see how well the genres can group writing styles
genre = dat$Genre
super.genre = double(length(genre))
super.genre[genre == 1 | genre == 2 | genre == 3] = 1
super.genre[genre == 4 | genre == 5 | genre == 6] = 2
super.genre[genre == 7] = 3
super.genre[genre == 8 | genre == 9] = 4
super.genre[genre == 10 | genre == 11 | genre == 12 | genre == 13 |
    genre == 14 | genre == 15] = 5
G = NULL
G$cutree = super.genre
G$counts = table(super.genre)

knearest(x, G, 10, floor(sqrt(min(G$counts))))
# An error rate of about 0.36, compared with the 5-cluster error
# rate of about 0.10. So writing styles aren't compared best using
# the super genre.
