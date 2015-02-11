library(foreign)
dat_sf = data.frame(read.spss("~/files/data/637/mormon_sf_sample.SAV"))
dat_slc = read.csv("~/files/data/637/mormon_slc_sample.csv")

logit = function(x)
    log(x/(1-x))
logistic = function(x)
    1/(1+exp(-x))

### data cleaning
# initialize
x = data.frame(matrix(0, nrow(dat_slc) + nrow(dat_sf), 12))
colnames(x) = c("CITY", "LDS", "PRNTACTV", "SACRMTG", "AGE",
    "INCOME", "EDUC", "PRVPRAYR", "READSCRP", "MARITAL",
    "SEX", "FRIEND")
# 256, 28, 251, 252, 271, 36, 33, 270, 269, 3, 4, 6

# city indicator (0 for SF, 1 for SLC)
x$CITY = c(rep(0, nrow(dat_sf)), rep(1, nrow(dat_slc)))

# born in the church, or convert (7)
temp_sf = ifelse(dat_sf[,7] == "Born in the LDS Church, but no longer affiliated" |
    dat_sf[,7] == "Life-long member and still affiliated", "birth", "convert")
temp_slc = ifelse(dat_slc[,7] == 1 | dat_slc[,7] == 2, "birth", "convert")
x$LDS = c(temp_sf, temp_slc)

# parents activity level (256)
x$PRNTACTV = c(dat_sf[,256], dat_slc[,256])

# sacrament meeting attendance (28)
x$SACRMTG = c(dat_sf[,28], dat_slc[,28])

# age, continu-ized (251)
x$AGE = c(dat_sf[,251], dat_slc[,251])
for (i in 1:length(x$AGE)){
    if (is.na(x$AGE[i]))
        next
    if (x$AGE[i] == 1) x$AGE[i] = 20
    if (x$AGE[i] == 2) x$AGE[i] = 25
    if (x$AGE[i] == 3) x$AGE[i] = 30
    if (x$AGE[i] == 4) x$AGE[i] = 35
    if (x$AGE[i] == 5) x$AGE[i] = 40
    if (x$AGE[i] == 6) x$AGE[i] = 45
    if (x$AGE[i] == 7) x$AGE[i] = 50
    if (x$AGE[i] == 8) x$AGE[i] = 55
    if (x$AGE[i] == 9) x$AGE[i] = 60
    if (x$AGE[i] == 10) x$AGE[i] = 65
    }

# income continu-ized (252) (smaller)
x$INCOME = c(dat_sf[,252], dat_slc[,252])
for (i in 1:length(x$INCOME)){
    if (is.na(x$INCOME[i]))
        next
    if (x$INCOME[i] == 1) x$INCOME[i] = 4000
    if (x$INCOME[i] == 2) x$INCOME[i] = 5000
    if (x$INCOME[i] == 3) x$INCOME[i] = 6000
    if (x$INCOME[i] == 4) x$INCOME[i] = 7000
    if (x$INCOME[i] == 5) x$INCOME[i] = 8000
    if (x$INCOME[i] == 6) x$INCOME[i] = 9000
    if (x$INCOME[i] == 7) x$INCOME[i] = 10000
    if (x$INCOME[i] == 8) x$INCOME[i] = 11000
    if (x$INCOME[i] == 9) x$INCOME[i] = 15000
    if (x$INCOME[i] == 10) x$INCOME[i] = 20000
    if (x$INCOME[i] == 11) x$INCOME[i] = 35000
    if (x$INCOME[i] == 12) x$INCOME[i] = 50000
    }

# education level (271)
x$EDUC = c(dat_sf[,271], dat_slc[,271])

# private prayer (36)
x$PRVPRAYR = c(dat_sf[,36], dat_slc[,36])

# scripture reading (33)
x$READSCRP = c(dat_sf[,33], dat_slc[,33])

# marital status (270)
x$MARITAL = c(dat_sf[,270], dat_slc[,270])

# sex (269)
x$SEX = 2 - c(dat_sf[,269], dat_slc[,269])

# # priesthood office (3)
# x$PSTH = c(dat_sf[,3], dat_slc[,3])
# x$PSTH[x$SEX == 0] = 0  # set priesthood to 0 if female
#                         # may need to deal with this
#                         # just ignore priesthood office?

# lds friends (4)
x$FRIEND = c(c(dat_sf[,4]) - 1, dat_slc[,4])

# # mission (6)
# x$MISSION = c(dat_sf[,6], dat_slc[,6]) - 1

# response
temp_sf = ifelse(dat_sf[,7] == "Born in the LDS Church, but no longer affiliated" |
    dat_sf[,7] == "Once converted, but no longer affiliated", 0, 1)
temp_slc = ifelse(dat_slc[,7] == 1 | dat_slc[,7] == 4, 0, 1)
y = c(temp_sf, temp_slc)

# mising values
miss.y = which(is.na(y))
miss.x = NULL
for (i in 1:ncol(x))
    miss.x = unique(c(miss.x, which(is.na(x[,i]))))
# every row that has at least one missing value
miss.all = c(miss.y, miss.x)

# counts for responses after removing missing values
table(y[-miss.all])

# check to make sure they are removed
sum(is.na(y[-miss.all]))
sum(is.na(x[-miss.all,]))

y = y[-miss.all]
x = x[-miss.all,]

# make all the x's factors (except age and income)
x[,names(x)] = lapply(x[,names(x)], factor)
x$AGE = c(x$AGE)
x$INCOME = c(x$INCOME)
x$FRIEND = c(x$FRIEND)

### model fitting
mod = glm(y ~ CITY + PRVPRAYR + LDS + SACRMTG + FRIEND +
    PRVPRAYR*CITY, data = x, family = binomial(link = "logit"))
summary(mod)

### confidence intervals on p's with various sets of x's
x.conf = x[1:6, -c(3, 5, 6, 7, 9, 10, 11)]
x.conf[,2] = "birth"
x.conf[,1] = c(0,0,0,1,1,1)
x.conf[,3] = c(1,4,6,1,4,6)
x.conf[,4] = c(1,2,6,1,2,6)
x.conf[,5] = c(5,2,0,5,2,0)
rownames(x.conf) = c("SF_High", "SF_Medium", "SF_Low",
    "SLC_High", "SLC_Medium", "SLC_Low")
x.conf[,names(x.conf)] = lapply(x.conf[,names(x.conf)], factor)
x.conf$FRIEND = c(x.conf$FRIEND)
x.conf[,5] = c(5,2,0,5,2,0)

p = predict(mod, se.fit = TRUE, newdata = x.conf)

conf.int = matrix(0, 6, 2)
for (i in 1:nrow(x.conf))
    conf.int[i,] = logistic(p$fit[i] + c(-1, 1)*qnorm(0.975)*p$se.fit[i])

rownames(conf.int) = rownames(x.conf)

library(xtable)
write.table(print(xtable(conf.int, digits = 3)), "./conf_int.txt", col.names = FALSE, quote = FALSE)
write.table(print(xtable(summary(mod))), "./coef_est.txt", col.names = FALSE, quote = FALSE)
