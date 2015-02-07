library(foreign)
dat_sf = data.frame(read.spss("~/files/data/637/mormon_sf_sample.SAV"))
dat_slc = read.csv("~/files/data/637/mormon_slc_sample.csv")

# data cleaning
x = matrix(0, nrow(dat_slc) + nrow(dat_sf), 1)
temp_sf = ifelse(dat_sf[,7] == "Born in the LDS Church, but no longer affiliated" |
    dat_sf[,7] == "Once converted, but no longer affiliated", 0, 1)
dat_slc[,7] = ifelse(dat_slc[,7] == 1 | dat_slc[,7] == 4, 0, 1)
dat_sf[,7] = ifelse(dat_sf[,7] == "Born in the LDS Church, but no longer affiliated" |
    dat_sf[,7] == "Once converted, but no longer affiliated", 0, 1)

# 7
256, 28, 251, 1, 252, 271, 36, 33, 270, 269, 3, 7
