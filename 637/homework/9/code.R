y = c(76, 160, 6, 25, 114, 181, 11, 48)
x = matrix(c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0), 8, 3, byrow = TRUE)
x[,1] = ifelse(x[,1] == 1, "male", "female")
x[,2] = ifelse(x[,2] == 1, "support", "oppose")
x[,3] = ifelse(x[,3] == 1, "support", "oppose")
x = data.frame(x)
names(x) = c("Gender", "Information", "Health")

### four models
mod.gh.gi = glm(y ~ Gender*Health + Gender*Information, data = x,
    family = poisson(link = "log"))
mod.gh.hi = glm(y ~ Gender*Health + Health*Information, data = x,
    family = poisson(link = "log"))
mod.gi.hi = glm(y ~ Gender*Information + Health*Information, data = x,
    family = poisson(link = "log"))
mod.gh.gi.hi = glm(y ~ Gender*Health + Gender*Information + Health*Information, data = x,
    family = poisson(link = "log"))


summary(mod.gh.gi)
summary(mod.gh.hi)
summary(mod.gi.hi)
summary(mod.gh.gi.hi)


mod.full = glm(y ~ . + .^2 + .^3, data = x, family = poisson(link = "log"))
summary(mod.full)

