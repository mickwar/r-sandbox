library(aod)

head(dja)

### quasi binomial
mod = glm(cbind(y, n-y) ~ group, data = dja, family = quasibinomial)
summary(mod)

n = nrow(dja)
p = 2 # intercept and group effect

(phi = summary(mod)$dispersion)

(n-p)*phi
qchisq(0.95, n-p) 

par(mfrow=c(2,2))
plot(mod, ask = FALSE)

### beta-binomial
mod = betabin(cbind(y, n-y) ~ group, ~1, data = dja, link = "logit")
summary(mod)

par(mfrow=c(2,2))
plot(summary(mod), ask = FALSE)


### quasi-poisson
mod = glm(y ~ group, offset = log(trisk), data = dja,
    family = quasipoisson)
summary(mod)

n = nrow(dja)
p = 2 # intercept and group effect

(phi = summary(mod)$dispersion)

(n-p)*phi
qchisq(0.95, n-p) 

par(mfrow=c(2,2))
plot(mod, ask = FALSE)

### negative binomial
library(MASS)
mod = glm(y ~ group, offset = log(trisk), data = dja,
    family = negative.binomial(theta = 1))
summary(mod)

n = nrow(dja)
p = 2 # intercept and group effect

(phi = summary(mod)$dispersion)

(n-p)*phi
qchisq(0.95, n-p) 

par(mfrow=c(2,2))
plot(mod, ask = FALSE)
