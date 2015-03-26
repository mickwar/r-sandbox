library(aod)

head(dja)

mod = glm(cbind(y, n-y) ~ group, data = dja, family = quasibinomial)
summary(mod)

n = nrow(dja)
p = 2 # intercept and group effect

(phi = summary(mod)$dispersion)

(n-p)* phi
qchisq(0.95, n-p) 
