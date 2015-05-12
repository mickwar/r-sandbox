set.seed(1)
auror = read.table("~/files/data/auror.txt")
auror = auror[-which(auror[,3] == 0),]
n = nrow(auror)

# earnwke (weekly earnings in galleons)
# senior (senior if an auror prior to the second fall of voldemort,
#         so a junior if the person became an auror after harry potter
#         destroyed voldemort)
# union (union if member of a death eaters union)

pairs(auror)

# mod = lm(earnwke ~ age + race + senior + married + union +
#     female + hogwarts + hhouse, data = auror)

summary(mod)
train = sort(sample(n, floor(n*0.5)))
test = (1:n)[-train]

null.mod = glm(earnwke ~ 1, data = auror, subset = train, family = Gamma(link = "log"))
full.mod = glm(earnwke ~ ., data = auror, subset = train, family = Gamma(link = "log"))

mod = step(null.mod, scope = list("lower"=null.mod, "upper"=full.mod),
    direction = "both", k = log(floor(n*0.5)))
summary(mod)

pred = predict.glm(mod, newdata = auror[test,], se.fit = TRUE)
pred = predict.glm(mod, newdata = auror[test,], type = "response")
mean((pred - auror[test, 3])^2)


pred = exp(cbind(pred$fit - qnorm(0.975) * pred$se.fit,
    pred$fit, pred$fit + qnorm(0.975) * pred$se.fit))

# mse
mean((pred[,2] - auror[test, 3])^2)

# coverage
mean(auror[test, 3] >= pred[,1] & auror[test, 3] <= pred[,3])
