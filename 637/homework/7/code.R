b1 = c(0.883, 0.419, 0.342)
b2 = c(-0.758, 0.105, 0.271)

x = c(1, 1, 1)

(exp(t(x) %*% b1))/(1 + exp(t(x) %*% b1) + exp(t(x) %*% b2))

dat = c("democrat", 132, "male", "white",
    "republican", 176, "male", "white",
    "independent", 127, "male", "white",
    "democrat", 42, "male", "black",
    "republican", 6, "male", "black",
    "independent", 12, "male", "black",
    "democrat", 172, "female", "white",
    "republican", 129, "female", "white",
    "independent", 130, "female", "white",
    "democrat", 56, "female", "black",
    "republican", 4, "female", "black",
    "independent", 15, "female", "black")

dat = (matrix(dat, 12, 4, byrow = TRUE))


multinom(factor(dat[,1]) ~ factor(dat[,3]) + factor(dat[,4]), weights = as.numeric(dat[,2]))
