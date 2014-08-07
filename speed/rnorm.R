set.seed(1)

n = 100
a = double(n)

for (i in 1:n){
    a[i] = system.time(rnorm(10000000))[3]
    }

write.table(a, "./time_rnorm.txt", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
