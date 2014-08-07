set.seed(1)

n = 100
a = double(n)

for (i in 1:n){
    x = double(1000000)
    a[i] = system.time(for (j in 1:1000000){ x[j] = 1.0 })[3]
    }

write.table(a, "./time_assign.txt", quote = FALSE,
    row.names = FALSE, col.names = FALSE)

# is this a good way to time assignment?
