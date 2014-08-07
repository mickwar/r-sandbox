set.seed(1)

n = 100
a = double(n)

for (i in 1:n){
    a[i] = system.time({x = matrix(runif(1000000), 1000, 1000);
        solve(t(x) %*% x)})[3]
    }

write.table(a, "./time_solve.txt", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
