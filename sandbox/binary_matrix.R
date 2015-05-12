### Convert a factor into a binary matrix with the number of
### columns equal to the number of factor levels. The columns
### in the matrix correspond to the order of the factor levels

### For example, the level of the first element in the factor
### is indicated by a 1 in the first column of the matrix, and
### all subsequent elements having that same factor level will
### have a 1 in the first column and 0 in the others.

binary.matrix = function(x){
    level = unique(x)
    n = length(x)
    k = length(level)
    out = matrix(0, n, k)
    for (i in 1:k)
        out[which(x == level[i]), i] = 1
    return (out)
    }

### example usage
set.seed(4)
n = 100
k = 6
af1 = factor(sample(k, n, replace=TRUE))
ex = binary.matrix(af1)

head(af1)
head(ex)

# full
table(af1)
ex
