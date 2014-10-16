### 3d

to.f = function(lam, nu.e, nu.h)
    (1-sqrt(lam))/sqrt(lam) * (nu.e - 1)/nu.h

(lam1 = 0.000645 / 0.00527159)
(lam2 = 0.0653 / 0.06994306)
(lam3 = 0.1379 / 0.21546479)

q = 2
p = 2
n = 60

# sowdate
k = 4
nu.h = k - 1
nu.e = k*(n-1)

(F.stat = to.f(lam1, nu.e, nu.h))
pf(F.stat, 2*nu.h, 2*(nu.e - 1), lower.tail = FALSE)

# variety
k = 3
nu.h = k - 1
nu.e = k*(n-1)

(F.stat = to.f(lam2, nu.e, nu.h))
pf(F.stat, 2*nu.h, 2*(nu.e - 1), lower.tail = FALSE)

# interaction
nu.h = (4 - 1) * (3 - 1)
nu.e = 4*3*(n-1)

(F.stat = to.f(lam3, nu.e, nu.h))
pf(F.stat, 2*nu.h, 2*(nu.e - 1), lower.tail = FALSE)

### 3e
D = diag(sqrt(c(11.896, 14.404, 13.656, 7245.6)))

# sowdate a.star
D %*% c(0.2089, -0.05317, 0.17493, -0.0042939)

# variety a.star
D %*% c(0.27131, -0.02672, 0.05558, 0.0017)

# interaction a.stars
D %*% c(0.28266, -0.04609, 0.02654, 0.00025815)
D %*% c(-0.06546287, -0.11341882, 0.20091054, 0.00350205)

### 4c
dat = matrix(scan("./MANDIBLE.DAT"), 18, 11, byrow=TRUE)

xbar = apply(dat[,3:11], 2, mean)

T.mat = matrix(c(-1, 0, 1, 1, -2, 1), 2, 3, byrow = TRUE)
(T.star = kronecker(T.mat, t(c(1, 1, 1))))

N = nrow(dat)
nu.e = N - 2
nu.h = 1
p = 9

S1 = var(dat[1:9,3:11])
S2 = var(dat[10:18,3:11])
Spl = 8*(S1 + S2)/(16)

# T^2 (test on linear and quadratic time trends)
N * t(T.star %*% xbar) %*% solve(T.star %*% (Spl %*% t(T.star))) %*% 
    (T.star %*% xbar)

# wilks lambda (test on time by group trends)
E = Spl * nu.e

x = dat[,3:11]
total = 0
for (i in 1:18)
    total = total + x[i,] %*% t(x[i,])
H = total - E - N * xbar %*% t(xbar)

(lam = det(T.star %*% E %*% t(T.star)) /
    det(T.star %*% (E + H) %*% t(T.star)))

(F.stat = (1-lam)/lam * (nu.e - p + 1)/p)
pf(F.stat, p, nu.e - p + 1, lower.tail = FALSE)

### 4d
# critical f value
qf(0.95, p, nu.e - p + 1, lower.tail = FALSE)
