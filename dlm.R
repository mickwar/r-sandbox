get.vals = function(dat, delta, m0, C0, n0, d0){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    n = length(y)
    p = length(m0)

    Rt = array(0, c(n+1, p, p))
    qt = double(n+1)
    At = matrix(0, n+1, p)
    ft = double(n+1)
    et = double(n+1)
    mt = matrix(0, n+1, p)
    nt = double(n+1)
    dt = double(n+1)
    st = double(n+1)
    Ct = array(0, c(n+1, p, p))

    mt[1,] = m0
    Ct[1,,] = C0
    nt[1] = n0
    dt[1] = d0
    st[1] = dt[1] / nt[1]

    for (i in 2:(n+1)){
#       Rt[i,,] = G.mat %*% Ct[i-1,,] %*% t(G.mat) + (1-delta)/delta*Ct[i-1,,]
        Rt[i,,] = 1/delta * G.mat %*% Ct[i-1,,] %*% t(G.mat)
        qt[i] = t(F.vec) %*% Rt[i,,] %*% F.vec + st[i-1]
        At[i,] = (Rt[i,,] %*% F.vec)/qt[i]   
        ft[i] = t(F.vec) %*% G.mat %*% mt[i-1,]
        et[i] = y[i-1] - ft[i]
        mt[i,] = G.mat %*% mt[i-1,] + At[i,] * et[i]
        nt[i] = nt[i-1] + 1
        dt[i] = dt[i-1] + st[i-1]*(et[i]^2) / qt[i]
        st[i] = st[i-1] + st[i-1] / nt[i]*( (et[i]^2) / qt[i] - 1)
        Ct[i,,] = (st[i] / st[i-1])*(Rt[i,,] - At[i,] %*% t(At[i,]) * qt[i])
        }

    
    return (list("Rt"=Rt[-1,,], "qt"=qt[-1], "At"=At[-1,],
        "ft"=ft[-1], "et"=et[-1], "mt"=mt[-1,], "nt"=nt[-1],
        "dt"=dt[-1], "st"=st[-1], "Ct"=Ct[-1,,], "delta"=delta,
        "m0"=m0, "C0"=C0, "n0"=n0, "d0"=d0))
    }
obs.pred.dens = function(dat, params){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    df = c(params$n0, params$nt[-length(params$nt)])
    mu = params$ft
    sig2 = params$qt
    out = lgamma((df+1)/2) - lgamma(df/2) - 1/2*log(df*pi*sig2) -
        (df+1)/2*log(1+(y - mu)^2 / (df*sig2))
    return (sum(out))
    }
smooth = function(dat, params){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    G.inv = solve(G.mat)
    n = length(y)
    p = length(F.vec)
    delta = params$delta
    out.at = matrix(0, n, p)
    out.Rt = array(0, c(n, p, p))
    out.at[n,] = params$mt[n,]
    out.Rt[n,,] = params$Ct[n,,]
    out.ft = double(n)
    out.qt = double(n)
    out.ft[n] = t(F.vec) %*% out.at[n,]
    out.ft[n] = t(F.vec) %*% out.at[n,]
    for (i in (n-1):1){
        out.at[i,] = (1-delta)*params$mt[i,] + delta*(G.inv %*% out.at[i+1,])
        out.Rt[i,,] = (1-delta)*params$Ct[i,,] + delta^2*(G.inv %*% out.Rt[i+1,,] %*% t(G.inv))
        out.ft[i] = t(F.vec) %*% out.at[i,]
        out.qt[i] = t(F.vec) %*% out.Rt[i,,] %*% F.vec + params$st[n]
        }
    return (list("at.s"=out.at, "Rt.s"=out.Rt, "ft.s"=out.ft, "qt.s"=out.qt))
    }
forecast = function(dat, params, h = 12){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    n = length(y)
    p = length(F.vec)
    delta = params$delta

#   at.p = rep(params$mt[n], h)
#   Rt.p = (1/delta^(1:h))*params$Ct[n]
#   ft.p = rep(params$mt[n], h)
#   qt.p = Rt.p + params$St[n]
#   return (list("at.p"=at.p, "Rt.p"=Rt.p,
#       "ft.p"=ft.p, "Qt.p"=Qt.p))

    at.p = matrix(0, h+1, p)
    Rt.p = array(0, c(h+1, p, p))
    ft.p = double(h+1)
    qt.p = double(h+1)

    at.p[1,] = params$mt[n,]
    Rt.p[1,,] = params$Ct[n,,]
    for (i in 2:(h+1)){
        at.p[i,] = G.mat %*% at.p[i-1,]
        Rt.p[i,,] = G.mat %*% (Rt.p[i-1,,] + (1-delta)/delta*params$Ct[n,,]) %*% t(G.mat)
        ft.p[i] = t(F.vec) %*% at.p[i,]
        qt.p[i] = t(F.vec) %*% Rt.p[i,,] %*% F.vec + params$st[n]
        }
    return (list("at.p"=at.p[-1,], "Rt.p"=Rt.p[-1,,],
        "ft.p"=ft.p[-1], "qt.p"=qt.p[-1]))
    }
