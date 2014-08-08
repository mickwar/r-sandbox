colfun = function(x, n){
	out = double(3)
	out[1] = -1 + x * 2/n
	out[2] = 1 - 2/n * abs(x - n/2)
	out[3] = 1 - x * 2/n
	for (i in 1:3){
		if (out[i] > 1 || out[i] < 0)
			out[i] = 0
		}
	return(out)
	}

rgb2 = function(x)
	rgb(x[1], x[2], x[3])

n=1000
xx=seq(0,1,length=n+1)

plot(0, type='n', xlim=c(0, 1))
for (i in 0:n){
	points(xx[i+1], 0, col=rgb2(colfun(i, n)))
	}

#######

color.manip = function(color, ...){
	R = color[1]
	G = color[2]
	B = color[3]
	M = max(color)
	m = min(color)
	actions = list(...)
	if (is.null(actions[["hue"]])){
		hue = atan2((sqrt(3)*(G-B))*0.5, (2*R-G-B)*0.5)
		hue = if (hue < 0) 2*pi+hue else hue
	} else {
		hue = actions[["hue"]]
		}
	if (is.null(actions[["sat"]])){
		sat = if (M > 0) (M-m)/M else 0
	} else {
		sat = actions[["sat"]]
		}
	if (is.null(actions[["val"]])){
		val = M
	} else {
		val = actions[["val"]]
		}
	C = val * sat
	H = hue / (pi/3)
	X = C * (1-abs(H %% 2 - 1))
	out = double(3)
	if (H >= 0 && H < 1)
		out = c(C, X, 0)
	if (H >= 1 && H < 2)
		out = c(X, C, 0)
	if (H >= 2 && H < 3)
		out = c(0, C, X)
	if (H >= 3 && H < 4)
		out = c(0, X, C)
	if (H >= 4 && H < 5)
		out = c(X, 0, C)
	if (H >= 5 && H < 6)
		out = c(C, 0, X)
	m = val - C
	out = out + m
	if (is.null(actions[["return"]]))
		return (out)
	return(list("RGB"=out,"HSV"=c(hue, sat, val)))
	}

color.manip(c(0.3,0.5,0.7), return=1, sat = 1)

A = c(0.0000, 0.0000, 0.7000)
B = c(0.9333, 0.3764, 0.0117)
C = c(0.0000, 0.7000, 0.0000)
D = c(0.8901, 0.0901, 0.0509)
E = c(0.1019, 0.1019, 0.1019)
F = c(0.0000, 0.5450, 0.5450)
G = c(0.8039, 0.6784, 0.0000)
H = c(0.2941, 0.0000, 0.5098)

aa = c(0.1137, 0.5254, 0.9333)
bb = c(1.0000, 0.6588, 0.1411)
cc = c(0.3000, 0.9000, 0.0000)
dd = c(0.9607, 0.4705, 0.3529)
ee = c(0.7294, 0.7294, 0.7294)
ff = c(0.0000, 0.9333, 0.9333)
gg = c(1.0000, 1.0000, 0.0000)
hh = c(0.7490, 0.3725, 1.0000)

AI = color.manip(A, return=1)
BI = color.manip(B, return=1)
CI = color.manip(C, return=1)
DI = color.manip(D, return=1)
EI = color.manip(E, return=1)
FI = color.manip(F, return=1)
GI = color.manip(G, return=1)
HI = color.manip(H, return=1)

newA = color.manip(A, sat = AI$HSV[2]/2)
newB = color.manip(B, sat = BI$HSV[2]/2)
newC = color.manip(C, sat = CI$HSV[2]/2)
newD = color.manip(D, sat = DI$HSV[2]/2)
newE = color.manip(E, val = EI$HSV[3]*6)
newF = color.manip(F, sat = FI$HSV[2]/2)
newG = color.manip(G, sat = GI$HSV[2]/2)
newH = color.manip(H, sat = HI$HSV[2]/2)


newA
newB
newC
newD
newE
newF
newG
newH

n = 300
all = matrix(0, n, 3)
all[1, ] = color.manip(c(0.1, 0.9, 0.5), val=1, sat=0.1)
for (i in 2:n)
	all[i, ] = color.manip(all[i-1, ], val=1-i/n, sat=i/n)

plot(1:n, 1:n, col=rgb(all), pch=20)
