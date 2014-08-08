###	Algorithm 1: Go through each rows set of permutations one by one
###		until a solution is found
###	Algorithm 2: Start at the row or column with the least number of
###		permutations. Go through that row/column permutation set
###		and find the points of intersection. Update the board with
###		the common points. Go to next row/column. Keep updating and
###		recylcing through past rows/columns until the board is done.
###	Algorithm 3: Same as Algorithm 2 except a change in how we go through
###		permutations for a row/column. Instead of going through every
###		permutation, select one segment in a row/column, move all other
###		segments as far away as possible (with the given board) and
###		going through the permutations for that segment and keep all
###		pieces in common. This should reduce the number of iterations
###		significantly.


random.board=function(x,y,p=0.5){ matrix(sample(c(0,1),x*y,T,c(1-p,p)),x,y) }

rows.cols=function(board, rows=TRUE){
	out=list()
	if (rows) rc=c(1,2) else rc=c(2,1)
	n=dim(board)[rc[1]]
	m=dim(board)[rc[2]]
	at=0
	if (rc[2]==1) board=t(board)
	for (i in 1:n){
		count=0
		part=NULL
		for (at in 1:m){
			if (at==which.max(board[i,])){
				if (board[i,at]==1)
					count=count+1
				board[i,at]=-1
				if (at+1 != which.max(board[i,]) || at==m){
					part[length(part)+1]=count
					count=0
					}
				}
			}
		out[[length(out)+1]]=part
		}
	return(out)
	}

### Making this run faster will mean solving faster
permutate=function(x, size, permutation.number){
	y=permutation.number
	n=length(x)
	zeros=integer(n)
	out=integer(size)
	MAX = size + 1 - sum(x) - n
	if (n==1){
		zeros[1]=y
		y=0
		if ( zeros[1] > MAX )
			zeros[1] = MAX
		}
	while (y > 0){
		y = y - 1
		zeros[n] = zeros[n] + 1
		for (i in n:2){
			if ( zeros[i] > MAX - sum(zeros[1:(i-1)]) ){
				zeros[i] = 0
				zeros[i-1] = zeros[i-1] + 1
				}
			}
			if ( zeros[1] > MAX ){
				zeros[1] = MAX
				break
				}
		}
	zeros=zeros+1
	zeros[1]=zeros[1]-1
	if (x[1] != 0)
		out[(zeros[1]+1):(zeros[1]+x[1])]=1
	if (n > 1){
		for (i in 2:n){
			out[(sum(zeros[1:i])+sum(x[1:(i-1)])+1):
				(sum(zeros[1:i])+sum(x[1:(i-1)])+x[i])]=1
			}
		}
	return (out*2-1)
	}

legal=function(board, rows, cols){
	nrow=length(rows)
	ncol=length(cols)

	# Row correctness check
	for (i in 1:nrow){
		temp_row = NULL
		adding=0
		for (j in 1:ncol){
			if (board[i,j]==1)
				adding=adding+1
			if (adding>0 && (board[i,j]==0 || j==ncol)){
				temp_row[length(temp_row)+1]=adding
				adding=0
				temp_row[length(temp_row)+1]=0
				}
			if (j==ncol){
				if (length(temp_row)==0)
					temp_row[1:2]=c(0,0)
				if (temp_row[length(temp_row)]==0)
					temp_row=temp_row[-length(temp_row)]
				}
			}
		# Make sure there is the right number of 1's
		if (sum(temp_row) != sum(rows[[i]]))
			return (FALSE)
		
		# Make sure the row is properly divided
		if (!all(temp_row[temp_row>0]==rows[[i]]))
			return (FALSE)
		}

	# Column correctness check
	for (j in 1:ncol){
		temp_col = NULL
		adding=0
		for (i in 1:nrow){
			if (board[i,j]==1)
				adding=adding+1
			if (adding>0 && (board[i,j]==0 || i==nrow)){
				temp_col[length(temp_col)+1]=adding
				adding=0
				temp_col[length(temp_col)+1]=0
				}
			if (i==nrow){
				if (length(temp_col)==0)
					temp_col[1:2]=c(0,0)
				if (temp_col[length(temp_col)]==0)
					temp_col=temp_col[-length(temp_col)]
				}
			}
		# Make sure there is the right number of 1's
		if (sum(temp_col) != sum(cols[[j]]))
			return (FALSE)
		
		# Make sure the row is properly divided
		if (!all(temp_col[temp_col>0]==cols[[j]]))
			return (FALSE)
		}
	return (TRUE)
	}

plot.during=function(out, init=FALSE, SIZE=9){
	require(rgl)
	x=1:dim(out)[1]
	y=1:dim(out)[2]
	black <- red <- matrix(0,max(x)*max(y),3)
	black[,2] <- red[,2] <- rep(x,max(y))
	red[,3]=NA
	for (i in 1:(max(x)*max(y))){
		black[i,1] <- red[i,1] <- y[ceiling(i/max(x))]
		black[i,3] = out[black[i,2],black[i,1]]
		if (black[i,3] <= 0){
			if (black[i,3]==-1)
				red[i,3]=1
			black[i,3]=NA		
			}
		}
	black[,2] = black[length(black[,2]):1,2]
	red[,2] = red[length(red[,2]):1,2]
	if (!init){
		if (!all(is.na(black[,3])))
			points3d(black,size=SIZE)
		if (!all(is.na(red[,3])))
			points3d(red,size=SIZE,col='red')
	} else {
		plot3d(matrix(c(1,1,1),1,3),xlim=c(0.5,max(c(x,y))+0.5),ylim=c(0.5,max(c(x,y))+0.5),
			zlim=c(1,1),type='n')
		}
	}

solver2=function(rows, cols, wait=TRUE){

	# Add option to give initial board
	# Add solving method that uses contradiction

	# Initialize the board
	LENrow = length(rows)
	LENcol = length(cols)
	board = matrix(0,LENrow, LENcol)
	prevboard = board

	# Get max permutation number
	max.perm=integer(LENrow + LENcol)
	rows.updated=list()
	cols.updated=list()
	for (i in 1:LENrow){
		max.perm[i]=choose(LENcol-sum(rows[[i]])+1, length(rows[[i]]))
		rows.updated[[length(rows.updated)+1]] = integer(LENcol)
		}
	for (i in 1:LENcol){
		max.perm[i + LENrow]=choose(LENrow-sum(cols[[i]])+1, length(cols[[i]]))
		cols.updated[[length(cols.updated)+1]]=integer(LENrow)	
		}

	good.perm=list()
	first.update=integer(LENrow + LENcol)
	for (i in 1:length(max.perm)){
		good.perm[[length(good.perm)+1]]=0:(max.perm[i]-1)
		first.update[i]=1
		}

	at=0
	position = 0

	out=c(sum(max.perm)/max(LENrow,LENcol), 1)

#	print(paste0("Difficulty: ",sum(max.perm)/max(LENrow,LENcol)))

#	plot.during(board, TRUE)
#	if (wait)
#		readline()
	while ( !legal(floor((board+1)/2),rows,cols) ){
		position = position + 1
		if (position > length(max.perm)){
			print("Solver could not find a solution.")
			out[2]=0
			break
			}

	# Need to figure out how to remove certain permutation numbers
	# after being given new information without haaving to iterate
	# through all the known good permutations. Look at how choose(n, k)
	# works with certain parts in vector moving

	# See if I can determine whether or not a change will happen
	# by looking at length, sum, parts, max number of permutations, and
	# what is given

	# Put priority on the ones with more given

		at=order(max.perm)[position]
		if (at <= LENrow){
			row_bool=TRUE
		} else {
			row_bool=FALSE
			}

		if ( row_bool && (!all(board[at,]==rows.updated[[at]]) || first.update[at]) ){
#			print(c(at, LENrow, max.perm[at]))
			first.update[at]=0
			end=integer(LENcol)
			change=0
			for (x in good.perm[[at]]){
				compare = permutate(rows[[at]],LENcol,x)
				# Check for what is known in the board
				for (c in 1:LENcol){
					if (board[at,c] != 0){
						if (board[at,c] != compare[c]){
							good.perm[[at]]=good.perm[[at]][-which(good.perm[[at]]==x)]
							break
							}
						}
					if (c==LENcol){
						end = end + compare
						change=change+1
						}
					}
				}
			for (c in 1:LENcol){
				if (abs(end[c])==change && change != 0){
					board[at,c]=end[c]/change
					}
				}
			rows.updated[[at]] = board[at,]
			}
		if (!row_bool && (!all(board[,at-LENrow]==cols.updated[[at-LENrow]]) || first.update[at]) ){
#			print(c(at, LENrow, max.perm[at]))
			first.update[at]=0
			end=integer(LENrow)
			change=0
			for (x in good.perm[[at]]){
				compare = permutate(cols[[at-LENrow]],LENrow,x)
				# Check for what is known in the board
				for (r in 1:LENrow){
					if (board[r,at-LENrow] != 0){
						if (board[r,at-LENrow] != compare[r]){
							good.perm[[at]]=good.perm[[at]][-which(good.perm[[at]]==x)]
							break
							}
						}
					if (r==LENrow){
						end = end + compare
						change=change+1
						}
					}
				}
			for (r in 1:LENrow){
					if (abs(end[r])==change && change != 0){
					board[r,at-LENrow]=end[r]/change
					}
				}
			cols.updated[[at-LENrow]] = board[,at-LENrow]
			}
		temp = max.perm[at]
		max.perm[at] = length(good.perm[[at]])

		if ( !all(prevboard == board) ){
#			plot.during(board)
			prevboard=board
			if ( max.perm[at] != temp ){
				position=length(max.perm[max.perm==1])
				}
			}
		}
#	plot.during(board)
	return(out)
	}


### Making this run faster will mean solving faster
permutate3=function(x, part, given, permutation.number){
	y=permutation.number
	n=length(x)
	(1:n)[-part]

	}


solver3=function(rows, cols, wait=TRUE){

	# Initialize the board
	LENrow = length(rows)
	LENcol = length(cols)
	board = matrix(0,LENrow, LENcol)
	prevboard = board

	# Get max permutation number
	max.perm=integer(LENrow + LENcol)
	rows.updated=list()
	cols.updated=list()
	for (i in 1:LENrow){
		max.perm[i]=choose(LENcol-sum(rows[[i]])+1, length(rows[[i]]))
		rows.updated[[length(rows.updated)+1]] = integer(LENcol)
		}
	for (i in 1:LENcol){
		max.perm[i + LENrow]=choose(LENrow-sum(cols[[i]])+1, length(cols[[i]]))
		cols.updated[[length(cols.updated)+1]]=integer(LENrow)	
		}

	good.perm=list()
	first.update=integer(LENrow + LENcol)
	for (i in 1:length(max.perm)){
		good.perm[[length(good.perm)+1]]=0:(max.perm[i]-1)
		first.update[i]=1
		}

	at=0
	position = 0

	out=c(sum(max.perm)/max(LENrow,LENcol), 1)

#	print(paste0("Difficulty: ",sum(max.perm)/max(LENrow,LENcol)))

#	plot.during(board, TRUE)
#	if (wait)
#		readline()
	while ( !legal(floor((board+1)/2),rows,cols) ){
		position = position + 1
		if (position > length(max.perm)){
			print("Solver could not find a solution.")
			out[2]=0
			break
			}

	# Need to figure out how to remove certain permutation numbers
	# after being given new information without haaving to iterate
	# through all the known good permutations. Look at how choose(n, k)
	# works with certain parts in vector moving

	# See if I can determine whether or not a change will happen
	# by looking at length, sum, parts, max number of permutations, and
	# what is given

	# Put priority on the ones with more given

		at=order(max.perm)[position]
		if (at <= LENrow){
			row_bool=TRUE
		} else {
			row_bool=FALSE
			}

		if ( row_bool && (!all(board[at,]==rows.updated[[at]]) || first.update[at]) ){
#			print(c(at, LENrow, max.perm[at]))
			first.update[at]=0
			end=integer(LENcol)
			change=0
			for (x in good.perm[[at]]){
				compare = permutate(rows[[at]],LENcol,x)
				# Check for what is known in the board
				for (c in 1:LENcol){
					if (board[at,c] != 0){
						if (board[at,c] != compare[c]){
							good.perm[[at]]=good.perm[[at]][-which(good.perm[[at]]==x)]
							break
							}
						}
					if (c==LENcol){
						end = end + compare
						change=change+1
						}
					}
				}
			for (c in 1:LENcol){
				if (abs(end[c])==change && change != 0){
					board[at,c]=end[c]/change
					}
				}
			rows.updated[[at]] = board[at,]
			}
		if (!row_bool && (!all(board[,at-LENrow]==cols.updated[[at-LENrow]]) || first.update[at]) ){
#			print(c(at, LENrow, max.perm[at]))
			first.update[at]=0
			end=integer(LENrow)
			change=0
			for (x in good.perm[[at]]){
				compare = permutate(cols[[at-LENrow]],LENrow,x)
				# Check for what is known in the board
				for (r in 1:LENrow){
					if (board[r,at-LENrow] != 0){
						if (board[r,at-LENrow] != compare[r]){
							good.perm[[at]]=good.perm[[at]][-which(good.perm[[at]]==x)]
							break
							}
						}
					if (r==LENrow){
						end = end + compare
						change=change+1
						}
					}
				}
			for (r in 1:LENrow){
					if (abs(end[r])==change && change != 0){
					board[r,at-LENrow]=end[r]/change
					}
				}
			cols.updated[[at-LENrow]] = board[,at-LENrow]
			}
		temp = max.perm[at]
		max.perm[at] = length(good.perm[[at]])

		if ( !all(prevboard == board) ){
#			plot.during(board)
			prevboard=board
			if ( max.perm[at] != temp ){
				position=length(max.perm[max.perm==1])
				}
			}
		}
#	plot.during(board)
	return(out)
	}


times=double(100)
difficulty=matrix(0,100,2)

for (i in 1:100){
	print(i)
	test=random.board(sample(15:25,1),sample(15:25,1),runif(1,0.6,0.8))
	rows=rows.cols(test)
	cols=rows.cols(test,F)

	times[i] = system.time(
		difficulty[i,] <- solver(rows, cols))[3]
	}

plot(difficulty)
plot(difficulty[difficulty[,2]==1,1],times[difficulty[,2]==1],
	ylim=c(min(times),max(times)),xlim=c(min(difficulty[,1]),max(difficulty[,1])))
plot(difficulty[difficulty[,2]==0,1],times[difficulty[,2]==0],
	ylim=c(min(times),max(times)),xlim=c(min(difficulty[,1]),max(difficulty[,1])))
plot(difficulty[,1],times,
	ylim=c(min(times),max(times)),xlim=c(min(difficulty[,1]),max(difficulty[,1])))

Y=times[c(1:71,73:100)]
X=difficulty[c(1:71,73:100),]
plot(X[,1],Y)

plot(difficulty[,1],log(times))
X=difficulty[,1]
Y=log(times)
mod=lm(Y~X)
plot(X,Y)
abline(coef(mod)[1],coef(mod)[2])

(out=solver(Row8.29, Col8.29))

common=function(vec, given){
	i=0
	out=NULL
	LEN=length(given)
	max=choose(LEN - sum(vec) + 1,length(vec))
	for (i in 0:(max-1)){
		perm=permutate(vec, LEN, i)
		for (r in 1:LEN){
			if (given[r] != 0){
				if (perm[r] != given[r]){
					break
					}
				}
			if (r==LEN){
				out[length(out)+1]=i
				}
			}
		}
	return(out)
	}

### Max permutations: 36
(permutate(c(1,1),10,0)+1)/2
(permutate(c(1,1),10,1)+1)/2
(permutate(c(1,1),10,2)+1)/2
(permutate(c(1,1),10,3)+1)/2
(permutate(c(1,1),10,4)+1)/2
(permutate(c(1,1),10,5)+1)/2
(permutate(c(1,1),10,6)+1)/2
(permutate(c(1,1),10,7)+1)/2

(permutate(c(1,1),10,8)+1)/2
(permutate(c(1,1),10,9)+1)/2
(permutate(c(1,1),10,10)+1)/2
(permutate(c(1,1),10,11)+1)/2
(permutate(c(1,1),10,12)+1)/2
(permutate(c(1,1),10,13)+1)/2
(permutate(c(1,1),10,14)+1)/2

(permutate(c(1,1),10,15)+1)/2
(permutate(c(1,1),10,16)+1)/2
(permutate(c(1,1),10,17)+1)/2
(permutate(c(1,1),10,18)+1)/2
(permutate(c(1,1),10,19)+1)/2
(permutate(c(1,1),10,20)+1)/2

(permutate(c(1,1),10,21)+1)/2
(permutate(c(1,1),10,22)+1)/2
(permutate(c(1,1),10,23)+1)/2
(permutate(c(1,1),10,24)+1)/2
(permutate(c(1,1),10,25)+1)/2

(permutate(c(1,1),10,26)+1)/2
(permutate(c(1,1),10,27)+1)/2
(permutate(c(1,1),10,28)+1)/2
(permutate(c(1,1),10,29)+1)/2

(permutate(c(1,1),10,30)+1)/2
(permutate(c(1,1),10,31)+1)/2
(permutate(c(1,1),10,32)+1)/2

(permutate(c(1,1),10,33)+1)/2
(permutate(c(1,1),10,34)+1)/2

(permutate(c(1,1),10,35)+1)/2

# Groups = 8 = 10 - 2 = Length - Sum?

# permutations 0:7 are			1 0 X X X X X X X X
# permutations 8:14 are			0 1 0 X X X X X X X
# permutations 15:20 are		X 0 1 0 X X X X X X
# permutations 21:25 are		X X 0 1 0 X X X X X
# permutations 26:29 are		X X X 0 1 0 X X X X
# permutations 30:32 are		X X X X 0 1 0 X X X
# permutations 33:34 are		X X X X X 0 1 0 X X
# permutation 35 is			X X X X X X 0 1 0 1

# 8:35, are					0 X X X X X X X X X (skip group 1)
# 0:7, 14:35 are				X 0 X X X X X X X X (skip group 2)
# 1:14, 21:35 are				X X 0 X X X X X X X (skip group 3,
										1st in group 1)
# 0,2:7,9:14,15:20,26:35 are		X X X 0 X X X X X X (skip group 4,
										2nd in group 1,
										1st in group 2)
# 0:1,3:8,10:14,16:25,30:35 are	X X X X 0 X X X X X (skip group 5,
										3rd in group 1,
										2nd in group 2,
										1st in group 3)
						X X X X X 0 X X X X (skip group 6,
										4th in group 1,
										3rd in group 2,
										2nd in group 3,
										1st in group 4)
						X X X X X X 0 X X X (skip group 7,
										5th in group 1,
										4th in group 2,
										3rd in group 3,
										2nd in group 4,
										1st in group 5)
						X X X X X X X 0 X X (skip group 8,
										6th in group 1,
										5th in group 2,
										4th in group 3,
										3rd in group 4,
										2nd in group 5,
										1st in group 6)
						X X X X X X X X 0 X (skip 2nd to last element in each group:
										7th in group 1,
										6th in group 2,
										5th in group 3,
										4th in group 4,
										3rd in group 5,
										2nd in group 6,
										1st in group 7,
						X X X X X X X X X 0 (skip last elemnt in each group:
										8th in group 1,
										7th in group 2,
										6th in group 3,
										5th in group 4,
										4th in group 5,
										3rd in group 6,
										2nd in group 7,
										1st in group 8,


# 7,14,20,25,29,32,34,35 are	X X X X X X X X X 1 (last from each group)
			implies	X X X X X X X X 0 1




(permutate(c(1,2),10,0)+1)/2
(permutate(c(1,2),10,1)+1)/2
(permutate(c(1,2),10,2)+1)/2
(permutate(c(1,2),10,3)+1)/2
(permutate(c(1,2),10,4)+1)/2
(permutate(c(1,2),10,5)+1)/2
(permutate(c(1,2),10,6)+1)/2

(permutate(c(1,2),10,7)+1)/2
(permutate(c(1,2),10,8)+1)/2
(permutate(c(1,2),10,9)+1)/2
(permutate(c(1,2),10,10)+1)/2
(permutate(c(1,2),10,11)+1)/2
(permutate(c(1,2),10,12)+1)/2

(permutate(c(1,2),10,13)+1)/2
(permutate(c(1,2),10,14)+1)/2
(permutate(c(1,2),10,15)+1)/2
(permutate(c(1,2),10,16)+1)/2
(permutate(c(1,2),10,17)+1)/2

(permutate(c(1,2),10,18)+1)/2
(permutate(c(1,2),10,19)+1)/2
(permutate(c(1,2),10,20)+1)/2
(permutate(c(1,2),10,21)+1)/2

(permutate(c(1,2),10,22)+1)/2
(permutate(c(1,2),10,23)+1)/2
(permutate(c(1,2),10,24)+1)/2

(permutate(c(1,2),10,25)+1)/2
(permutate(c(1,2),10,26)+1)/2

(permutate(c(1,2),10,27)+1)/2

# Groups = 7 = 10 - 3 = Length - Sum?

Given:
0 X X X X X X X X -- Skip group 1
X 0 X X X X X X X -- Skip group 2
X X 0 X X X X X X -- Skip group 3 and
				group 1's [1]
X X X 0 X X X X X -- Skip group 4 and
				group 1's [1], [2]
				group 2's [1]
X X X X 0 X X X X -- Skip group 5 and
				group 1's [2], [3]
				group 2's [1], [2]
				group 3's [1]
X X X X X 0 X X X -- Skip group 6 and
				group 1's [3], [4]
				group 2's [2], [3]
				group 3's [1], [2]
				group 4's [1]




