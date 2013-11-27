require(SparkR)
require(MASS)
require('Matrix')
require(Rcpp)
require(svd)

args <- commandArgs(trailing = TRUE)

if (length(args) < 1) {
	print("Usage: DFC.R <master>[<slices>] <slices> <masked_file> <iterations>")
	q("no")
}

# Takes a list of columns, makes a matrix and applies SGD algorithm
factorCols <- function(itersAndMat) {
	#require('Matrix') # this is very important for some reason, probably should understand it
	t1 <- proc.time()
	iters <- itersAndMat[[1]][[1]]
	mat <- itersAndMat[[1]][[2]]
	UV <- apgBase(mat,iters)
	subproblemTime <- as.numeric((proc.time() - t1)["elapsed"])
	out <- list(UV, subproblemTime)
	list(out)
}

# Takes a list of factors of submatrices and projects them
# onto the column space of the first submatrix
# The factors (U,V) should be m-by-r and n-by-r respectively
dfcProject <- function(factorList) {
	t1 <- proc.time()
	U_1 <- factorList[[1]][[1]]
	V_1 <- factorList[[1]][[2]]
	#pseudoinverses
	U_1pinv <- ginv(U_1) 
	V_1pinv <- ginv(V_1)
	numParts <- length(factorList)
	partSize <- dim(U_1)[1] %/% numParts
	r <- dim(U_1)[2]
	# To be returned
	X_A <- U_1
	X_B <- V_1

	for (pair in tail(factorList,-1)){
		U_i <- pair[[1]]
		V_i <- pair[[2]]
		# We want to have U_1*Vhat_i = U_i*V_i, so we basically just solve
		Vhat_i <- t(((V_1pinv %*% V_1)%*%(U_1pinv %*% U_i)) %*% t(V_i))
		X_B <- rBind(X_B, Vhat_i)
	}
	projTime <- as.numeric((proc.time() - t1)["elapsed"])
	list(X_A, X_B,projTime)
}

dfcRandProject <- function(factorList) {
	t1 <- proc.time()
	V_1 <- factorList[[1]][[2]]
	slices <- length(factorList)
	partSize <- dim(V_1)[1]
	n <- slices * partSize # assume OK
	k <- dim(V_1)[2]
	
	# random projection default parameters
	p <- 5
	q <- 2
	
	# Random Gaussian matrix, break into chunks for simpler processing
	G <- Matrix(rnorm(n*(k+p),mean = 0,sd = 1),n,k+p)
	Glist <- lapply(1:slices, function(i) G[(1 + floor((i-1)*n/slices)):floor(i*n/slices),,drop=FALSE])
	
	# Initial QR factorization
	# Y = AG then factor Y = QR
	Ylist <- mapply(function(UV,G) UV[[1]] %*% (t(UV[[2]]) %*% G),factorList,Glist,SIMPLIFY = F)
	Y <- do.call(cBind,Ylist)
	QR <- qr(Y)
	Q <- qr.Q(QR)
	for (j in 1:q) {
		# Yhat = A'*Q then factor Yhat = Qhat*Rhat
		YhatList <- lapply(factorList, function(UV) UV[[2]] %*% (t(UV[[1]]) %*% Q))
		Yhat <- do.call(rBind,YhatList)
		QRhat <- qr(Yhat)
		Qhat <- qr.Q(QRhat)
		QhatList <- lapply(1:slices, function(i) Qhat[(1 + floor((i-1)*n/slices)):floor(i*n/slices),,drop=FALSE])
		
		# Y = A*Qhat then factor Y = Q*R
		Ylist <- mapply(function(UV,Qhat) UV[[1]] %*% (t(UV[[2]]) %*% Qhat),factorList,QhatList,SIMPLIFY =F)
		Y <- do.call(cBind,Ylist)
		QR <- qr(Y)
		Q <- qr.Q(QR)
	}
	
	# Take only the first k columns of Q
	Q <- Q[,1:k]
	
	# Finally project (Q*Q^+)*M
	Qpinv <- ginv(Q)
	Vlist <- lapply(factorList, function(UV) (Qpinv %*% UV[[1]]) %*% t(UV[[2]]))
	V <- t(do.call(cBind,Vlist))
	randprojTime <- as.numeric((proc.time() - t1)["elapsed"])
	list(Q,V,randprojTime) 
}

apgBase <- function(mat,maxiter) {
	# load required packages
	require('Matrix')
	require(Rcpp)
	require(svd)
	# Load and compile the fast C++ code
	sourceCpp('maskUV.cpp')
	# should figure out how to do all of the above better
	######## Set Initial Parameters #####################################
	m <- dim(mat)[1]
	n <- dim(mat)[2]
	IIJJ <- which(mat != 0,arr.ind = T) # list of nonzero indices
	p <- length(IIJJ) # number of nonzero entries
	II <- IIJJ[,1] # nonzero row indices
	JJ <- IIJJ[,2] # nonzero col indices
	L <- 1 # Lipschitz constant for 1/2*||Ax - b||_2^2
	t <- 1 # speed time [BC13] parameter
	told <- t
	beta <- 0 # beta = (told - 1)/t
	num_sv <- 5 # number of SV to look at
	num_pos_sv <-5 # initial number of positive singular values
	
	U <- Matrix(0,m,1) # Factor of X
	Uold <- U
	V <- Matrix(0,n,1) # Factor of X
	Vold <- V
	mX <- sparseMatrix(m,n,x=0) # Sparse matrix containing predicted values
	mXold <- mX # mX of previous iteration
	mY <- mX # Sparse matrix "average" of Xold and X
	
	mu0 <- norm(mat,type="F")
	mu <- 0.1*mu0
	muTarget <- 10^(-4)*mu0
	cat("mu :", mu, "\n")
	######################################################################

	for(iter in 1:maxiter) {
		cat("iteration: ",iter,"\n")
		# Get query access to G = Y - 1/L*Grad
		Grad <- mY - mat
		f <- function(z) as.numeric((1+beta)*(U %*% (t(V) %*% z)) - beta*(Uold %*% (t(Vold) %*% z)) - 1/L*(Grad %*% z)) # query oracle to Gk
		tf <- function(z) as.numeric((1+beta)*(V %*% (t(U) %*% z)) - beta*(Vold %*% (t(Uold) %*% z)) - 1/L*(t(Grad) %*% z))
		G <- extmat(f, tf, m, n)
		
		# Compute partial SVD
		svd <- propack.svd(G, neig = num_sv)
		
		# Update Params
		Uold <- U
		Vold <- V
		mXold <- mX
		told <- t
		
		s <- svd$d
		Shlf <- sqrt(s[which(s > mu/L)])
		if(num_sv == num_pos_sv) {
			num_sv <- num_pos_sv + 5
		}
		else {
			num_sv <- num_pos_sv + 1
		}
		cat("num sv: ",num_sv,"\n")
		# update number of positive singular values of X^k AFTER the above test
		num_pos_sv <- length(Shlf)
		
		Sig <- diag(x = Shlf,num_pos_sv,num_pos_sv)
		
		U <- (svd$u[,1:num_pos_sv] %*% Sig)
		V <- (svd$v[,1:num_pos_sv] %*% Sig)
		mX <- sparseMatrix(i=II,j=JJ,x=maskUV(as.matrix(U),as.matrix(V),II,JJ)) # Call into C++ code
		t <- (1+sqrt(1+4*t^2))/2
		beta <- (told - 1)/t
		mY <- (1+beta)*mX - beta*mXold	
		mu <- max(0.7*mu,muTarget)
		cat("mu: ",mu,"\n")
	}
	cat("U: ", dim(U),"\n")
	cat("V: ", dim(V),"\n")
	cat("RMSE for submatrix: ",errorCal(mat,U,V),"\n")
	list(U,V)
}
	

# Base stochastic gradient descent algorithm for matrix completion
sgdBase <- function(mat) {
	# Set Parameters
	m <- dim(mat)[1]
	n <- dim(mat)[2]
	lrate <- .04 # learning rate
	k <- .04 # parameter used to minimize over-fitting
	min_impr <- .001 # min improvement
	init <- 0.2 # initial value for features
	rank <- 10 # rank of feature vector
	min_itrs <- 10

	# Initialize
	minval <- min(mat)
	maxval <- max(mat)
	row_feats <- matrix(rnorm(rank*m,mean=0,sd = 0.2/sqrt(sqrt(rank))),rank,m)
	col_feats <- matrix(rnorm(rank*m,mean = 0,sd = 0.2/sqrt(sqrt(rank))),rank,n)
	rmse <- 2.0 # set rmse
	rmse_prev <- 2.0 # set previous rmse

	# Find nonzero entries
	nonzero_rowscols <- which(mat != 0,arr.ind = T)
	nonzero_rows <- nonzero_rowscols[,1]
	nonzero_cols <- nonzero_rowscols[,2]
	nonzero_entries <- mat[nonzero_rowscols]
	num_nonzeros <- length(nonzero_entries)
	

	
	for(i in 1:rank) {
		cat("rank: ", i, "\n")
		t <- 0
		impr <- 0.0
		while(t < min_itrs || impr > min_impr) {
			sq_err <- 0.0
			pred_time <- system.time(0)
			err_time <- system.time(0)
			update_time <- system.time(0)
			
			for(j in 1:num_nonzeros) {
				# find predicted val
				#start_pred <- proc.time()
				predval <- t(row_feats[,nonzero_rows[j] ]) %*% col_feats[,nonzero_cols[j] ]
				#pred_time <- pred_time + proc.time() - start_pred
				
				# apply cut off
				if(predval < minval) { predval <- minval }
				if(predval > maxval) { predval <- maxval }
				
				# Find Error
				#start_err <- proc.time()
				err <- nonzero_entries[j] - predval
				sq_err <- sq_err + err*err + k/2.0 * ((row_feats[i, nonzero_rows[j] ])^2) * ((col_feats[i, nonzero_cols[j] ])^2)
				#err_time <- err_time + proc.time() - start_err
				
				# Update row and col features
				#start_update <- proc.time()
				new_row_feat <- (1-lrate*k)*row_feats[i, nonzero_rows[j] ] + lrate*err*col_feats[i, nonzero_cols[j] ]
				new_col_feat <- (1-lrate*k)*col_feats[i, nonzero_cols[j] ] + lrate*err*row_feats[i, nonzero_rows[j] ]
				row_feats[i, nonzero_rows[j] ] <- new_row_feat
				col_feats[i, nonzero_cols[j] ] <- new_col_feat
				#update_time <- update_time + proc.time() - start_update
			}
			# Calculate RMSE
			rmse_prev <- rmse
			rmse <- sqrt(sq_err/num_nonzeros)
			#cat("pred_time: ",pred_time,"\n")
			#cat("err_time: ",err_time,"\n")
			#cat("update_time: ",update_time,"\n")
			#cat("fortime: ", fortime,"\n")
			cat("root mean squared error: ",rmse)
			cat("\n")
			impr <- rmse_prev - rmse
			t <- t + 1
		}
	}
	cat("RMSE for submatrix: ",rmse,"\n")
	list(row_feats,col_feats)
}

# Calculate the root mean squared error of the matrix 
# factorization UV' on the non-zero entries of mat
errorCal <- function(mat, U, V){
	# Find nonzero entries
	IIJJ <- which(mat != 0,arr.ind = T)
	numNonzero <- length(IIJJ)
	II <- IIJJ[,1]
	JJ <- IIJJ[,2]
	mX <- sparseMatrix(i=II,j=JJ,x=maskUV(as.matrix(U),as.matrix(V),II,JJ)) # Call into C++ code
	# Frobenius norm/sqrt(num_nonzero) = root mean squared error
	rmse <- norm(mat - mX,type = 'F')/sqrt(numNonzero)
	rmse
}

# Divide factor combine
dfc <- function(mat, sc, slices, iters) {
	sourceCpp('maskUV.cpp')
	# pick a random permutation of the columns
	cols <- dim(mat)[2]
	#sampleCols <- sample(cols)
	
	# make the matrix into a list of column chunks
	listMat <- lapply(1:slices, function(i) list(iters,mat[,(1 + floor((i-1)*cols/slices)):floor(i*cols/slices),drop=FALSE]))
	
	# distribute the column slices with spark
	# might need to pass in desired num slices here?
	
	t1 <- proc.time()
	subMatRDD <- parallelize(sc,listMat,slices)
	overhead <- as.numeric((proc.time() - t1)["elapsed"])
	
	# factor each slice
	factorsRDD <- lapplyPartition(subMatRDD,factorCols)
	
	# collect the results and project them onto the first column slice
	factorList <- collect(factorsRDD)
	matrixList <- lapply(seq(1,length(factorList)), function(i) factorList[[i]][[1]])
	subTimeList <- lapply(seq(1,length(factorList)), function(i) factorList[[i]][[2]])
	subTime <- max(unlist(subTimeList))
	cat("Time for subproblems: \n")
	print(subTimeList)

	# collect the results and project them onto the first column slice
	result <- dfcRandProject(matrixList)
	projTime <- result[[3]]
	cat("Time for collection: ",projTime,"\n")

	# get the error
	#print(result[[1]]%*%t(result[[2]]))
	error <- errorCal( mat, result[[1]], result[[2]])
	list(error,overhead,subTime,projTime)
}

sc <- sparkR.init(args[[1]], "DFCR")

slices <- ifelse(length(args) > 1, as.integer(args[[2]]),2)


# Test matrix init
#dims <- 100
#r <- 10
#testU <- matrix(rnorm(r*dims,mean = 4,sd = 1),r,dims)
#testV <- matrix(rnorm(r*dims,mean = 4,sd = 1),r,dims)
#testG <- matrix(rnorm(dims*dims,mean = 0,sd = .05),dims,dims)
#testM <- t(testV) %*% testU + testG

#zeros <- sample.int(dims^2, size=floor(0.7*dims^2)) 

#zeroedM <- testM
#for(i in zeros){
	#x = i %% dims
	#y = i %/% dims
	#zeroedM[x,y] <- 0
#}

# Read matrix from file
maskedFile <- args[[3]]
iterations <- args[[4]]
#trueVFile <- args[[5]]

maskedM <- readMM(maskedFile)
dims <- dim(maskedM)[1]
revealedEntries <- nnzero(maskedM)
#trueU <- read(trueUFile)
#trueV <- read(trueVFile)

# Run DFC
t1 <- proc.time()
outs <- dfc(maskedM, sc, slices, iterations)
totalTime <- as.numeric((proc.time() - t1)["elapsed"])
cat("RMSE for the entire matrix: ",outs[[1]],"\n")
cat("Total time for DFC: ",totalTime,"\n")
#cat("Average magnitude of entries of M: ",sum(abs(maskedM))/revealedEntries,"\n")
outs
