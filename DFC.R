require(SparkR)
require(MASS)
require('Matrix')

args <- commandArgs(trailing = TRUE)

if (length(args) < 1) {
	print("Usage: DFC.R <master> [<slices>] <masked_file> <U_file> <V_file>")
	q("no")
}


# Takes a list of columns, makes a matrix and applies SGD algorithm
factorCols <- function(colList,rows) {
	require('Matrix')
	print(colList)
	M <- do.call(cBind,colList)
	UV <- sgdBase(M)
	list(UV)
}

# Takes a list of factors of submatrices and projects them
# onto the column space of the first submatrix
dfcProject <- function(factorList) {
	U_1 <- factorList[[1]][[1]]
	V_1 <- factorList[[1]][[2]]
	#pseudoinverses
	U_1pinv <- ginv(U_1) 
	V_1pinv <- ginv(V_1)
	numParts <- length(factorList)
	partSize <- dim(U_1)[2] %/% numParts
	r <- dim(U_1)[1]
	# To be returned
	X_A <- U_1
	X_B <- V_1

	for (pair in tail(factorList,-1)){
		U_i <- pair[[1]]
		V_i <- pair[[2]]
		# We want to have U_1*Vhat_i = U_i*V_i, so we basically just solve
		Vhat_i <- ((V_1 %*% V_1pinv)%*%(t(U_1pinv) %*% t(U_i)))%*%V_i
		X_B <- cBind(X_B, Vhat_i)
	}
	list(X_A, X_B)

}
# Base stochastic gradient descent algorithm for matrix completion
sgdBase <- function(mat) {
	# Set Parameters
	m <- dim(mat)[1]
	n <- dim(mat)[2]
	lrate <- .004 # learning rate
	k <- .04 # parameter used to minimize over-fitting
	min_impr <- .01 # min improvement
	init <- .1 # initial value for features
	rank <- 10 # rank of feature vector
	min_itrs <- 10

	# Initialize
	minval <- min(mat)
	maxval <- max(mat)
	row_feats <- matrix(init,rank,m)
	col_feats <- matrix(init,rank,n)
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
			for(j in 1:num_nonzeros) {
				# find predicted val
				predval <- t(row_feats[,nonzero_rows[j] ]) %*% col_feats[,nonzero_cols[j] ]
				
				# apply cut off
				if(predval < minval) { predval <- minval }
				if(predval > maxval) { predval <- maxval }
				
				# Find Error
				err <- nonzero_entries[j] - predval
				sq_err <- sq_err + err*err
				
				# Update row and col features
				new_row_feat <- (1-lrate*k)*row_feats[i, nonzero_rows[j] ] + lrate*err*col_feats[i, nonzero_cols[j] ]
				new_col_feat <- (1-lrate*k)*col_feats[i, nonzero_cols[j] ] + lrate*err*row_feats[i, nonzero_rows[j] ]
				row_feats[i, nonzero_rows[j] ] <- new_row_feat
				col_feats[i, nonzero_cols[j] ] <- new_col_feat
			}
			# Calculate RMSE
			rmse_prev <- rmse
			rmse <- sqrt(sq_err/num_nonzeros)
			#cat("root mean squared error: ",rmse)
			#cat("\n")
			impr <- rmse_prev - rmse
			t <- t + 1
		}
	}
	cat("RMSE for submatrix: ",rmse,"\n")
	list(row_feats,col_feats)
}

errorCal <- function(mat, row_pred, col_pred){
	# Find nonzero entries
	nonzero_rowscols <- which(mat != 0,arr.ind = T)
	nonzero_rows <- nonzero_rowscols[,1]
	nonzero_cols <- nonzero_rowscols[,2]
	nonzero_entries <- mat[nonzero_rowscols]
	num_nonzeros <- length(nonzero_entries)
	
	sq_err <- 0.0	
	for(j in 1:num_nonzeros) {
		predval <- t(row_pred[,nonzero_rows[j] ]) %*% col_pred[, nonzero_cols[j] ]
		err <- nonzero_entries[j] - predval
		sq_err <- sq_err + err*err
	}
	rmse <- sqrt(sq_err/num_nonzeros)
	rmse
}

# Divide factor combine
dfc <- function(mat, sc, slices) {
	# pick a random permutation of the columns
	cols <- dim(mat)[2]
	#sampleCols <- sample(cols)
	
	# make the matrix into a list of sparse columns
	listMat <- lapply(1:cols, function(i) mat[,i,drop=FALSE])
	
	# distribute the column slices with spark
	# might need to pass in desired num slices here?
	
	subMatRDD <- parallelize(sc,listMat,slices)
	
	# factor each slice
	factorsRDD <- lapplyPartition(subMatRDD,factorCols)
	
	# collect the results and project them onto the first column slice
	factorList <- collect(factorsRDD)

	# collect the results and project them onto the first column slice
	result <- dfcProject(factorList)

	# get the error
	error <- errorCal( mat, result[[1]], result[[2]])
	error
}

sc <- sparkR.init(args[[1]], "DFCR")

slices <- ifelse(length(args) > 1, as.integer(args[[2]]),2)
cat("\n\n",slices,"\n\n")


# Test matrix init
dims <- 100
r <- 10
testU <- matrix(rnorm(r*dims,mean = 4,sd = 1),r,dims)
testV <- matrix(rnorm(r*dims,mean = 4,sd = 1),r,dims)
testG <- matrix(rnorm(dims*dims,mean = 0,sd = .05),dims,dims)
testM <- t(testV) %*% testU + testG

zeros <- sample.int(dims^2, size=floor(0.7*dims^2)) 

zeroedM <- testM
for(i in zeros){
	x = i %% dims
	y = i %/% dims
	zeroedM[x,y] <- 0
}

# Read matrix from file
maskedFile <- args[[3]]
#trueUFile <- args[[4]]
#trueVFile <- args[[5]]

maskedM <- readMM(maskedFile)
#trueU <- read(trueUFile)
#trueV <- read(trueVFile)

error <- dfc(maskedM, sc, slices)
cat("RMSE for the entire matrix: ",error,"\n")
cat("Average magnitude of entries of M: ",sum(abs(maskedM))/(.3*dims^2),"\n")
