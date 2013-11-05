#require(SparkR)
#
#args <- commandArgs(trailing = TRUE)
#
#if (length(args) < 1) {
#	print("Usage: pi <master> [<slices>]")
#	q("no")
#}
#
#sc <- sparkR.init(args[[1]], "DFCR")

#slices <- ifelse(length(args) > 1, as.integer(args[[2]]),2)

dims <- 30
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

cat(testU,"\n\n")
cat(testV,"\n\n")
cat(testM,"\n\n")

sgdBase <- function(mat) {
	# Set Parameters
	m <- dim(mat)[1]
	n <- dim(mat)[2]
	lrate <- .004 # learning rate
	k <- .04 # parameter used to minimize over-fitting
	min_impr <- .0001 # min improvement
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
			cat("root mean squared error: ",rmse)
			cat("\n")
			impr <- rmse_prev - rmse
			t <- t + 1
		}
	}
	list(row_feats,col_feats)
}

feats <- sgdBase(zeroedM)
row_feats <- feats[[1]]
col_feats <- feats[[2]]
dif <- testM - (t(row_feats) %*% col_feats)
cat("dif: ",dif)
cat("\n")
