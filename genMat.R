library(MASS)
library('Matrix')


args <- commandArgs(trailing = TRUE)

if (length(args) < 1) {
	print("Usage: genMat.R <dimensions> <rank> <noisyness> <masked_proportion> <output_file_name>")
	q("no")
}

# Test matrix init
dim <- as.integer(args[[1]])
r <- as.integer(args[[2]])
noise <- as.double(args[[3]])
maskProp <- as.double(args[[4]])
trueU <- Matrix(rnorm(r*dim,mean = 0,sd = 1.0/sqrt(r)),r,dim)
trueV <- Matrix(rnorm(r*dim,mean = 0,sd = 1.0/sqrt(r)),r,dim)
# noiseMat <- Matrix(rnorm(dims*dims,mean = 0,sd = noise),dims,dims)
# noisyM <- t(trueV) %*% trueU + noiseMat


unMasked <- sample.int(dim^2, size=floor((1-maskProp)*dim^2)) 

is <-unlist(lapply(unMasked, function(i) ((i-1) %% dim)+1)) 
js <-unlist(lapply(unMasked, function(j) ((j-1)%/% dim)+1))
maskedM <- sparseMatrix(i=is, j=js, dims=c(dim,dim), index1=TRUE, x=1.0)

for(i in 1:length(is)){
	x <- is[[i]]
	y <- js[[i]]
	maskedM[x,y] <- t(trueV[,x]) %*% trueU[,y] + rnorm(1, mean=0, sd=noise)
}

maskedName <- paste(args[[5]], "_masked.out", collapse=" ")
trueNameU <- paste(args[[5]], "_trueU.out", collapse=" ")
trueNameV <- paste(args[[5]], "_trueV.out", collapse=" ")

writeMM(maskedM, file=maskedName)
write.matrix(trueU, file=trueNameU)
write.matrix(trueV, file=trueNameV)
warnings()
