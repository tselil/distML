library(MASS)
library('Matrix')
require(Rcpp)


args <- commandArgs(trailing = TRUE)

if (length(args) < 1) {
	print("Usage: genMat.R <dimensions> <rank> <noisyness> <masked_proportion> <output_file_name>")
	q("no")
}

# Load and compile the fast C++ code
sourceCpp('maskUV.cpp')

# Test matrix init
dim <- as.integer(args[[1]])
r <- as.integer(args[[2]])
noise <- as.double(args[[3]])
maskProp <- as.double(args[[4]])
trueU <- Matrix(rnorm(r*dim,mean = 0,sd = 1.0/sqrt(sqrt(r))),dim,r)
trueV <- Matrix(rnorm(r*dim,mean = 0,sd = 1.0/sqrt(sqrt(r))),dim,r)
# noiseMat <- Matrix(rnorm(dims*dims,mean = 0,sd = noise),dims,dims)
# noisyM <- t(trueV) %*% trueU + noiseMat


unMasked <- sample.int(dim^2, size=floor((1-maskProp)*dim^2)) 

is <-unlist(lapply(unMasked, function(i) ((i-1) %% dim)+1)) 
js <-unlist(lapply(unMasked, function(j) ((j-1)%/% dim)+1))
rand <- rnorm(length(is), mean=0, sd=sqrt(noise))

maskedM <- sparseMatrix(i=is, j=js, dims=c(dim,dim), index1=TRUE, x=rand)
# Call into C++ code to get relevant entries of UV'
maskedUV <- sparseMatrix(i=is,j=js,x=maskUV(as.matrix(trueU),as.matrix(trueV),is,js))
maskedM <- maskedUV + maskedM

maskedName <- paste(args[[5]], "_masked.out", sep= "")
trueNameU <- paste(args[[5]], "_trueU.out", sep = "")
trueNameV <- paste(args[[5]], "_trueV.out", sep = "")

writeMM(maskedM, file=maskedName)
write.matrix(trueU, file=trueNameU)
write.matrix(trueV, file=trueNameV)
warnings()
