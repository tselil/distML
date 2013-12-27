First, run R. Within the R interpreter, run:

install.packages("Rcpp")
install.packages("svd")

You are now ready to run DFC! Type the following into the command line:



We have included a sample matrix, called "example.mat". 
# Read matrix from file
maskedFile <- args[[3]]
iterations <- args[[4]]
if(length(args) > 5) {
	outfile <- args[[6]]
}
