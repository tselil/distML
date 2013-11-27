require(SparkR)
require(MASS)
require('Matrix')
require(Rcpp)
require(svd)

args <- commandArgs(trailing = TRUE)

if (length(args) < 1) {
	print("Usage: DFC.R <master> <lower_slices> <upper_slices> <masked_file> <lower_iterations> <upper_iterations> <output_file>")
	q("no")
}
outfile <- paste(args[[7]],".out",sep="")
data <- list()
for(i in seq(as.numeric(args[[2]]),as.numeric(args[[3]]))) { # num slices
	for(j in seq(as.numeric(args[[5]]),as.numeric(args[[6]]),by = 5)) { # num iterations
		master <- paste(args[[1]],"[",i,"]",sep="")
		commandArgs <- function(trailing) list(master,i,args[[4]],j)
		outs <- list(i,j,unlist(source("DFC.R")))
		data <- rbind(data, outs)
	}
}
write(data, file=outfile)
		
		
