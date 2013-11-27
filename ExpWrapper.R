require(SparkR)
require(MASS)
require('Matrix')
require(Rcpp)
require(svd)

upperArgs <- commandArgs(trailing = TRUE)

if (length(upperArgs) < 1) {
	print("Usage: ExpWrapper.R <master> <lower_slices> <upper_slices> <masked_file> <lower_iterations> <upper_iterations> <output_file>")
	q("no")
}
outfile <- paste(upperArgs[[7]],".out",sep="")
master <- paste(upperArgs[[1]],"[",upperArgs[[3]],"]",sep="")
data <- c()
for(i in seq(as.numeric(upperArgs[[2]]),as.numeric(upperArgs[[3]]))) { # num slices
	for(j in seq(as.numeric(upperArgs[[5]]),as.numeric(upperArgs[[6]]),by = 5)) { # num iterations
		ins <- list(master,i,upperArgs[[4]],j)
		commandArgs <- function(trailing) ins
		outs <- c(i,j,unlist(source("DFC.R")))
		data <- rbind(data, outs)
		print(data)
	}
}
write.matrix(data, file=outfile, sep="\t")
		
		
