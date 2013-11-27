require(SparkR)
require(MASS)
require('Matrix')
require(Rcpp)
require(svd)

upperArgs <- commandArgs(trailing = TRUE)

if (length(upperArgs) < 1) {
	print("Usage: DFC.R <master> <lower_slices> <upper_slices> <masked_file> <lower_iterations> <upper_iterations> <output_file>")
	q("no")
}
outfile <- paste(upperArgs[[7]],".out",sep="")
data <- list()
for(i in seq(as.numeric(upperArgs[[2]]),as.numeric(upperArgs[[3]]))) { # num slices
	for(j in seq(as.numeric(upperArgs[[5]]),as.numeric(upperArgs[[6]]),by = 5)) { # num iterations
		master <- paste(upperArgs[[1]],"[",i,"]",sep="")
		commandArgs <- function(trailing) list(master,i,upperArgs[[4]],j)
		outs <- list(i,j,unlist(source("DFC.R")))
		data <- rbind(data, outs)
	}
}
write(data, file=outfile)
		
		
