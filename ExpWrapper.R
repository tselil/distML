require(SparkR)
require(MASS)
require('Matrix')
require(Rcpp)
require(svd)

upperArgs <- commandArgs(trailing = TRUE)

if (length(upperArgs) < 1) {
	print("Usage: ExpWrapper.R <master> <lower_slices> <upper_slices> <masked_file> <lower_iterations> <upper_iterations> <iteration_increment> <output_file>")
	q("no")
}
outfile <- paste(upperArgs[[7]],".out",sep="")
if (upperArgs[[1]] == "local"){
	master <- paste(upperArgs[[1]],"[",upperArgs[[3]],"]",sep="")
} else {
	master <- upperArgs[[1]]
}
data <- c()
for(i in seq(as.numeric(upperArgs[[2]]),as.numeric(upperArgs[[3]]))) { # num slices
	for(j in seq(as.numeric(upperArgs[[5]]),as.numeric(upperArgs[[6]]),by = as.numeric(upperArgs[[7]]))) { # num iterations
		ins <- list(master,i,upperArgs[[4]],j)
		commandArgs <- function(trailing) ins
		outs <- c(i,j,unlist(source("DFC.R")))
		data <- rbind(data, outs)
		colnames(data) <- c("slices","iterations","RMSE","overhead","subproblem time","projection time","visible")
		print(data)
	}
}
write.matrix(data, file=outfile, sep="\t")
		
		
