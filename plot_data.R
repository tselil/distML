filename <- "results/outputfile200090_2_masked.out"
slices <- seq(2,10,by=2)

data <- lapply(seq(1:length(slices)), function(i) read.table(paste(filename,".slices",slices[i],".results.out",sep=""),header=T,sep="\t"))
x <- lapply(data, function(d) d$"subproblem.time"+d$"projection.time")
y <- lapply(data, function(d) d$"RMSE")
plot_colors <- c("red","blue","green","purple","black")
#png(filename="figure.png",height = 295,width=300,bg="white")
plot(x[[1]],y[[1]],type="o",col=plot_colors[1])

for(i in seq(1,length(x))) {
	print(x[[i]])
	lines(x[[i]],y[[i]],type="o",col=plot_colors[i])
}

#dev.off()
