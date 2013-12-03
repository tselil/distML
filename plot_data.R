file <- "200090_8"
filename <- paste("results/outputfile",file,"_masked.out",sep="")
title <- "RMSE vs. Time for a 2000 x 2000 MC Problem"
title2 <- "RMSE vs. Money for a 2000 x 2000 MC Problem"
outputfile <- paste(file,"RvT_graph.pdf",sep="")
outputfile2 <- paste(file,"RvD_graph.pdf",sep="")
slices <- seq(2,10,by=2)

pdf(outputfile)
data <- lapply(seq(1:length(slices)), function(i) read.table(paste(filename,".slices",slices[i],".results.out",sep=""),header=T,sep="\t"))
x <- lapply(data, function(d) d$"subproblem.time"+d$"projection.time")
y <- lapply(data, function(d) d$"RMSE")
plot_colors <- c("red","blue","green","purple","black")
#png(filename="figure.png",height = 295,width=300,bg="white")
plot(x[[1]],y[[1]],type="o",col=plot_colors[1],xlab="Time",ylab="RMSE")
title(main=title, col.main="black",font.main=4)

for(i in seq(1,length(x))) {
	lines(x[[i]],y[[i]],type="o",col=plot_colors[i])
}
dev.off()

pdf(outputfile2)
x <- lapply(data, function(d) (d$"subproblem.time" + d$"projection.time")*d$"slices")
plot(x[[1]],y[[1]],type="o",col=plot_colors[1],xlab="Money",ylab="RMSE",xlim=c(0,1000))
title(main=title2,col.main="black",font.main=4)

for(i in seq(1,length(x))) {
	lines(x[[i]],y[[i]],type="o",col=plot_colors[i])
}
dev.off()
