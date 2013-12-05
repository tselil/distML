# PARAMS ###########################################
size <- 2000
masked <- 90
slices <- c(1,2,4,6,8,10)
trials <- c(2,3,4,5,6,7,8)
plot_colors <- c("red","blue","green","purple","black","orange")
######################################################

file <- paste(size,masked,sep="")
filename <- paste("results/outputfile",file,sep="")
title <- paste("Average RMSE vs. Time for a ",size," x ",size," MC Problem over ",length(trials)," Trials",sep="")
#title2 <- "RMSE vs. Money for a 2000 x 2000 MC Problem"
outputfile <- paste(file,"RvT_graph.pdf",sep="")

outputfile2 <- paste(file,"RvD_graph.pdf",sep="")

pdf(outputfile)
x <- list()
y <- list()
for(i in seq(1,length(slices))) {
	# x[[i]] and y[[i]] are lists such that x[[i]][[j]] contains the data for the jth run with i slices
	data <- lapply(seq(1:length(trials)), function(t) read.table(paste(filename,"_",trials[t],"_masked.out.slices",slices[i],".results.out",sep=""),header=T,sep="\t"))
	x[[i]] <- lapply(data, function(d) (d$"subproblem.time"+d$"projection.time"))
	y[[i]] <- lapply(data, function(d) d$"RMSE")
}
xav <- as.list(rep(0,length(slices)))
yav <- as.list(rep(0,length(slices)))
xsd <- as.list(rep(0,length(slices)))
ysd <- as.list(rep(0,length(slices)))

# Average the data for each slice. xav is a list such that xav[[i]] is the average data for i slices
# xsd is a list such that xsd[[i]] is a vector with the SD of xav[[i]]
for(i in seq(1,length(slices))) {
	for(j in seq(1,length(trials))) {
		xav[[i]] <- xav[[i]] + x[[i]][[j]]/length(trials)
		yav[[i]] <- yav[[i]] + y[[i]][[j]]/length(trials)
	}
	for(j in seq(1,length(trials))) {
		xsd[[i]] <- xsd[[i]] + (x[[i]][[j]]-xav[[i]])^2
		ysd[[i]] <- ysd[[i]] + (y[[i]][[j]]-yav[[i]])^2
	}
	xsd[[i]] <- sqrt(xsd[[i]])
	ysd[[i]] <- sqrt(ysd[[i]])
}

#png(filename="figure.png",height = 295,width=300,bg="white")
plot(xav[[1]],yav[[1]],type="o",col=plot_colors[1],xlab="Time",ylab="RMSE")

#errbar(xav[[1]],yav[[1]],yplus=ysd[[1]],yminus=ysd[[1]])
title(main=title, col.main="black",font.main=4)

for(i in seq(1,length(x))) {
	lines(xav[[i]],yav[[i]],type="o",col=plot_colors[i])
	#Error bars
	segments(xav[[i]], yav[[i]]-ysd[[i]],xav[[i]], yav[[i]]+ysd[[i]],lwd=2)
	epsilon = 2
	segments(xav[[i]]-epsilon,yav[[i]]-ysd[[i]],xav[[i]]+epsilon,yav[[i]]-ysd[[i]],lwd=2)
	segments(xav[[i]]-epsilon,yav[[i]]+ysd[[i]],xav[[i]]+epsilon,yav[[i]]+ysd[[i]],lwd=2)
}
legend(250,.25,c("1 Slice","2 Slices", "4 Slices", "6 Slices", "8 Slices", "10 Slices"),cex=.8,col=plot_colors,pch=21,lty=1)
dev.off()

#pdf(outputfile2)
#
#x <- lapply(data, function(d) (d$"subproblem.time" + d$"projection.time")*d$"slices")
#plot(x[[1]],y[[1]],type="o",col=plot_colors[1],xlab="Money",ylab="RMSE",xlim=c(0,1000))
#title(main=title2,col.main="black",font.main=4)
#
#for(i in seq(1,length(x))) {
#	lines(x[[i]],y[[i]],type="o",col=plot_colors[i])
#}
#legend(800,.25,c("1 Slice","2 Slices", "4 Slices", "6 Slices", "8 Slices", "10 Slices"),cex=.8,col=plot_colors,pch=21,lty=1)
#dev.off()
