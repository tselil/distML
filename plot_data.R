# PARAMS ###########################################
size <- 2000
masked <- 90
slices <- c(1,2,4,6,8,10)
parts <- c(2,3,4,5,6,7,8,9)
test <- c(10)
II <- length(slices) # Slices range
JJ <- length(parts) # Trials range
plot_colors <- c("red","blue","green","purple","black","orange")
######################################################
file <- paste(size,masked,sep="")
#filename <- paste("results/movielens/movielens10M_part1",file,sep="")
title <- paste("Average RMSE vs. Time for Gaussian Matrices over ",length(parts)," Trials",sep="")
#title2 <- "RMSE vs. Money for a 2000 x 2000 MC Problem"
outputfile <- paste("2000_90_1_to_9_RvT_graph.pdf",sep="")

outputfile2 <- paste(file,"RvD_graph.pdf",sep="")

# filenames[[I]][[J]] is the filename of data for slices[I] on trials[J]
#filenames <- lapply(seq(1,length(slices)), function(i) lapply(seq(1,length(parts)), function(j) paste("results/movielens/movielens10M_part",parts[j],".mm.slices",slices[i],".results.out",sep="")))
filenames <- lapply(seq(1,length(slices)), function(i) lapply(seq(1,length(parts)), function(j) paste("results/outputfile200090_",parts[j],"_masked.out.slices",slices[i],".results.out",sep="")))


pdf(outputfile)
x <- list()
y <- list()
for(i in seq(1,II)) {
	# x[[i]] and y[[i]] are lists such that x[[i]][[j]] contains the data for the jth run with i slices
	data <- lapply(seq(1,JJ), function(j) read.table(filenames[[i]][[j]],header=T,sep="\t"))
	x[[i]] <- lapply(data, function(d) (d$"subproblem.time"+d$"projection.time"))
	y[[i]] <- lapply(data, function(d) d$"RMSE")
}

xav <- lapply(seq(1,II), function(i) c(0,0,0,0,0))
yav <- lapply(seq(1,II), function(i) c(0,0,0,0,0))
xsd <- lapply(seq(1,II), function(i) c(0,0,0,0,0))
ysd <- lapply(seq(1,II), function(i) c(0,0,0,0,0))

# Average the data for each slice. xav is a list such that xav[[i]] is the average data for i slices
# xsd is a list such that xsd[[i]] is a vector with the SD of xav[[i]]
for(i in seq(1,II)) {
	for(j in seq(1,JJ)) {
		xav[[i]] <- xav[[i]] + x[[i]][[j]]/JJ
		yav[[i]] <- yav[[i]] + y[[i]][[j]]/JJ
		
	}
	for(j in seq(1,JJ)) {
		xsd[[i]] <- xsd[[i]] + (x[[i]][[j]]-xav[[i]])^2
		ysd[[i]] <- ysd[[i]] + (y[[i]][[j]]-yav[[i]])^2
	}
	xsd[[i]] <- sqrt(xsd[[i]])
	ysd[[i]] <- sqrt(ysd[[i]])
}
print(xav)
print(yav)
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

legendcap = unlist(lapply(seq(1,II), function(i) paste(slices[i]," Slice",sep="")))
legend(200,.25,legendcap,cex=.8,col=plot_colors,pch=21,lty=1)
dev.off()

# Second Graph: Optimizer and Test Data
pdf(outputfile2)
x <- list()
y <- list()







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
