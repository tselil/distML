# PARAMS ###########################################
size <- 4000
masked <- 90
slices <- c(1,2,4,6,8)
parts <- c(1,2,3,4,5)
test <- c(10)
II <- length(slices) # Slices range
JJ <- length(parts) # Trials range
plot_colors <- c("red","blue","green","purple","dark blue","black")
######################################################
file <- paste(size,masked,sep="")
#filename <- paste("results/movielens/movielens10M_part1",file,sep="")
title <- paste("Average RMSE vs. Time for Gaussian Matrices over ",length(parts)," Trials",sep="")
#title2 <- "RMSE vs. Money for a 2000 x 2000 MC Problem"
outputfile <- paste("4000_90_1_to_5_RvT_graph.pdf",sep="")

outputfile2 <- paste(file,"RvD_graph.pdf",sep="")

# filenames[[I]][[J]] is the filename of data for slices[I] on trials[J]
filenames <- lapply(seq(1,length(slices)), function(i) lapply(seq(1,length(parts)), function(j) paste("results/gaussian/4k/outputfile400090_",parts[j],"_masked.out.slices",slices[i],".results.out",sep="")))
#filenames <- lapply(seq(1,length(slices)), function(i) lapply(seq(1,length(parts)), function(j) paste("results/outputfile200090_",parts[j],"_masked.out.slices",slices[i],".results.out",sep="")))


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
	if(length(xav[[i]]) != length(x[[i]][[1]])) {
		xav[[i]] = rep(0,length(x[[i]][[1]]))
		xsd[[i]] = rep(0,length(x[[i]][[1]]))
	}
	if(length(yav[[i]]) != length(y[[i]][[1]])) {
		yav[[i]] = rep(0,length(y[[i]][[1]]))
		ysd[[i]] = rep(0,length(y[[i]][[1]]))
	}
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
print(xsd)

#png(filename="figure.png",height = 295,width=300,bg="white")
plot(xav[[1]],yav[[1]],type="o",col=plot_colors[1],xlab="Time",ylab="RMSE",xlim=c(0,1000))

#errbar(xav[[1]],yav[[1]],yplus=ysd[[1]],yminus=ysd[[1]])
title(main=title, col.main="black",font.main=4)

for(i in seq(1,length(xav))) {
	lines(xav[[i]],yav[[i]],type="o",col=plot_colors[i])
	#Error bars
	segments(xav[[i]], yav[[i]]-ysd[[i]],xav[[i]], yav[[i]]+ysd[[i]],lwd=2)
	epsilon = 2
	segments(xav[[i]]-epsilon,yav[[i]]-ysd[[i]],xav[[i]]+epsilon,yav[[i]]-ysd[[i]],lwd=2)
	segments(xav[[i]]-epsilon,yav[[i]]+ysd[[i]],xav[[i]]+epsilon,yav[[i]]+ysd[[i]],lwd=2)
}

legendcap = unlist(lapply(seq(1,II), function(i) paste(slices[i]," Slice",sep="")))
legend(800,.25,legendcap,cex=.8,col=plot_colors,pch=21,lty=1)
dev.off()

# Second Graph: Optimizer and Test Data
opt_title <- "Trained Optimizer and Static Parameters on 4000 x 4000 Gaussian Matrix"
opt_output <- "Opt_Gaussian_4000.pdf"
pdf(opt_output)

filenames_to_plot <-
opt_filenames <- 
data <- lapply(seq(1,length(filenames_to_plot)), function(i) read.table(filenames_to_plot[i],header=T,sep="\t"))

# x[[i]] and y[[i]] contain the data for the ith filename
x <- lapply(data, function(d) d$"subproblem.time" + d$"projection.time")
y <- lapply(data, function(d) d$"RMSE")

opt_data <- lapply(seq(1,length(opt_filenames)), function(i) read.table(opt_filenames[i],header=T,sep="\t"))
print(opt_data)

opt_read <- function(i) {
	opt_data <- read.table(opt_filenames[i],header=T,sep="\t")
	c(opt_data$"subproblem.time"+opt_data$"projection.time",opt_data$"RMSE")
}

opt <- lapply(seq(1,length(opt_filenames)), function(i) opt_read(i))
opt_points <- do.call(rbind,opt)
x_opt <- opt_points[,1]
y_opt <- opt_points[,2]
plot(x[[1]],y[[1]],type="o",col=plot_colors[1],xlab="Time",ylab="RMSE")

for(i in seq(1,length(x))) {
	lines(x[[i]],y[[i]],type="o",col=plot_colors[i])
}
lines(x_opt,y_opt,type="o",col="black",lwd=3)
legend()
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
