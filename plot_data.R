# PARAMS ###########################################
size <- 4000
masked <- 90
slices <- c(1,2,4,6,8)

# Parts to Train on, and validate on. They SHOULD NOT overlap
parts <- c(1,2,3,4)
validate <- 5


test <- c(10)
II <- length(slices) # Slices range
JJ <- length(parts) # Trials range
plot_colors <- c("black","black","black","black","black","black")
margins <- c(6,6,6,6)

# Params for Gaussians
xlims <- c(0,700)
ylims <- c(.06,.1)
legendloc <- c(350,.1)

# Params for Movielens
#xlims <- c(0,800)
#ylims <- c(.4,1)
#legendloc <- c(400,1)

######################################################

# filenames[[I]][[J]] is the filename of data for slices[I] on trials[J]

# Params for Gaussians
filenames <- lapply(seq(1,length(slices)), function(i) lapply(seq(1,length(parts)), function(j) paste("results/gaussian/4k/outputfile400090_",parts[j],"_masked.out.slices",slices[i],".results.out",sep="")))
outputfile <- paste("Paper/Graphs/4000_90_1_to_4_RvT_graph.pdf",sep="")
title <- paste("Trained on Average RMSE vs. Time over ",length(parts)," Trials:\nGaussian Random Matrices",sep="")

# Params for Movielens
#filenames <- lapply(seq(1,length(slices)), function(i) lapply(seq(1,length(parts)), function(j) paste("results/movielens/movielens10M_part",parts[j],".mm.slices",slices[i],".results.out",sep="")))
#outputfile <- paste("Paper/Graphs/movielens10M_1_to_4_RvT_graph.pdf",sep="")
#title <- paste("Trained on Average RMSE vs. Time over ",length(parts)," Trials:\nMovielens10M Data",sep="")

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

#png(filename="figure.png",height = 295,width=300,bg="white")
par(bg = "white",mar=margins)
plot_symbols=c(15,16,17,18,19)
plot_line_types = c(2,3,2,2,2)
plot(xav[[1]],yav[[1]],type="o",col=plot_colors[1],xlab="Time",ylab="RMSE",xlim=xlims,ylim=ylims,cex.lab=2,pch=plot_symbols[1],lty=plot_line_types[1])

#errbar(xav[[1]],yav[[1]],yplus=ysd[[1]],yminus=ysd[[1]])
title(main=title, col.main="black",cex.main=1.7)

for(i in seq(1,length(xav))) {
	lines(xav[[i]],yav[[i]],type="o",col=plot_colors[i],pch=plot_symbols[i],lty=plot_line_types[i])
	#Error bars
	segments(xav[[i]], yav[[i]]-ysd[[i]],xav[[i]], yav[[i]]+ysd[[i]],lwd=2)
	epsilon = 2
	segments(xav[[i]]-epsilon,yav[[i]]-ysd[[i]],xav[[i]]+epsilon,yav[[i]]-ysd[[i]],lwd=2)
	segments(xav[[i]]-epsilon,yav[[i]]+ysd[[i]],xav[[i]]+epsilon,yav[[i]]+ysd[[i]],lwd=2)
}

legendcap = unlist(lapply(seq(1,II), function(i) paste(slices[i]," Slice",sep="")))
legend(legendloc[1],legendloc[2],c("1 Slice","2 Slices","4 Slices","6 Slices","8 Slices"),col=plot_colors,pch=plot_symbols,lty=plot_line_types,cex=1.7)
dev.off()

########################################################################################
# Second Graph: Optimizer and Test Data

# Params for Gaussians
opt_title <- "Fixed-parameter and Optimized Time vs. Error:\nGaussian Random Matrices"
opt_output <- "Paper/Graphs/Opt_Gaussian_4000.pdf"
filenames_to_plot <- lapply(seq(1,length(slices)), function(i) paste("results/gaussian/4k/outputfile400090_5_masked.out.slices",slices[i],".results.out",sep=""))
opt_pts <- c(100,125,150,175,200,250,300,400,700)
opt_filenames <- lapply(seq(1,length(opt_pts)), function(i) paste("optResults/4k_no",validate,"_",opt_pts[i],".results",sep=""))

# Params for Movielens
#opt_title <- "Fixed-parameter and Optimized Time vs. Error:\nMovielens10M Data"
#opt_output <- "Paper/Graphs/Opt_movielens.pdf"
#filenames_to_plot <- lapply(seq(1,length(slices)), function(i) paste("results/movielens/movielens10M_part5.mm.slices",slices[i],".results.out",sep=""))
#opt_pts <- c(100,125,150,175,200,250,300,400,700)
#opt_filenames <- lapply(seq(1,length(opt_pts)), function(i) paste("optResults/movie_no",validate,"_",opt_pts[i],".results",sep=""))

pdf(opt_output)
data <- lapply(seq(1,length(filenames_to_plot)), function(i) read.table(filenames_to_plot[[i]],header=T,sep="\t"))

# x[[i]] and y[[i]] contain the data for the ith filename
x <- lapply(data, function(d) d$"subproblem.time" + d$"projection.time")
y <- lapply(data, function(d) d$"RMSE")

opt_data <- lapply(seq(1,length(opt_filenames)), function(i) read.table(opt_filenames[[i]],header=T,sep="\t"))

opt_read <- function(i) {
	opt_data <- read.table(opt_filenames[[i]],header=T,sep="\t")
	c(opt_data$"subproblem.time"+opt_data$"projection.time",opt_data$"RMSE")
}

opt <- lapply(seq(1,length(opt_filenames)), function(i) opt_read(i))
opt_points <- do.call(rbind,opt)
x_opt <- opt_points[,1]
y_opt <- opt_points[,2]
par(bg = "white",mar=margins)
plot(x_opt,y_opt,type="o",col="black",lwd=3,xlab="Time",ylab="RMSE",xlim=xlims,ylim=ylims,cex.lab=2,lty=1,pch=4)

plot_symbols=c(15,16,17,18,19)
plot_line_types = c(2,3,2,2,2)
for(i in seq(1,length(x))) {
	lines(x[[i]],y[[i]],type="o",col=plot_colors[i],pch=plot_symbols[i],lty=plot_line_types[i])
}
legend(legendloc[1],legendloc[2],c("1 Slice","2 Slices","4 Slices","6 Slices","8 Slices","Optimizer"),col=plot_colors,pch=c(15,16,17,18,19,4),lty=c(2,3,2,2,2,1),lwd=c(1,1,1,1,1,3),cex=1.7)
title(opt_title,cex.main=1.7)
dev.off()

