# TODO: Add comment
# 
# Author: lgoff
###############################################################################

JSdist<-function(mat,...){
	res<-matrix(0,ncol=dim(mat)[2],nrow=dim(mat)[2])
	
#	col_js <- matrix(0,ncol=dim(mat)[2],nrow=1)
#	for(i in 1:dim(mat)[2]){
#	    col_js[,i] <- shannon.entropy(mat[,i])
#   }
	col_js<-apply(mat,MARGIN=2,shannon.entropy)
    #print(col_js)
	colnames(res)<-colnames(mat)
	rownames(res)<-colnames(mat)
	for(i in 1:dim(mat)[2]){
		for(j in i:dim(mat)[2]){
			a<-mat[,i]
			b<-mat[,j]
			JSdiv<-shannon.entropy((a+b)/2)-(col_js[i]+col_js[j])*0.5
			res[i,j] = sqrt(JSdiv)
			res[j,i] = sqrt(JSdiv)
		}
	}
	res<-as.dist(res,...)
	attr(res,"method")<-"JSdist"
	res
}

JSdistVec<-function(p,q){
	JSdiv<-shannon.entropy((p+q)/2)-(shannon.entropy(p)+shannon.entropy(q))*0.5
	JSdist<-sqrt(JSdiv)
	JSdist
}

JSdivVec<-function(p,q){
	JSdiv<-shannon.entropy((p+q)/2)-(shannon.entropy(p)+shannon.entropy(q))*0.5
	#JSdist<-sqrt(JSdiv)
	JSdiv
}

JSdistFromP<-function(mat,q){
	#row_js<-apply(mat,MARGIN=1,shannon.entropy)
	res<-apply(mat,MARGIN=1,function(p) {
				JSdistVec(p,q)
			}
	)
	res
}

makeprobsvec<-function(p){
	phat<-p/sum(p)
	phat[is.na(phat)] = 0
	phat
}

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <=0)
		return(Inf)
	p.norm<-p[p>0]/sum(p)
	-sum( log2(p.norm)*p.norm)
}

makeprobs<-function(a){
	colSums<-apply(a,2,sum)
	b<-t(t(a)/colSums)
	b[is.na(b)] = 0
	b
}

.plotmatrix <- function (data, hexbin=FALSE, mapping = aes())
#Modified from original ggplot2 plotmatrix
{
	grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
	grid <- subset(grid, x != y)
	all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
						xcol <- grid[i, "x"]
						ycol <- grid[i, "y"]
						data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
								x = data[, xcol], y = data[, ycol], data)
					}))
	all$xvar <- factor(all$xvar, levels = names(data))
	all$yvar <- factor(all$yvar, levels = names(data))
	densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
						data.frame(xvar = names(data)[i], yvar = names(data)[i], 
								x = data[, i])
					}))
	mapping <- plyr::defaults(mapping, aes_string(x = "x", y = "y"))
	class(mapping) <- "uneval"
	p <-ggplot(all) + facet_grid(xvar ~ yvar)#, scales = "free")
	
	if(hexbin){ 
		p<- p + geom_hex(mapping,size=1.5,na.rm = TRUE) 
	}else{
		p<- p + geom_point(mapping,alpha=0.2,size=0.8,na.rm=TRUE)
	}
	
	p<- p + stat_density(aes(x = x, 
					y = ..scaled.. * diff(range(x)) + min(x)), data = densities, 
			position = "identity", colour = "grey20", geom = "line")
	
	p
}

.hclustToJSON<-function(d,...){
	require(rjson)
	d$call<-NULL
	res<-toJSON(d,...)
	res
}


.dfToJSONarray <- function(dtf){
	clnms <- colnames(dtf)
	
	name.value <- function(i){
		quote <- '';
		if(class(dtf[, i])!='numeric'){
			quote <- '"';
		}
		
		paste('"', i, '" : ', quote, dtf[,i], quote, sep='')
	}
	
	objs <- apply(sapply(clnms, name.value), 1, function(x){paste(x, collapse=', ')})
	objs <- paste('{', objs, '}')
	
	res <- paste('[', paste(objs, collapse=', '), ']')
	
	return(res)
}

.specificity<-function(mat,logMode=T,pseudocount=1,relative=FALSE,...){
	if(logMode){
		mat<-log2(mat+pseudocount)
	}
	mat<-t(makeprobs(t(mat)))
	d<-diag(ncol(mat))
	res<-apply(d,MARGIN=1,function(q){
				JSdistFromP(mat,q)
			})
	colnames(res)<-paste(colnames(mat))
	
	if(relative){
		res<-res/max(res)
	}
	1-res
}


#multiplot <- function(..., plotlist=NULL, cols) {
#	require(grid)
#	
#	# Make a list from the ... arguments and plotlist
#	plots <- c(list(...), plotlist)
#	
#	numPlots = length(plots)
#	
#	# Make the panel
#	plotCols = cols                          # Number of columns of plots
#	plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
#	
#	# Set up the page
#	grid.newpage()
#	pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
#	vplayout <- function(x, y)
#		viewport(layout.pos.row = x, layout.pos.col = y)
#	
#	# Make each plot, in the correct location
#	for (i in 1:numPlots) {
#		curRow = ceiling(i/plotCols)
#		curCol = (i-1) %% plotCols + 1
#		print(plots[[i]], vp = vplayout(curRow, curCol ))
#	}
#	
#}

#THIS IS NOT MINE....I MUST REMOVE IT PRIOR TO SUBMISSION (For detailed GO analysis, check out clusterProfiler and goProfiles)
#ClusterProfiles <- function(geneClusters, onto="CC", level=3, orgPackage="org.Hs.eg.db") {
#	require(goProfiles)
#	require(plyr)
#	require(ggplot2)
#	clusterProfile <- llply(geneClusters, as.data.frame(basicProfile), onto=onto, level=level, orgPackage = orgPackage)
#	clusterProfile.df <- ldply(clusterProfile, rbind)
#	colnames(clusterProfile.df) <- c("Cluster", "Description", "GOID", "Frequency")
#	clusterProfile.df <- clusterProfile.df[clusterProfile.df$Frequency !=0,]
#	clusterProfile.df$Description <- as.character(clusterProfile.df$Description) ## un-factor
#	clusterProfile.df <- ddply(clusterProfile.df, .(Description), transform, Percent = Frequency/sum(Frequency), Total = sum(Frequency))
#	
#	x <- mdply(clusterProfile.df[, c("Description", "Total")], paste, sep=" (")
#	y <- sapply(x[,3], paste, ")", sep="")
#	clusterProfile.df$Description <- y		### label GO Description with gene counts.
#	clusterProfile.df <-  clusterProfile.df[, -6] ###drop the *Total* column##
#	mtitle <- paste(onto, "Ontology Distribution", sep = " ")
#	p <- ggplot(clusterProfile.df, aes(x = Cluster, y = Description, size = Percent))
#	p <- p + geom_point(colour="steelblue") + opts(title = mtitle) + xlab("") + ylab("")
#	p <- p + opts(axis.text.x = theme_text(colour="black", size="11", vjust = 1))
#	p <- p + opts(axis.text.y = theme_text(colour="black", size="11", hjust = 1))
#	result <- list(data=clusterProfile.df, p=p)
#	return(result)
#}