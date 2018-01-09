library(monocle)
library(reshape2)
library(ggplot2)
library(marray)
library(stringr)
library(heatmap.plus)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(Rtsne)
require(Heatplus)
require(Hmisc)
library(RColorBrewer)
library(gplots)
library(MASS)
library(vegan)
library(slackr)
library(gridExtra)
library(tsne)
library(lfa)
library(jackstraw)
library(ROCR)
library(SC3)
library(ggbiplot)


############
# Color Palette
############
celltype_colors<-c("darkgoldenrod1","blue4")

############
# Helper functions
############

hclust2 <- function(x, method="ward", ...)
  hclust(x, method="ward.D", ...)

dist2s <- function(x, ...)
  as.dist(1-cor(t(x), method="spearman"))

dist2p <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

dist4e <- function(x, ...)
  dist(1-x, method='euclidean');

dist4m <- function(x, ...)
  dist(1-x, method='manhattan');


myHeatCols<-maPalette(low="steelblue",mid="white",high="darkred",k=100)

myCoolerCols<-maPalette(low="black",mid="darkgreen",high="white",k=100)
#myCoolerCols<-maPalette(low="white",mid="red",high="darkred",k=100)

standardize <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

lookupGeneId<-function(eset,gene_names){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  res <- c(res,rownames(fData(eset))[rownames(fData(eset)) %in% gene_names])
  res <- unique(res)
  res
}

lookupGeneName<-function(eset,gene_id){
  res <- fData(eset[gene_id,])$gene_short_name
  #res <- unique(res)
  res
}

#source('hNSC_funcs.R')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# myHeatmap<-function(cds,geneset,pseudotime=FALSE){
#   sub<-cds[lookupGeneId(cds,geneset),]
#   if(pseudotime){
#     heatmap.2(log10(exprs(sub[,order(pData(sub)$Pseudotime,decreasing=FALSE)])+1),scale="none",trace="none",col=brewer.pal(9,"GnBu"),ColSideColors=gg_color_hue(5)[pData(sub[,order(pData(sub)$Pseudotime,decreasing=FALSE)])$State],labRow=fData(sub[,order(pData(sub)$Pseudotime,decreasing=FALSE)])$gene_short_name,Colv=FALSE,labCol=FALSE)#,distfun=function(x){JSdist(t(x))})
#   }else {
#     heatmap.2(log10(exprs(sub)+1),scale="none",trace="none",col=brewer.pal(9,"GnBu"),ColSideColors=gg_color_hue(5)[pData(sub)$State],labRow=fData(sub)$gene_short_name,labCol=FALSE)#,distfun=function(x){JSdist(t(x))})
#   }
# }

myHeatmap<-function(cds,geneset,logMode=FALSE){
  sub<-cds[lookupGeneId(cds,geneset),]
  if(logMode){
    heatmap.2(as.matrix(log2(exprs(sub)+1)),scale="none",trace="none",col=myCoolerCols,labRow=fData(sub)$gene_short_name,ColSideColors=gg_color_hue(3)[pData(sub)$region],labCol=FALSE,distfun=dist,hclustfun=hclust2,key = F)
  }else{
    heatmap.2(as.matrix(exprs(sub)),scale="none",trace="none",col=myCoolerCols,labRow=fData(sub)$gene_short_name,ColSideColors=gg_color_hue(3)[pData(sub)$region],labCol=FALSE,distfun=dist,hclustfun=hclust2,key = F)
  }
}

meltCDS<-function(cds,geneset,logMode=F){
  sub<-cds[lookupGeneId(cds,geneset),]
  sub.expr<-as.matrix(exprs(sub))
  if(logMode){
    sub.expr<-log10(sub.expr+1)
  }
  sub.expr.melt<-melt(sub.expr)
  colnames(sub.expr.melt)<-c("gene_id","cell_id","value")
  res<-merge(sub.expr.melt,pData(sub),by.x="cell_id",by.y="cell_id")
  res<-merge(res,fData(sub),by.x="gene_id",by.y="gene_id")
  res
}

myBarMap<-function(cds,geneset,facet_by="kmeans_tSNE_cluster",color_by="factor(kmeans_tSNE_cluster)",cluster="both",showSummary=T,...){
  sub.melt<-meltCDS(cds,geneset,...)
  facet_by_melt<-strsplit(facet_by,"\\+")[[1]]
  sub.melt.summary<-sub.melt %>%
    dplyr::group_by_(.dots=c("gene_short_name",facet_by_melt)) %>%
    dplyr::summarise(mean=mean(value),median=median(value),sd=sd(value),upper_bound=mean+sd,lower_bound=max(mean-sd,0))

  if(cluster %in% c("row","both",T)){
    sub.sum.mat<-sub.melt.summary %>%
      recast(as.formula(paste("gene_short_name ~",facet_by)),measure.var="mean",fun.aggregate=mean)
    sub.sum.hclust<-hclust2(dist(sub.sum.mat[,-1]))
    gene.order.idx<-order.dendrogram(as.dendrogram(sub.sum.hclust))
    gene.order<-sub.sum.mat$gene_short_name[gene.order.idx]
    sub.melt$gene_short_name<-factor(sub.melt$gene_short_name, levels=gene.order)
  }

  if(cluster %in% c("column","both",T)){
    sub.mat<-sub.melt %>%
      recast(as.formula("gene_short_name ~ cell_id"),measure.var="value",fun.aggregate=mean)
    sub.hclust<-hclust2(dist(t(sub.mat[,-1])))
    cell.order.idx<-order.dendrogram(as.dendrogram(sub.hclust))
    cell.order<-colnames(sub.mat[,-1])[cell.order.idx]
    #print(cell.order)
    sub.melt$cell_id<-factor(sub.melt$cell_id,levels=cell.order)
  }

  p<-ggplot(sub.melt)
  p<-p + geom_bar(aes_string(x="cell_id",y="value",fill=color_by,color=color_by),stat="identity")

  if(showSummary){
    p<-p + geom_hline(aes(yintercept=mean),data=sub.melt.summary,size=1.0)
    p<-p + geom_hline(aes(yintercept=upper_bound),data=sub.melt.summary,linetype="dashed")
    p<-p + geom_hline(aes(yintercept=lower_bound),data=sub.melt.summary,linetype="dashed")
  }
  p<-p +
    facet_grid(as.formula(paste("gene_short_name ~", facet_by)),scale="free",space="free_x",labeller=labeller(.default=label_both,gene_short_name=label_value)) +
    theme_bw() + guides(color=FALSE) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text.y = element_text(angle=0,hjust=0),
          strip.background = element_blank(),
          panel.background = element_rect(fill="white"),
          panel.margin = unit(0, "lines"),
          panel.grid = element_blank()
    )
  p
}

diffAUC <- function(x,y) {
  prediction.use=ROCR::prediction(c(x,y),c(rep(1,length(x)),rep(0,length(y))),0:1)
  perf.use=ROCR::performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}

per_gene_AUC<-function(dat1,dat2,geneset){
  myAUC=unlist(lapply(geneset,function(x)diffAUC(as.numeric(dat1[x,]),as.numeric(dat2[x,]))))
  myAUC[is.na(myAUC)]=0
  avg_diff=unlist(lapply(geneset,function(x)(mean(as.numeric(dat1[x,]))-mean(as.numeric(dat2[x,])))))
  toRet=data.frame(cbind(myAUC,avg_diff),row.names=geneset)
  toRet=toRet[rev(order(toRet$myAUC)),]
  return(toRet)
}

ROC_test<-function(cds,cluster_1,cluster_2,geneset,thresh.use=log(2)) {
  gene_ids<-lookupGeneId(cds,geneset)
  dat<-exprs(cds)
  dat.1<-apply(dat[gene_ids,pData(cds)$kmeans_tSNE_cluster %in% cluster_1],1,mean)
  dat.2<-apply(dat[gene_ids,pData(cds)$kmeans_tSNE_cluster %in% cluster_2],1,mean)
  total.diff<-abs(dat.1-dat.2)
  genes.diff<-names(which(total.diff>thresh.use))
  genes.use<-genes.diff[genes.diff%in%rownames(dat)]
  res<-per_gene_AUC(dat[,pData(cds)$kmeans_tSNE_cluster %in% cluster_1],dat[,pData(cds)$kmeans_tSNE_cluster %in% cluster_2],genes.use)
  res<-res[rev(order(abs(res$myAUC-0.5))),]
  res$power=abs(res$myAUC-0.5)*2
  return(res)
}

ROC_test_region<-function(cds,region_1,region_2,geneset,thresh.use=log(2)) {
  gene_ids<-lookupGeneId(cds,geneset)
  dat<-exprs(cds)
  dat.1<-apply(dat[gene_ids,pData(cds)$region %in% region_1],1,mean)
  dat.2<-apply(dat[gene_ids,pData(cds)$region %in% region_2],1,mean)
  total.diff<-abs(dat.1-dat.2)
  genes.diff<-names(which(total.diff>thresh.use))
  genes.use<-genes.diff[genes.diff%in%rownames(dat)]
  res<-per_gene_AUC(dat[,pData(cds)$region %in% region_1],dat[,pData(cds)$region %in% region_2],genes.use)
  res<-res[rev(order(abs(res$myAUC-0.5))),]
  res$power=abs(res$myAUC-0.5)*2
  return(res)
}



where.expressed<-function(cds,geneId){
  cpc<-round(exprs(cds[lookupGeneId(cds,geneId),]))
  cpc<- as.vector(cpc > 1)
  cpc<- unlist(lapply(cpc,FUN=function(x){if(x){1}else{0}}))
  cpc
}

jaccard_sim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  jaccard <- as.data.frame(m)

  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      jaccard[i,j]= length(which(x[,i] & x[,j])) / length(which(x[,i] | x[,j]))
      jaccard[j,i]=jaccard[i,j]
    }
  }
}

library(Matrix)
jaccard_distance <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)

  ## only non-zero values of common
  Aim = A[im]

  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )

  return( as.dist(1 - J ))
}

theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  #panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_blank()
  #axis.text.y = element_blank()
  #axis.title.x = element_blank(),
  #axis.title.y = element_blank()
)

myggheatmap<-function(cds,geneset=c("CD44","ACHE","ETV4","ISL1","ISL2","NEUROG2","NEUROD1","CHAT","CRIM1","DIAPH3","BCL11B","FEZF2","FEZF1","LHX1","CRYM","NTNG1","CDH13","CDH22","SOX5","BCL6","NETO1","PCP4","ITM2A","POU3F1"),logMode=T,rowK=5){
  sub<-cds[lookupGeneId(cds,geneset),]
  if(logMode){
    exprs(sub)<-log10(exprs(sub)+1.0)
  }

  row.hclust<-hclust(dist(exprs(sub)))

  row.cluster<-as.data.frame(cutree(row.hclust,k=rowK))
  colnames(row.cluster)<-c("rowCluster")
  row.order<-order.dendrogram(as.dendrogram(row.hclust))

  sub.melt<-melt(exprs(sub))

  colnames(sub.melt)<-c("gene_id","sample_id","fpkm")

  sub.melt<-merge(sub.melt,pData(sub),by="sample_id")

  sub.melt<-merge(sub.melt,fData(sub)[,c("gene_short_name","gene_id")],by="gene_id")

  sub.melt<-merge(sub.melt,row.cluster,by.x="gene_id",by.y="row.names")

  sub.melt$gene_short_name<-factor(sub.melt$gene_short_name,levels=fData(sub[row.order,])$gene_short_name)

  sub.melt$sample_id <- factor(sub.melt$sample_id,levels=pData(sub)$sample_id[order(pData(sub)$Pseudotime,decreasing=F)])

  p <- ggplot(sub.melt)
  p + geom_tile(aes(x=sample_id,y=gene_short_name,fill=fpkm)) + facet_grid(rowCluster~State,scales="free",space="free") + theme_change + scale_fill_gradient(low="white",high="darkred") + theme(axis.text.y=element_blank())
}

myCorheatmap<-function(cds,logMode=T,method="color",cor.method="pearson",addrect=NULL,order="hclust",hclust.method="ward",...){
  dat<-as.matrix(exprs(cds))
  if(logMode){
    dat<-log2(dat+1)
  }
  dat.cor<-cor(t(dat),method=cor.method)

  rownames(dat.cor)<-lookupGeneName(cds,rownames(dat.cor))
  colnames(dat.cor)<-lookupGeneName(cds,colnames(dat.cor))

  corrplot(dat.cor,method=method,hclust.method=hclust.method,order=order,addrect=addrect,...)
}


PCbiplot <- function(PC, x="PC1", y="PC2", color="black", shape=NULL) {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_point(alpha=.4, size=3)
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_point(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1,color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75,color="red")
  plot <- plot + theme_bw()
  plot
}

myTSNEPlot<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(tmp))
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,size="value")) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts()
    }else{
      p + geom_point(aes_string(color=color_by,size="value",shape=shape_by)) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by)) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by)) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }
  }
}

myTSNEPlotAlpha<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=3){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    #print(geneMeans)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    #print(genes) #
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=F)
    #print(head(genes)) #
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id",sort = F)
    tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)
    #print(head(tmp)) #
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,alpha="value", shape = shape_by),stroke=0,size=cell_size,pch=19) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1), name = "Scaled Expresssion")
    }else{
      p + geom_point(aes_string(color=color_by,alpha="value",shape=shape_by),stroke=0,size=cell_size,pch=19) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1), name = "Scaled Expression")
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by),size=cell_size, shape = shape_by) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }
  }
}

myTSNEPlotAlpha.triangle<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=3){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    #print(geneMeans)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    #print(genes) #
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=F)
    #print(head(genes)) #
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id",sort = F)
    tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)
    #print(head(tmp)) #
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,alpha="value"),shape = 17,size=cell_size,pch=19) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1), name = "Scaled Expression")
    }else{
      p + geom_point(aes_string(color=color_by,alpha="value"),shape = 17,size=cell_size,pch=19) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1), name = "Scaled Expression")
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by),size=cell_size, shape = 17) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }else{
      p + geom_point(aes_string(color=color_by,shape= 17),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts()
    }
  }
}

myTSNEPlotRainbow<-function(cds,red="Bdnf",green="Fos",blue="empty",logMode=T,shape_by=NULL,scaled=FALSE,cell_size=2,discrete=FALSE){
  tmp<-pData(cds)
  markers<-c(red,green,blue)
  genes<-exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)])
  if(logMode){
    genes<-log10(genes+1)
  }
  geneMeans<-rowMax(genes)
  if(scaled){
    genes<-genes/geneMeans
  }
  genes<-t(genes)
  colnames(genes)<-lookupGeneName(cds,colnames(genes))
  if(discrete){
    genes<-genes/rowSums(genes)
    genes[is.na(genes)]<-0
  }else{
    genes[is.na(genes)]<-1
  }
  genes<-as.data.frame(genes)
  #Map to Rgb
  if(blue=="empty"){
    genes$empty<-0
  }
  genes$plotColor<-as.vector(unlist(apply(genes[,c(red,green,blue)],1,function(x){rgb(x[1],x[2],x[3])})))
  genes$plotColor[rowSums(genes[,markers])==0]<-"#FFFFFF"
  genes$cell_id<-rownames(genes)
  #print(dim(genes))
  genes<-merge(tmp,genes,by.x=0,by.y="cell_id")
  #print(head(genes))
  titleString<-paste("Red=",red,": Green=",green,": Blue=",blue,sep="")
  p<-ggplot(genes,aes(x=tSNE1_pos,y=tSNE2_pos))
  if(is.null(shape_by)){
    p + geom_point(fill="white",color="black",stroke=0.25,size=cell_size) +
      geom_point(color=genes$plotColor,stroke=0,size=cell_size) +
      monocle:::monocle_theme_opts() + ggtitle(titleString)
  }else{
    p + geom_point(aes_string(shape=shape_by),fill="white",color="black",stroke=0.25,size=cell_size) +
      geom_point(aes_string(shape=shape_by),stroke=0,color=genes$plotColor,size=cell_size) +
      monocle:::monocle_theme_opts() + ggtitle(titleString)
  }
}

myTSNEPlotRainbow2<-function(cds,red="Bdnf",green="Fos",blue="empty",logMode=T,shape_by=NULL,scaled=FALSE,cell_size=2,discrete=TRUE){
  tmp<-pData(cds)
  markers<-c(red,green,blue)
  genes<-exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)])
  if(logMode){
    genes<-log10(genes+1)
  }
  geneMeans<-rowMax(genes)
  if(scaled){
    genes<-genes/geneMeans
  }
  genes<-t(genes)
  colnames(genes)<-lookupGeneName(cds,colnames(genes))
  if(discrete){
    #genes<-genes/rowSums(genes)
    #genes[is.na(genes)]<-0

    genes<-as.data.frame(genes)
    #Map to Rgb
    if(blue=="empty"){
      genes$empty<-0
    }
    genes$cell_id<-rownames(genes)
    #print(dim(genes))
    genes<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(genes))
    titleString<-paste("Red=",red,": Green=",green,": Blue=",blue,sep="")
    if(is.null(shape_by)){
      p<-ggplot(genes,aes(x=tSNE1_pos,y=tSNE2_pos))
      p + geom_point(fill="white",color = "black",stroke=0.25,size=cell_size,shape=21) +
        geom_point(aes_string(alpha=red),color="red",fill="red",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=green),color="green",fill="green",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=blue),color="blue",fill="blue",stroke=0,size=cell_size) +
        scale_alpha(range=c(0,0.5)) + guides(alpha=FALSE) +
        monocle:::monocle_theme_opts() + ggtitle(titleString)
    }else{
      p<-ggplot(genes,aes_string(x='tSNE1_pos',y='tSNE2_pos',shape=shape_by))
      p + geom_point(aes_string(shape=shape_by),color="white",stroke=0.25,size=cell_size) +
        geom_point(aes_string(alpha=genes$red),color="red",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=genes$green),color="green",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=genes$blue),color="blue",stroke=0,size=cell_size) +
        scale_alpha(range=c(0,0.5)) + guides(alpha=FALSE) +
        monocle:::monocle_theme_opts() + ggtitle(titleString)
    }
  }else{
    genes[is.na(genes)]<-1

    genes<-as.data.frame(genes)
    #Map to Rgb
    if(blue=="empty"){
      genes$empty<-0
    }
    genes$plotColor<-as.vector(unlist(apply(genes[,c(red,green,blue)],1,function(x){rgb(x[1],x[2],x[3])})))
    genes$plotColor[rowSums(genes[,markers])==0]<-"#FFFFFF"
    genes$cell_id<-rownames(genes)
    #print(dim(genes))
    genes<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(genes))
    titleString<-paste("Red=",red,": Green=",green,": Blue=",blue,sep="")
    p<-ggplot(genes,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(fill="white",color="black",stroke=0.25,size=cell_size) +
        geom_point(color=genes$plotColor,stroke=0,size=cell_size) +
        monocle:::monocle_theme_opts() + ggtitle(titleString)
    }else{
      p + geom_point(aes_string(shape=shape_by),fill="white",color="black",stroke=0.25,size=cell_size) +
        geom_point(aes_string(shape=shape_by),stroke=0,color=genes$plotColor,size=cell_size) +
        monocle:::monocle_theme_opts() + ggtitle(titleString)
    }
  }
}

##########
# Geneset import
##########
#Arking.Neuron<-read.delim("http://arkinglab.org/upload/GeneLists/NEURO.upper",header=F,stringsAsFactors=F)$V1
#Arking.Astro<-read.delim("http://arkinglab.org/upload/GeneLists/ASTRO.upper",header=F,stringsAsFactors=F)$V1
#Arking.Synaptic<-read.delim("http://arkinglab.org/upload/GeneLists/SynapticProteins.upper",header=F,stringsAsFactors=F)$V1
CellCycleGenes<-c("ACD","ACTR1A","AHCTF1","AKAP9","ALMS1","ANAPC1","ANAPC10","ANAPC11","ANAPC2","ANAPC4","ANAPC5","ANAPC7","APITD1","ATM","ATR","ATRIP","AURKA","AURKB","AZI1","B9D2","BIRC5","BRCA1","BTRC","BUB1","BUB1B","BUB3","CASC5","CCDC99","CCNA1","CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2","CCNH","CDC14A","CDC16","CDC20","CDC23","CDC25A","CDC25B","CDC25C","CDC26","CDC26P1","CDC27","CDC45","CDC6","CDC7","CDCA8","CDK1","CDK2","CDK4","CDK5RAP2","CDK6","CDK7","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDT1","CENPA","CENPC1","CENPH","CENPI","CENPJ","CENPK","CENPL","CENPM","CENPN","CENPO","CENPP","CENPQ","CENPT","CEP135","CEP164","CEP192","CEP250","CEP290","CEP41","CEP57","CEP63","CEP70","CEP72","CEP76","CETN2","CHEK1","CHEK2","CKAP5","CKS1B","CLASP1","CLIP1","CNTRL","CSNK1D","CSNK1E","CUL1","DBF4","DCTN1","DCTN2","DCTN3","DHFR","DHFRP1","DIDO1","DKC1","DNA2","DSN1","DYNC1H1","DYNC1I2","DYNLL1","DYRK1A","E2F1","E2F2","E2F3","E2F4","E2F5","ERCC6L","FBXO5","FEN1","FGFR1OP","FKBP6","GINS1","GINS2","GINS4","GMNN","GORASP1","H2AFX","H2AFZ","HAUS2","HDAC1","HIST1H2AB","HIST1H2AC","HIST1H2AD","HIST1H2AE","HIST1H2AJ","HIST1H2BA","HIST1H2BB","HIST1H2BC","HIST1H2BD","HIST1H2BE","HIST1H2BF","HIST1H2BG","HIST1H2BH","HIST1H2BI","HIST1H2BJ","HIST1H2BK","HIST1H2BL","HIST1H2BM","HIST1H2BN","HIST1H2BO","HIST1H4A","HIST1H4B","HIST1H4C","HIST1H4D","HIST1H4E","HIST1H4F","HIST1H4H","HIST1H4I","HIST1H4J","HIST1H4K","HIST1H4L","HIST2H2AA3","HIST2H2AA4","HIST2H2AC","HIST2H2BE","HIST2H4A","HIST2H4B","HIST3H2BB","HIST3H3","HIST4H4","HJURP","HSP90AA1","HSPA2","HUS1","INCENP","ITGB3BP","KIF18A","KIF20A","KIF23","KIF2A","KIF2B","KIF2C","KNTC1","LIG1","LIN37","LIN52","LIN54","LIN9","LMNA","LMNB1","LOC440577","LOC440917","LOC441488","LOC645084","LOC647654","LOC648152","LOC649620","LOC650621","LOC651610","LOC651763","LOC651921","LOC652826","LOC729964","LOC730418","LOC730594","MAD1L1","MAD2L1","MAPRE1","MAX","MCM10","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM8","MDM2","MIS12","MIS18A","MIS18BP1","MLF1IP","MNAT1","MYBL2","MYC","NDC80","NDEL1","NEDD1","NEK2","NHP2","NINL","NPM1","NSL1","NUDC","NUF2","NUMA1","NUP107","NUP133","NUP37","NUP43","NUP85","OFD1","OIP5","ORC1","ORC2","ORC3","ORC4","ORC5","ORC6","PAFAH1B1","PCM1","PCNA","PCNT","PKMYT1","PLK1","PLK4","PMF1","POLA1","POLA2","POLD1","POLD2","POLD3","POLD4","POLE","POLE2","POT1","PPP1CC","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B","PPP2R2A","PPP2R3B","PPP2R5A","PPP2R5B","PPP2R5C","PPP2R5D","PPP2R5E","PRIM1","PRIM2","PRKACA","PRKAR2B","PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMA8","PSMB1","PSMB10","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","PSMD10","PSMD11","PSMD12","PSMD13","PSMD14","PSMD2","PSMD3","PSMD4","PSMD5","PSMD6","PSMD7","PSMD8","PSMD9","PSME1","PSME2","PSME4","PSMF1","PTTG1","RAD1","RAD17","RAD21","RAD9A","RANBP2","RANGAP1","RB1","RBBP4","RBBP4P1","RBBP7","RBL1","RBL2","RCC2","REC8","RFC2","RFC3","RFC4","RFC5","RFWD2","RFWD2P1","RPA1","RPA2","RPA3","RPA4","RPS27","RPS27A","RPS27AP11","RRM2","RSF1","RUVBL1","RUVBL2","SDCCAG8","SEC13","SEH1L","SGOL1","SGOL2","SKA1","SKA2","SKA2L","SKP1","SKP2","SMARCA5","SMC1A","SMC1B","SMC3","SPC24","SPC25","SSNA1","STAG1","STAG2","STAG3","SUN2","SYCP1","SYCP2","SYCP3","SYNE1","SYNE2","TAOK1","TERF1","TERF2","TERF2IP","TERT","TEX12","TFDP1","TINF2","TK2","TP53","TUBA1A","TUBA4A","TUBB","TUBB4A","TUBB4B","TUBBP2","TUBG1","TUBG2","TUBGCP2","TUBGCP3","TUBGCP5","TUBGCP6","TYMS","UBA52","UBE2C","UBE2D1","UBE2E1","UBE2I","WEE1","WRAP53","XPO1","YWHAE","YWHAG","ZW10","ZWILCH","ZWINT")
NSC.genes<-c("NEUROG2","NES","MSI1","SOX2","ASCL1","SOX1","SOX9",
             "FABP7","MSI2","CD133","FGFR4","GLUT1","NEUROD1",
             "RC2","RC1","PAX6","TBR2","EOMES","DCX","BLBP",
             "NCAM","MAP2","TUBB3","FZD9")
TCA.Cycle.genes<-c("ACO2","ADHFE1","BSG","CS","D2HGDH","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","L2HGDH","LDHA","LDHB","LOC283398","LOC646675","LOC646677","LOC650667","LOC650674","LOC650883","LOC651820","MDH2","NNT","OGDH","PDHA1","PDHB","PDHX","PDK1","PDK2","PDK3","PDK4","PDP1","PDP2","PDPR","SDHA","SDHB","SDHC","SDHD","SLC16A1","SLC16A3","SLC16A8","SUCLA2","SUCLA2P1","SUCLG1","SUCLG2")

Neurotransmitter.genes<-c("ABAT","ACHE","ANXA9","BRS3","TSPO","CCKAR","CCKBR","CHAT","CHRM1","CHRM2","CHRM3","CHRNA1","CHRNA2","CHRNA3","CHRNA4","CHRNA5","CHRNA6","CHRNA7","CHRNB1","CHRNB2","CHRNB4","CHRND","CHRNE","CHRNG","COMT","DRD1","DRD2","DRD3","GABRA1","GABRA2","GABRA3","GABRB1","GABRB2","GABRD","GABRE","GABRG1","GABRG2","GABRP","GABRQ","GABRR1","GABRR2","GAD1","GALR1","GALR2","GALR3","GCH1","GCHFR","GLRA1","GLRA2","GLRA3","QRFPR","NPFFR1","MCHR1","PROKR1","PROKR2","NPFFR2","GPR83","GRIA1","GRIN1","GRPR","HCRTR2","HTR1B","HTR2A","HTR3A","HTR3B","MAOA","NMBR","NMUR1","NMUR2","NPY1R","NPY2R","PHOX2A","NPY4R","PRLHR","SLC5A7","SORCS1","SORCS2","SSTR1","SSTR2","SSTR3","SSTR4","TACR1","TACR2","TPH1","B2M","HPRT1","RPL13A","GAPDH","ACTB","HGDC","RTC","RTC","RTC","PPC","PPC","PPC")

Reactome.Apoptosis.genes<-c("ACIN1","ADD1","AKT1","APAF1","APC","APPL1","ARHGAP10","BAD","BAK1","BAX","BBC3","BCAP31","BCL2","BCL2L1","BCL2L11","BID","BIRC2","BMF","BMX","CASP10","CASP3","CASP6","CASP7","CASP8","CASP9","CDH1","CFLAR","CTNNB1","CYCS","DAPK1","DAPK2","DAPK3","DBNL","DCC","DFFA","DFFB","DIABLO","DNM1L","DSG1","DSG2","DSG3","DSP","DYNLL1","DYNLL2","E2F1","FADD","FAS","FASLG","FNTA","GAS2","GSN","GZMB","H1F0","HIST1H1A","HIST1H1B","HIST1H1C","HIST1H1D","HIST1H1E","HMGB1","HMGB2","KPNA1","KPNB1","LMNA","LMNB1","LOC441488","LOC647859","LOC652460","LOC652826","MAGED1","MAPK8","MAPT","MST4","NMT1","OCLN","PAK2","PKP1","PLEC","PMAIP1","PPP3R1","PRKCD","PRKCQ","PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMA8","PSMB1","PSMB10","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","PSMD10","PSMD11","PSMD12","PSMD13","PSMD14","PSMD2","PSMD3","PSMD4","PSMD5","PSMD6","PSMD7","PSMD8","PSMD9","PSME1","PSME2","PSME4","PSMF1","PTK2","RIPK1","ROCK1","ROCK1P1","RPS27A","RPS27AP11","SATB1","SPTAN1","STK24","TFDP1","TJP1","TJP2","TNF","TNFRSF10B","TNFRSF1A","TNFSF10","TP53","TRADD","TRAF2","UBA52","UNC5A","UNC5B","VIM","XIAP","YWHAB")

GO.Neurogenesis.genes<-c("AGRN","ALS2","AMIGO1","APOE","ARTN","ATP2B2","AZU1","BAI1","BAIAP2","BRSK2","BTG4","CDK5","CDK5R1","CDK6","CIT","CLN5","CNTN4","CYFIP1","DPYSL5","DTX1","EIF2B1","EIF2B2","EIF2B3","EIF2B4","EIF2B5","FARP2","FEZ1","FEZ2","GDNF","GHRL","GLI2","KAL1","KCNIP2","KLK8","KRT2","LAMB1","LDB1","LMX1B","LRRC4C","LST1","MAP1S","MAPT","MDGA1","MDGA2","NF1","NF2","NLGN1","NPTN","NRCAM","NRP1","NRP2","NRTN","NRXN1","NRXN3","NTNG1","NTNG2","OPHN1","OTX2","PARD3","PARD6B","PAX2","PCSK9","PICK1","POU4F1","POU6F2","PPT1","RACGAP1","RND1","ROBO1","ROBO2","RTN1","RTN4","RTN4RL1","RTN4RL2","S100B","SEMA3B","SEMA4F","SERPINF1","SHH","SIAH1","SLIT1","SLIT2","SMARCA1","SOD1","SPON2","TGFB2","THY1","TRAPPC4","UBB","UNC5C","VWC2","YWHAG","YWHAH")

ESC.genes<-c("FOXD3",
             "GATA6",
             "GBX2",
             "NANOG",
             "NR5A2",
             "NR6A1",
             "POU5F1",
             "SOX2",
             "TFCP2L1",
             "UTF1",
             "ZFP42COMMD3",
             "CRABP2",
             "EDNRB",
             "FGF4",
             "FGF5",
             "GABRB3",
             "GAL",
             "GRB7",
             "HCK",
             "IFITM1",
             "IL6ST",
             "KIT",
             "LEFTY1",
             "LEFTY2",
             "LIFR",
             "NODAL",
             "NOG",
             "NUMB",
             "PTEN",
             "SFRP2",
             "TDGF1FGF4",
             "FGF5",
             "GDF3",
             "LEFTY1",
             "LEFTY2",
             "NODAL",
             "TDGF1BRIX1",
             "CD9",
             "DIAPH2",
             "DNMT3B",
             "IFITM2",
             "IGF2BP2",
             "LIN28A",
             "PODXL",
             "REST",
             "SEMA3A",
             "TERTFOXA2",
             "GATA4",
             "PTF1ACDX2",
             "EOMES",
             "GCM1",
             "KRT1AFP",
             "SERPINA1FN1",
             "LAMA1",
             "LAMB1",
             "LAMC1",
             "SOX17T",
             "WT1DES",
             "MYF5",
             "MYOD1HBB",
             "HBZCOL1A1",
             "RUNX2NES",
             "NEUROD1",
             "PAX6CD34",
             "CDH5",
             "FLT1",
             "PECAM1DDX4",
             "SYCP3GCG",
             "IAPP",
             "INS",
             "PAX4",
             "PDX1",
             "SSTOLIG2",
             "TA")

Gage_activity_padj_lt_0.01<-c("Arc","Pim1","Ptgs2","Nr4a2","Neat1","Arl4d","Dclk1","Plce1","Baz1a","Fosb","Vgf","Plk2","Atf3","Gpr3","Samd4","Zfp516","Pcdh8","Sgk1","Dpysl5","Siglec1","Lphn3","Midn","Tll1","Phf21b","B230319C09Rik","Egr4","Tsc22d2","Prdm5","Lingo1","Csrnp1","Pcsk1","Ap2b1","Nptx2","Gm25411","Inhba","Baiap2","Bdnf","Spry2","Sertad1","Frmd6","Nr4a3","Nr4a1","Kdm6b","Zdbf2","Blnk","Ddx51","Rnd3","Dmxl1","Calcoco1","Zfp157","Rasl11a","Jdp2","Fbxo33","Homer1","Tacc1","Tmco1","Sik1","Cenpa","Pola2","Fos","Bicd1","Ccrn4l","Rheb","Ankrd33b","Ccdc50","Rfwd2","Nek11","Cstf2","Wdr90","Syne1","Klhdc2","Synpo","Slc25a3","Gm10800","Npas4","Etv3","3110057O12Rik","Cdan1","Hspbap1","Armc9","Xndc1","Gm26870","Gm17024","Dnajc1","Tex29","1700016P03Rik","A730017C20Rik","Cltc","Arsi","Nedd9","Kcnmb4os1","Cep162","Snord16a","Gm12227","Telo2","Sft2d3","Gm3617","Gm26772","Pvrl3","5530601H04Rik","AF357399","Epha5","Tdp1","Rbbp7","Zfand6","Epha7","Ssh2","Clasp1","Prex1","Iqce","Abhd2","Dusp6","B3galt5","Mir673","Snora7a","Bcl7a","Inpp1","Prim2","Mir22hg","Slc6a6","Rasal1","Mir181a-2","Gm28229","Ksr2","Derl2","Hnrnpll","Gm29491","Dusp1","Hltf","Pde4a","Pmepa1","Mir1954","Tmem167","Epha10","Gm12216","Myc","Trip10","Gm17056","Tmem251","Usp45","Ell2","Fmnl1","Amigo3","Dusp5","Dot1l","Rpa1","4931415C17Rik","Glmp","Gm8741","Brip1os","Lmbr1","Ext2","Gm13677","Fam83d","Gm24224","Tet3","Gm13524","Rln3","Gm23202","Ralgps2","Golga5","Exoc1","Capn3","Gm14654","9330159M07Rik","Taok2","Ranbp2","Gm23744","Arid5a","Blvrb","Aldh9a1","Snord12","Mbnl2","Snord65","Zfyve26","Otud3","DQ267101","Scarna3a","Gm27747","Gprc5a","Gpatch2","Oscp1","Csrnp3","Tor1a","Gtf2h2","Shprh","Fhad1os2","Malat1","Tspan33","Gm26767","Blmh","Shc4","1110007C09Rik","Fbxo41","Gem","Gm3076","Tanc1","Mt1","5430405H02Rik","Topbp1","Pld4","Snord49a","Fsd1","Arhgap12","Tas1r3","Fndc3a","Gm16436","Gm6345","Gm10717","Rnpep","Osbpl6","Mir3072","B4galt3","Rgs2","Lmo7","Gadd45g","Gm22358","Adal","Tmem134","Cdkl1","Chst10","Gm6061","Smtn","Gm14239","Gm9917","Armcx1","Gemin8","Tiam1","Fgfr1","Slc26a2","Gm17087","Bai1","Serinc2","Tiam2","Plk3","RP23-459L15.7","Jakmip3","Mcl1","Timm10b","Mthfd2","Ccdc97","Glt8d2","Fam13a","Fam57a","Med14","Bri3","Slc39a3","Hemk1","Gm15201","Pgap1","RP23-382I12.1","Egr3","Kcnj2","Pdk1","Klhl3","Egfl7","Gm28289","Upf1","RP23-14P23.9","Arhgap20","Lrrcc1","Dcaf8","S100a10","Prkab2","Fscn2","Htr1a","Kdm5a","Gm4991","Herc2","Tgfbrap1","Alg6","Acad8","Spag4","Hmga1","Gm26703","Uck1","Klhl5","Farsa","Chrd","Rapgef6","Fam76b","Tmem2","Sds","Tubgcp6","Ncan","Adamts17","Mapk4","Tsen15","Cgref1","Ict1","Paip2b","Clec10a","Ctsz","A830018L16Rik","Gabpb2","Acbd6","Igsf9b","Slc2a1","Map7","Gadd45b","4932438A13Rik","A630023P12Rik","Pkn1","Gm11611","Rbm15","Pbx1","Ift20","Asap2","Ide","Nav1","P4htm","Mapkapk2","Igf2bp2","Brix1","Gm26775","Nub1","Cenpm","Ick","Stxbp4","1500004A13Rik","Orc5","Atg4c","Ptprt","Gm13468","Nek6","Rprl2","Irf3","Gm5601","Mtrr","Mobp","Slc6a17","Gm11794","Dixdc1","Usp30","Stmn4","Gm26530","Fbxl22","Trip12","Rgs7bp","Pam","Flot1","Gm11954","Hspb8","Pms2","Rnu3b4","Rnu3b2","1700102P08Rik","Gpr61","Tarbp2","Npy","Uchl5","Slc6a8","Drc1","Cnot11","Scn3b","Siah2","Gm10801","Scg2","Gm10705","Polg2","Tnfrsf23","Acvr1b","1700019D03Rik","Ifnar1","Kif3c","Cep78","Junb","Cc2d1a","Zfp414","Chgb","Ccnt2","Tyro3","Gm26683","Rell2","Slc16a1","Alpk1","Gm28294","RP24-458F14.4","Skil","Ppm1f","Setdb1","Abca9","Gm20938","Mtfmt","Ankrd32","Dyrk3","Zfp239","Abr","Gprc5c","Lrrc58","D630045J12Rik","RP23-459L15.8","Rrp1b","Stab2","Sbf2","Fbxw7","Gldc","Mysm1","Gpr82","Dirc2","Ptprn","Gls2","Gm15796","Sema3e","Slc6a18","Map3k5","Gm13684","Zfp354c","Zfp831","Ralbp1","Lrrc24","C030039L03Rik","1700064H15Rik","Rasal3","Gm26617","1700001D01Rik","Schip1","Bloc1s2","Nit1","Serpina10","Tmem101","Pak6","A830012C17Rik","Bend4","Gm17229","Kazn","Dcun1d2","Zfp804b","Tsen2")

Gage_FOS_up<-c("Arc","Pim1","Ptgs2","Nr4a2","Neat1","Arl4d","Dclk1","Plce1","Baz1a","Fosb","Vgf","Plk2","Atf3","Gpr3","Samd4","Zfp516","Pcdh8","Sgk1","Dpysl5","Siglec1","Lphn3","Midn","Tll1","Phf21b","B230319C09Rik","Egr4","Tsc22d2","Lingo1","Csrnp1","Pcsk1","Ap2b1","Nptx2","Gm25411","Inhba","Baiap2","Bdnf","Spry2","Sertad1","Frmd6","Nr4a3","Nr4a1","Kdm6b","Zdbf2","Blnk","Ddx51","Rnd3","Dmxl1","Rasl11a","Jdp2","Fbxo33","Homer1","Tacc1","Sik1","Cenpa","Pola2","Fos","Ccrn4l","Rheb","Ankrd33b","Rfwd2","Nek11","Cstf2","Wdr90","Syne1","Synpo","Slc25a3","Gm10800","Npas4","Etv3","Hspbap1","Gm26870","Gm17024","Dnajc1","Tex29","1700016P03Rik","Cltc","Arsi","Nedd9","Kcnmb4os1","Snord16a","Gm12227","Gm3617","Gm26772","Pvrl3","AF357399","Epha5","Rbbp7","Epha7","Ssh2","Prex1","Abhd2","Dusp6","Mir673","Snora7a","Prim2","Mir22hg","Slc6a6","Mir181a-2","Gm28229","Hnrnpll","Gm29491","Dusp1","Pde4a","Pmepa1","Mir1954","Epha10","Myc","Trip10","Gm17056","Tmem251","Ell2","Fmnl1","Amigo3","Dusp5","Dot1l","Rpa1","4931415C17Rik","Gm8741","Brip1os","Gm13677","Fam83d","Gm24224","Tet3","Gm13524","Rln3","Gm23202","Gm14654","9330159M07Rik","Ranbp2","Gm23744","Arid5a","Blvrb","Snord12","Mbnl2","Snord65","Otud3","DQ267101","Scarna3a","Gm27747","Gprc5a","Fhad1os2","Malat1","Gm26767","Blmh","Shc4","1110007C09Rik","Gem","Tanc1","Mt1","Pld4","Snord49a","Tas1r3","Gm16436","Gm6345","Gm10717","Osbpl6","Mir3072","Rgs2","Lmo7","Gadd45g","Gm22358","Gm6061","Smtn","Gm14239","Gm9917","Gemin8","Tiam1","Fgfr1","Gm17087","Bai1","Serinc2","Plk3","RP23-459L15.7","Mcl1","Mthfd2","Med14","Gm15201","Pgap1","RP23-382I12.1","Egr3","Kcnj2","Gm28289","Upf1","S100a10","Prkab2","Fscn2","Kdm5a","Gm4991","Herc2","Spag4","Hmga1","Gm26703","Rapgef6","Tmem2","Sds","Ncan","Mapk4","Cgref1","Clec10a","Ctsz","Igsf9b","Slc2a1","Gadd45b","4932438A13Rik","A630023P12Rik","Gm11611","Rbm15","Asap2","Mapkapk2","Igf2bp2","Brix1","Gm26775","Cenpm","Ptprt","Gm13468","Rprl2","Gm5601","Mobp","Slc6a17","Gm11794","Stmn4","Gm26530","Fbxl22","Trip12","Rgs7bp","Pam","Hspb8","Rnu3b4","Rnu3b2","1700102P08Rik","Npy","Uchl5","Slc6a8","Cnot11","Siah2","Gm10801","Scg2","Polg2","Tnfrsf23","Junb","Chgb","Ccnt2","Tyro3","Gm26683","Rell2","Slc16a1","Alpk1","Gm28294","RP24-458F14.4","Skil","Abca9","Gm20938","Dyrk3","Abr","Lrrc58","RP23-459L15.8","Rrp1b","Stab2","Fbxw7","Gldc","Gpr82","Dirc2","Ptprn","Gls2","Gm15796","Sema3e","Slc6a18","Map3k5","Gm13684","1700064H15Rik","Rasal3","Gm26617","1700001D01Rik","Schip1","Serpina10","Pak6","A830012C17Rik","Bend4","Gm17229","Zfp804b")
Gage_FOS_dn<-c("Prdm5","Calcoco1","Zfp157","Tmco1","Bicd1","Ccdc50","Klhdc2","3110057O12Rik","Cdan1","Armc9","Xndc1","A730017C20Rik","Cep162","Telo2","Sft2d3","5530601H04Rik","Tdp1","Zfand6","Clasp1","Iqce","B3galt5","Bcl7a","Inpp1","Rasal1","Ksr2","Derl2","Hltf","Tmem167","Gm12216","Usp45","Glmp","Lmbr1","Ext2","Ralgps2","Golga5","Exoc1","Capn3","Taok2","Aldh9a1","Zfyve26","Gpatch2","Oscp1","Csrnp3","Tor1a","Gtf2h2","Shprh","Tspan33","Fbxo41","Gm3076","5430405H02Rik","Topbp1","Fsd1","Arhgap12","Fndc3a","Rnpep","B4galt3","Adal","Tmem134","Cdkl1","Chst10","Armcx1","Slc26a2","Tiam2","Jakmip3","Timm10b","Ccdc97","Glt8d2","Fam13a","Fam57a","Bri3","Slc39a3","Hemk1","Pdk1","Klhl3","Egfl7","RP23-14P23.9","Arhgap20","Lrrcc1","Dcaf8","Htr1a","Tgfbrap1","Alg6","Acad8","Uck1","Klhl5","Farsa","Chrd","Fam76b","Tubgcp6","Adamts17","Tsen15","Ict1","Paip2b","A830018L16Rik","Gabpb2","Acbd6","Map7","Pkn1","Pbx1","Ift20","Ide","Nav1","P4htm","Nub1","Ick","Stxbp4","1500004A13Rik","Orc5","Atg4c","Nek6","Irf3","Mtrr","Dixdc1","Usp30","Flot1","Gm11954","Pms2","Gpr61","Tarbp2","Drc1","Scn3b","Gm10705","Acvr1b","1700019D03Rik","Ifnar1","Kif3c","Cep78","Cc2d1a","Zfp414","Ppm1f","Setdb1","Mtfmt","Ankrd32","Zfp239","Gprc5c","D630045J12Rik","Sbf2","Mysm1","Zfp354c","Zfp831","Ralbp1","Lrrc24","C030039L03Rik","Bloc1s2","Nit1","Tmem101","Kazn","Dcun1d2","Tsen2")

myBoxplot.subset <- function(cds,markers=NULL,logMode=T,color_by="color", scaled = FALSE){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log2(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)
    #print(head(tmp))

    p<-ggplot(tmp,aes(factor(subset.cluster),value, fill = subset.cluster))
    p + geom_boxplot() + facet_wrap('gene_short_name', scales = "free_y") + theme_bw() + ylab("log2(Transcripts+1)")
  }
}

myBoxplot.cluster <- function(cds,markers=NULL,logMode=T,color_by="color", scaled = FALSE){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log2(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    tmp$gene_short_name <- factor(tmp$gene_short_name, levels = markers)
    #print(head(tmp))

    p<-ggplot(tmp,aes(factor(kmeans_tSNE_cluster),value, fill = kmeans_tSNE_cluster))
    p + geom_boxplot() + facet_wrap('gene_short_name', scales = "free_y") + theme_bw() + ylab("log2(Transcripts+1)")
  }
}

cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix("", mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

