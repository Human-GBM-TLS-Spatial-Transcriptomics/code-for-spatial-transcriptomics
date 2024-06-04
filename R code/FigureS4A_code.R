####################################################################################################################################################################################
###################################################################### code for figureS4 ####################################################################################### 
####################################################################################################################################################################################

#Figure S4A
svg("FigureS4A_UMAP_aglined_pc30_r0.6.svg",width = 8,height = 6)
DimPlot(nano.obj.with.label.transfer,group.by = "SCT_snn_res.0.6",label=T,raster=FALSE,pt.size = 0.1)
dev.off()



####gene signature score
library(UCell)
markers<-list()
markers$Endothelial<-c("PECAM1","VWF","CDH5","CD34","KDR",'ANGPT2')
markers$Bcell<-c("CD79A","MS4A1")
markers$Fibroblast<-c("FN1","DCN","LUM")
markers$Mural<-c("ACTA2","PDGFRB","RGS5","NOTCH3")
markers$Neutrophil<-c("S100A8","S100A9","CSF3R")
markers$Tcell<-c("CD2","CD3E")
markers$Tumor<-c("GFAP","SOX2","SOX9")
markers$Myeloid<-c("C1QA","C1QB","CSF1R")
markers$erythrocyte<-c("HBB","HBA1","CD36")



marker_score<-AddModuleScore_UCell(nano.obj.with.label.transfer,features = markers)


library(stringr)
library(ggplot2)
library(viridis)

#Figure S4B-S4J
svg("Endothelial_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Endothelial_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Endothelial\nSignature\nScore"))
dev.off()

svg("Bcell_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Bcell_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("B cell\nSignature\nScore"))
dev.off()

svg("Mural_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Mural_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Mural cell\nSignature\nScore"))
dev.off()

svg("Tumor_gene_signature_UMAP(new).svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Tumor_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Tumor cell\nSignature\nScore"))
dev.off()

svg("Neutrophil_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Neutrophil_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Neutrophil\nSignature\nScore"))
dev.off()

svg("Fibroblast_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Fibroblast_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Fibroblast\nSignature\nScore"))
dev.off()

svg("Myeloid_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Myeloid_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Myeloid cell\nSignature\nScore"))
dev.off()

svg("Tcell_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "Tcell_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("T cell\nSignature\nScore"))
dev.off()

svg("Erythrocyte_gene_signature_UMAP.svg",width = 8,height = 6)
FeaturePlot(marker_score,features = "erythrocyte_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Erythrocyte\nSignature\nScore"))
dev.off()

#Figure S4K
svg("FigureS4K_CXCL12_expression_in_situ.svg",width = 18,height = 6)
ImageFeaturePlot(nano.obj.with.label.transfer, fov = c("STGG1_fov5","STGU1_fov8"), 
                 cells = WhichCells(nano.obj,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                 features = c("CXCL12"), min.cutoff = "q5",max.cutoff = "q95")
dev.off()

#Figure S4L
STGG1.fov5.id<-as.vector(nano.obj@images$STGG1_fov5@boundaries$centroids@cells)
STGU1.fov8.id<-as.vector(nano.obj@images$STGU1_fov8@boundaries$centroids@cells)
barplot.obj<-subset_opt(nano.obj,cells = c(STGG1.fov5.id,STGU1.fov8.id))
Idents(barplot.obj)<-"celltype.new"
barplot.obj.f<-subset(barplot.obj,idents = c("Endothelial","Fibroblast","Myeloid cell","T cell","B cell"))

DefaultAssay(barplot.obj.f)<-"SCT"
barplot.normalized=exp(GetAssayData(object = barplot.obj.f))-1
barplot.normalized=round(barplot.normalized, 0)
dim(barplot.normalized)
[1]  989 3534

barplot.cells2cluster=barplot.obj.f@meta.data$celltype.new
names(barplot.cells2cluster)=rownames(barplot.obj.f@meta.data)
table(barplot.cells2cluster)
B cell  Endothelial   Fibroblast Myeloid cell       T cell 
1167          174          461         1038          694 

barplot.nk=length(unique(barplot.cells2cluster))
barplot.cluster.k.order=c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")

barplot.dat.nk.n= as.numeric(unlist(lapply(barplot.cluster.k.order, function(x){ sum(barplot.cells2cluster==x) })))

#set cluster distance and position
t.space=20
barplot.dat.nk.space=c(0, unlist(lapply(barplot.dat.nk.n, function(x){c(rep(0, x-1), t.space)})))
barplot.dat.nk.space=rev(rev(barplot.dat.nk.space)[-1])
## and center position of each cluster
barplot.dat.nk.lab.p=cumsum(barplot.dat.nk.n)-(barplot.dat.nk.n/2)

barplot.order.info<-data.frame("Cluster"=barplot.obj.f@meta.data$celltype.new,"Sample"=barplot.obj.f@meta.data$Sample,"id"=rownames(barplot.obj.f@meta.data))
barplot.order.info<-barplot.order.info[order(match(barplot.order.info$Cluster,c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")),barplot.order.info$Sample),]

barplot.dat.cell.order= unlist(lapply(barplot.cluster.k.order, function(x){
  names(barplot.cells2cluster[barplot.cells2cluster==x]) }) )

barplot.norm.ordered=barplot.normalized[,barplot.dat.cell.order]
barplot.norm.ordered=as.matrix(barplot.norm.ordered)
barplot.n.c=ncol(barplot.norm.ordered)

#test barplot
t.g="GPR183";
barplot.color<-c(rep("cyan",694),rep("magenta",1167),rep("red",174),rep("#619CFF",461),rep("#AC844B",1038))
#g.title=paste(t.g, " (K=", cef.gene.k[t.g,6], ")", sep = "")
barplot(as.numeric((barplot.norm.ordered[t.g,])), space=barplot.dat.nk.space, col=barplot.color,main=t.g,las=2, border = NA,ylab="counts")
barplot(rep(0, barplot.n.c), col="grey",  space=barplot.dat.nk.space,add=TRUE, axes=F)
cluster.name<-c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")

for(p in 1:barplot.nk){   
  t.labels=cluster.name[p]
  axis(1, tcl=1,barplot.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.4,las=1,tick = F)
}

svg("FigureS4L_2FOV_VCAM1_barplot.svg",8,10)
par(mfrow=c(5,1))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
g.select<-"CXCL12"
for(i in 1:length(g.select)){
  t.g=g.select[i];#g.select means the gene list you want to use
  barplot(as.numeric((barplot.norm.ordered[t.g,])), space=barplot.dat.nk.space, col=barplot.color,main=t.g,las=2, border = NA,ylab="Normalised Counts")
  barplot(rep(0, barplot.n.c), col="grey",  space=barplot.dat.nk.space,add=TRUE, axes=F)
  #axis(1, tcl=-0.2,-10, paste("Cluster#", "Cell#",sep="\n"),  padj = 0.1, col = "white",cex.axis=1)
  for(p in 1:barplot.nk){   
    t.labels=cluster.name[p]
    axis(1, tcl=1,barplot.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.4,las=1,tick = F)
  }
}  
dev.off()


#barplot color by fov
barplot.order.info<-data.frame("Cluster"=barplot.obj.f@meta.data$celltype.new,"Sample"=barplot.obj.f@meta.data$Sample,"id"=rownames(barplot.obj.f@meta.data))
barplot.order.info<-barplot.order.info[order(match(barplot.order.info$Cluster,c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")),barplot.order.info$Sample),]

barplot.dat.cell.order1= barplot.order.info$id
barplot.norm.ordered1=barplot.normalized[,barplot.dat.cell.order1]
barplot.norm.ordered1=as.matrix(barplot.norm.ordered1)
barplot.n.c1=ncol(barplot.norm.ordered1)


table(barplot.obj.f@meta.data$Sample,barplot.obj.f@meta.data$celltype.new)
B cell Endothelial Fibroblast Myeloid cell T cell
STGG1     60         148        362          445    398
STGU1   1107          26         99          593    296

barplot.color1<-c(rep("red",398),rep("darkblue",296),
                  rep("red",60),rep("darkblue",1107),
                  rep("red",148),rep("darkblue",26),
                  rep("red",362),rep("darkblue",99),
                  rep("red",445),rep("darkblue",593))
#test
#g.title=paste(t.g, " (K=", cef.gene.k[t.g,6], ")", sep = "")
barplot(as.numeric((barplot.norm.ordered1[t.g,])), space=barplot.dat.nk.space, col=barplot.color1,main=t.g,las=2, border = NA,ylab="counts")
barplot(rep(0, barplot.n.c1), col="grey",  space=barplot.dat.nk.space,add=TRUE, axes=F)
cluster.name<-c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")

for(p in 1:barplot.nk){   
  t.labels=cluster.name[p]
  axis(1, tcl=1,barplot.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.4,las=1,tick = F)
}

svg("FigureS4L_2FOV_CXCL12_barplot_color_by_sample.svg",8,10)
par(mfrow=c(5,1))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
for(i in 1:length(g.select)){
  t.g=g.select[i];#g.select means the gene list you want to use
  barplot(as.numeric((barplot.norm.ordered1[t.g,])), space=barplot.dat.nk.space, col=barplot.color1,main=t.g,las=2, border = NA,ylab="Normalised Counts")
  barplot(rep(0, barplot.n.c1), col="grey",  space=barplot.dat.nk.space,add=TRUE, axes=F)
  #axis(1, tcl=-0.2,-10, paste("Cluster#", "Cell#",sep="\n"),  padj = 0.1, col = "white",cex.axis=1)
  for(p in 1:barplot.nk){   
    t.labels=cluster.name[p]
    axis(1, tcl=1,barplot.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.4,las=1,tick = F)
  }
}  
dev.off()