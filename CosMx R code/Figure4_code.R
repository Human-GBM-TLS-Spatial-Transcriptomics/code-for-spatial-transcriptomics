####################################################################################################################################################################################
###################################################################### code for figure4 ####################################################################################### 
####################################################################################################################################################################################


#fraction plot theme
fraction_plot_theme <- theme(axis.title=element_text(face="bold", size=10,colour = 'white'), 
                 axis.text=element_text(face="bold", size=7,colour = 'white'), 
                 axis.line = element_line(size=0.5, colour = 'white'), 
                 axis.line.y = element_blank(),axis.line.x = element_line(colour = "white"),  
                 axis.ticks.y = element_blank(), axis.title.y = element_blank(),axis.title.x = element_text(colour = "white"),
                 axis.text.y = element_blank(),axis.ticks.x = element_line(colour = "white"),
                 panel.grid = element_blank(),
                 panel.background = element_rect(fill="black"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank(),
                 plot.background = element_rect(fill = "black"))  

######Figure 4J
svg("Figure4J_UMAP_all_annotated_new_color_white.svg",width = 8,height = 6)
DimPlot(nano.obj.with.label.transfer,raster = F,label = F,cols = c("red","#D39200","#AC844B","#00BA38","#00C19F","cyan","magenta","#DB72FB","#619CFF"))
dev.off()

######Figure 4K
Idents(nano.obj.with.label.transfer)<-"celltype.new"

p1<-ImageDimPlot(nano.obj.with.label.transfer, fov = "P240_fov5",size = 0.5, coord.fixed = FALSE,axes = T,
                 cells = WhichCells(nano.obj.with.label.transfer,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                 cols=c("red","#AC844B","cyan","magenta","white"),border.color = "black")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))


P240.fov5.id<-as.vector(nano.obj.with.label.transfer@images$P240_fov5@boundaries$centroids@cells)
P240.fov5.obj<-subset_opt(nano.obj.with.label.transfer, cells =P240.fov5.id)

DimPlot(P240.fov5.obj)
table(P240.fov5.obj@meta.data$celltype.new)
B cell  Endothelial  Erythrocyte   Fibroblast   Mural cell Myeloid cell   Neutrophil       T cell   Tumor cell 
60          148           83          362          110          445          109          398           57 

P240.fov5.cell.frac<-data.frame(table(P240.fov5.obj@meta.data$celltype.new))
P240.fov5.cell.frac$FOV<-"T-TLS"
P240.fov5.cell.frac<-P240.fov5.cell.frac[-c(3,5,7,9),]

colnames(P240.fov5.cell.frac) = c('Celltype','Number of cell','Group')
P240.fov5.cell.frac<-P240.fov5.cell.frac[order(P240.fov5.cell.frac$`Number of cell`,decreasing = F),]
P240.fov5.cell.frac$Celltype<-factor(P240.fov5.cell.frac$Celltype,levels = P240.fov5.cell.frac$Celltype)


p2<-ggplot(P240.fov5.cell.frac)+geom_col(aes(x=Celltype,y=P240.fov5.cell.frac$`Number of cell`,fill=Celltype),width = 0.8)+
  fraction_plot_theme+coord_flip()+ylab("Number of cells")+NoLegend()+theme(plot.margin = margin(0,0,0,-1))+
  scale_fill_manual(values = c("magenta","red","white","cyan","#AC844B"))

svg("Figure4K_P240_fov5_with_cell_fraction.svg",width = 8,height = 6)
p1+p2+plot_layout(widths = c(5, 1))
dev.off()

######Figure 4L
p3<-ImageDimPlot(nano.obj.with.label.transfer, fov = "P237_fov8",size = 0.5, coord.fixed = FALSE,axes = T,
                 cells = WhichCells(nano.obj.with.label.transfer,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                 cols=c("red","#AC844B","cyan","magenta","white"),border.color = "black")+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))

P237.fov8.id<-as.vector(nano.obj@images$P237_fov8@boundaries$centroids@cells)
P237.fov8.obj<-subset_opt(nano.obj, cells =P237.fov8.id)

DimPlot(P237.fov8.obj)
table(P237.fov8.obj@meta.data$celltype.new)
B cell  Endothelial  Erythrocyte   Fibroblast   Mural cell Myeloid cell   Neutrophil       T cell   Tumor cell 
1107           26            1           99           13          593            3          296          880 

P237.fov8.cell.frac<-data.frame(table(P237.fov8.obj@meta.data$celltype.new))
P237.fov8.cell.frac<-P237.fov8.cell.frac[-c(3,5,7,9),]
P237.fov8.cell.frac$FOV<-"B-TLS"

colnames(P237.fov8.cell.frac) = c('Celltype','Number of cell','Group')
P237.fov8.cell.frac<-P237.fov8.cell.frac[order(P237.fov8.cell.frac$`Number of cell`,decreasing = F),]
P237.fov8.cell.frac$Celltype<-factor(P237.fov8.cell.frac$Celltype,levels = P237.fov8.cell.frac$Celltype)

p4<-ggplot(P237.fov8.cell.frac)+geom_col(aes(x=Celltype,y=P237.fov8.cell.frac$`Number of cell`,fill=Celltype),width = 0.8)+
  fraction_plot_theme+coord_flip()+ylab("Number of cells")+NoLegend()+theme(plot.margin = margin(0,0,0,-1))+
  scale_fill_manual(values = c("red","white","cyan","#AC844B","magenta"))

svg("Figure4L_P237_fov8_with_cell_fraction.svg",width = 8,height = 6)
p3+p4+plot_layout(widths = c(5, 1))
dev.off()

######Figure 4M

svg("Figure4M_IL7R_expression_in_situ.svg",width = 18,height = 6)
ImageFeaturePlot(nano.obj.with.label.transfer, fov = c("P240_fov5","P237_fov8"), 
                 cells = WhichCells(nano.obj,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                 features = c("IL7R"), min.cutoff = "q5",max.cutoff = "q95")
dev.off()

######Figure 4N
barplot.obj<-subset_opt(nano.obj,cells = c(P240.fov5.id,P237.fov8.id))
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

#barplot.order.info<-data.frame("Cluster"=barplot.obj.f@meta.data$celltype.new,"Sample"=barplot.obj.f@meta.data$Sample,"id"=rownames(barplot.obj.f@meta.data))
#barplot.order.info<-barplot.order.info[order(match(barplot.order.info$Cluster,c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")),barplot.order.info$Sample),]

barplot.dat.cell.order= unlist(lapply(barplot.cluster.k.order, function(x){
  names(barplot.cells2cluster[barplot.cells2cluster==x]) }) )

barplot.norm.ordered=barplot.normalized[,barplot.dat.cell.order]
barplot.norm.ordered=as.matrix(barplot.norm.ordered)
barplot.n.c=ncol(barplot.norm.ordered)

barplot.color<-c(rep("cyan",694),rep("magenta",1167),rep("red",174),rep("#619CFF",461),rep("#AC844B",1038))

#test bar plot
t.g="CCL5";
#g.title=paste(t.g, " (K=", cef.gene.k[t.g,6], ")", sep = "")
barplot(as.numeric((barplot.norm.ordered[t.g,])), space=barplot.dat.nk.space, col=barplot.color,main=t.g,las=2, border = NA,ylab="counts")
barplot(rep(0, barplot.n.c), col="grey",  space=barplot.dat.nk.space,add=TRUE, axes=F)
cluster.name<-c("T cell","B cell","Endothelial","Fibroblast","Myeloid cell")

for(p in 1:barplot.nk){   
  t.labels=cluster.name[p]
  axis(1, tcl=1,barplot.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.4,las=1,tick = F)
}


g.select<-c("IL7R")
svg("Figure4N_2FOV_gene_barplot.svg",8,10)
par(mfrow=c(5,1))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
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


