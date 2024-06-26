####################################################################################################################################################################################
###################################################################### code for figure6 ####################################################################################### 
####################################################################################################################################################################################

#############################################################Figure 6A 6D
P240.subset.obj<-subset_opt(nano.obj.with.label.transfer, cells =nano.obj.with.label.transfer@images$P240@boundaries$centroids@cells)
Idents(P240.subset.obj)<-"LT_full_gene"

tumor.fov.id<-c(nano.obj.with.label.transfer@images$P240_fov9@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov10@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov11@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov12@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov13@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov14@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov15@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov16@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov17@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov18@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov19@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov20@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov21@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov24@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov25@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov26@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov27@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov28@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov29@boundaries$centroids@cells,
                nano.obj.with.label.transfer@images$P240_fov31@boundaries$centroids@cells)

TLS.fov.id<-c(nano.obj.with.label.transfer@images$P240_fov1@boundaries$centroids@cells,
              nano.obj.with.label.transfer@images$P240_fov2@boundaries$centroids@cells,
              nano.obj.with.label.transfer@images$P240_fov5@boundaries$centroids@cells,
              nano.obj.with.label.transfer@images$P240_fov6@boundaries$centroids@cells,
              nano.obj.with.label.transfer@images$P240_fov7@boundaries$centroids@cells)

PN.fov.id<-c(nano.obj.with.label.transfer@images$P240_fov3@boundaries$centroids@cells,
             nano.obj.with.label.transfer@images$P240_fov22@boundaries$centroids@cells,
             nano.obj.with.label.transfer@images$P240_fov23@boundaries$centroids@cells)

P240.subset.obj<-subset_opt(P240.subset.obj,cells = c(tumor.fov.id,TLS.fov.id,PN.fov.id))

region<-as.vector(colnames(P240.subset.obj))
region[region%in% tumor.fov.id]<-"Tumor"

region[region%in% TLS.fov.id]<-"TLS"

region[region%in% PN.fov.id]<-"PN"

table(region)
PN   TLS Tumor 
3865  6362 25151 

P240.subset.obj@meta.data$region<-region

P240.Tcell.obj<-subset(P240.subset.obj,idents = c("CD4.T.cells","CD8.T.cells"))

DefaultAssay(P240.Tcell.obj)<-"SCT"
P240.Tcell.obj.normalized=exp(GetAssayData(object = P240.Tcell.obj))-1
P240.Tcell.obj.normalized=round(P240.Tcell.obj.normalized, 0)
dim(P240.Tcell.obj.normalized)
[1]  989 2104

P240.Tcell.obj.cells2cluster=P240.Tcell.obj@meta.data$region
names(P240.Tcell.obj.cells2cluster)=rownames(P240.Tcell.obj@meta.data)
table(P240.Tcell.obj.cells2cluster)
PN   TLS Tumor 
497  1046   561 

P240.Tcell.obj.cells2cluster1=paste0(P240.Tcell.obj@meta.data$region,"_",P240.Tcell.obj@meta.data$LT_full_gene)
names(P240.Tcell.obj.cells2cluster1)=rownames(P240.Tcell.obj@meta.data)
table(P240.Tcell.obj.cells2cluster1)
PN_CD4.T.cells    PN_CD8.T.cells   TLS_CD4.T.cells   TLS_CD8.T.cells Tumor_CD4.T.cells 
406                91               848               198               341 
Tumor_CD8.T.cells 
220 

P240.Tcell.obj.nk=length(unique(P240.Tcell.obj.cells2cluster))
P240.Tcell.obj.cluster.k.order=c("Tumor","PN","TLS")
P240.Tcell.obj.cluster.k.order1=c("Tumor_CD4.T.cells","Tumor_CD8.T.cells","PN_CD4.T.cells","PN_CD8.T.cells","TLS_CD4.T.cells","TLS_CD8.T.cells")

P240.Tcell.obj.dat.nk.n= as.numeric(unlist(lapply(P240.Tcell.obj.cluster.k.order, function(x){ sum(P240.Tcell.obj.cells2cluster==x) })))

#set cluster distance and position
t.space=50
P240.Tcell.obj.dat.nk.space=c(0, unlist(lapply(P240.Tcell.obj.dat.nk.n, function(x){c(rep(0, x-1), t.space)})))
P240.Tcell.obj.dat.nk.space=rev(rev(P240.Tcell.obj.dat.nk.space)[-1])
## and center position of each cluster
P240.Tcell.obj.dat.nk.lab.p=cumsum(P240.Tcell.obj.dat.nk.n)-(P240.Tcell.obj.dat.nk.n/2)

P240.Tcell.obj.dat.cell.order= unlist(lapply(P240.Tcell.obj.cluster.k.order1, function(x){
  names(P240.Tcell.obj.cells2cluster1[P240.Tcell.obj.cells2cluster1==x]) }) )

P240.Tcell.obj.norm.ordered=P240.Tcell.obj.normalized[, P240.Tcell.obj.dat.cell.order]
P240.Tcell.obj.norm.ordered=as.matrix(P240.Tcell.obj.norm.ordered)
P240.Tcell.obj.n.c=ncol(P240.Tcell.obj.norm.ordered)



P240.Tcell.obj.color<-c(rep("#08A045",341),rep("orange",220),rep("#08A045",406),rep("orange",91),rep("#08A045",848),rep("orange",198))

#test plot
t.g="CCR7";
#g.title=paste(t.g, " (K=", cef.gene.k[t.g,6], ")", sep = "")
barplot(as.numeric((P240.Tcell.obj.norm.ordered[t.g,])), space=P240.Tcell.obj.dat.nk.space, col=P240.Tcell.obj.color,main=t.g,las=2, border = NA,ylab="counts")
barplot(rep(0, P240.Tcell.obj.n.c), col="grey",  space=P240.Tcell.obj.dat.nk.space,add=TRUE, axes=F)
cluster.name<-c("Tumor","PN","TLS")

for(p in 1:P240.Tcell.obj.nk){   
  t.labels=cluster.name[p]
  axis(1,tcl=-0.2,P240.Tcell.obj.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.5,las=1,tick = F)
}

svg("Figure6AD_P240_Tcell_barplot.svg",8,10)
g.select<-c("CCR7","IL7R")
par(mfrow=c(5,1))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
for(i in 1:length(g.select)){
  t.g=g.select[i];#g.select means the gene list you want to use
  barplot(as.numeric((P240.Tcell.obj.norm.ordered[t.g,])), space=P240.Tcell.obj.dat.nk.space, col=P240.Tcell.obj.color,main=t.g,las=2, border = NA,ylab="Normalised Counts")
  barplot(rep(0, P240.Tcell.obj.n.c), col="grey",  space=P240.Tcell.obj.dat.nk.space,add=TRUE, axes=F)
  #axis(1, tcl=-0.2,-10, paste("Cluster#", "Cell#",sep="\n"),  padj = 0.1, col = "white",cex.axis=1)
  for(p in 1:P240.Tcell.obj.nk){   
    t.labels=cluster.name[p]
    axis(1, tcl=1,P240.Tcell.obj.dat.nk.lab.p[p]+t.space*(p-1),t.labels, padj = -2,  cex.axis=0.7,hadj = 0.5,las=1,tick = F)
  }
}  

dev.off()

#############################################################Figure 6B 6C
P240.fov22.id<-as.vector(c(P240.subset.obj@images$P240_fov22@boundaries$centroids@cells))
P240.fov22.obj<-subset_opt(P240.subset.obj, cells =P240.fov22.id)

Idents(P240.fov22.obj)<-"LT_full_gene"

P240.fov22.meta.data<-data.frame(P240.fov22.obj@meta.data)
P240.fov22.meta.data.f<-P240.fov22.meta.data[,c(10,11)]
P240.fov22.meta.data.f$celltype<-as.vector(Idents(P240.fov22.obj))
P240.fov22.meta.data.f<-cbind(paste0("X",rownames(P240.fov22.meta.data.f)),P240.fov22.meta.data.f)
colnames(P240.fov22.meta.data.f)<-c("cell","x","y","celltype")

P240.fov22.counts<-data.frame(P240.fov22.obj@assays$SCT@data)

# create SpaTalk data
P240.fov22.SpaT.obj <- createSpaTalk(st_data = as.matrix(P240.fov22.counts),
                                     st_meta =P240.fov22.meta.data.f[,-4],
                                     species = "Human",
                                     if_st_is_sc = T,
                                     spot_max_cell = 1,
                                     celltype = P240.fov22.meta.data.f$celltype)

# Filter LRIs with downstream targets
P240.fov22.SpaT.obj <- find_lr_path(object = P240.fov22.SpaT.obj, lrpairs = lrpairs, pathways = pathways)

# Infer cell-cell communications from all celltype pairs
P240.fov22.SpaT.obj <- dec_cci_all(object = P240.fov22.SpaT.obj,use_n_cores = 10)
save(P240.fov22.SpaT.obj,file = "P240_fov22_cell-cell_interatcion.RData")

P240.fov22.LR_pair<-P240.fov22.SpaT.obj@lrpair
write.xlsx(P240.fov22.LR_pair,file = "P240_fov22_cell-cell_interaction.xlsx",rowNames=T,colNames=T)

svg("Figure6B_P240_fov22_CCR7_interaction.svg",width = 7,height = 4)
plot_lrpair(object = obj,
            celltype_sender = 'Endothelial',
            ligand = 'CCL19',
            celltype_receiver = 'CD4.T.cells',
            receptor = 'CCR7',
            if_plot_density = F,
            size = 1.5,
            color = c("red","#08A045","lightgray"),
            arrow_length = 0.05)
dev.off()

svg("Figure6C_P240_fov22_CCR7_interaction.svg",width = 7,height = 4)
plot_lrpair(object = obj,
            celltype_sender = 'Myeloid_cell',
            ligand = 'CCL19',
            celltype_receiver = 'CD4.T.cells',
            receptor = 'CCR7',
            if_plot_density = F,
            size = 1.5,
            color = c("#AC844B","#08A045","lightgray"),
            arrow_length = 0.05)
dev.off()


#############################################################Figure 6E 6F 6J 6K

#load interaction anlysis results of P240_fov5 from Figure5 code
P240.fov5.SpaT.obj <- load("P240_fov5_cell-cell_interatcion.RData")

svg("Figure6E_P240_fov5_IL7(Mye)_interaction.svg",width = 7,height = 4)
plot_lrpair(object = P240.fov5.SpaT.obj,
            celltype_sender = 'Myeloid_cell',
            ligand = 'IL7',
            celltype_receiver = 'CD4.T.cells',
            receptor = 'IL7R',
            if_plot_density = F,
            size = 1.5,
            color = c("#AC844B","#08A045","lightgray"),
            arrow_length = 0.05)
dev.off()

svg("Figure6F_P240_fov5_IL7(FB)_interaction.svg",width = 7,height = 4)
plot_lrpair(object = P240.fov5.SpaT.obj,
            celltype_sender = 'Fibroblast',
            ligand = 'IL7',
            celltype_receiver = 'CD4.T.cells',
            receptor = 'IL7R',
            if_plot_density = F,
            size = 1.5,
            color = c("#619CFF","#08A045","lightgray"),
            arrow_length = 0.05)
dev.off()

svg("Figure6J_P240_fov5_CXCL12(Mye)_CXCR4(CD4)_interaction.svg",width = 7,height = 4)
plot_lrpair(object = P240.fov5.SpaT.obj,
            celltype_sender = 'Myeloid_cell',
            ligand = 'CXCL12',
            celltype_receiver = 'CD4.T.cells',
            receptor = 'CXCR4',
            if_plot_density = F,
            size = 1.5,
            color = c("#AC844B","#08A045","lightgray"),
            arrow_length = 0.05)
dev.off()

svg("Figure6K_P240_fov5_CXCL12(Mye)_CXCR4(B)interaction.svg",width = 7,height = 4)
plot_lrpair(object = P240.fov5.SpaT.obj,
            celltype_sender = 'Myeloid_cell',
            ligand = 'CXCL12',
            celltype_receiver = 'B_cell',
            receptor = 'CXCR4',
            if_plot_density = F,
            size = 1.5,
            color = c("#AC844B","magenta","lightgray"),
            arrow_length = 0.05)
dev.off()

#############################################################Figure 6G

svg("P240_fov5_IL7(Mye)_IL7R(CD4)_downstream.svg",width = 7,height = 4)
plot_lr_path(object = P240.fov5.SpaT.obj,                
             celltype_sender = 'Myeloid_cell',
             ligand = 'IL7',
             celltype_receiver = 'CD4.T.cells',
             receptor = 'IL7R',
             color = c("#AC844B","#08A045"))+labs(title = "Downstream targets and TFs")
dev.off()

#############################################################Figure 6H

svg("P240_fov5_IL7(Fib)_IL7R(CD4)_downstream.svg",width = 7,height = 4)
plot_lr_path(object = P240.fov5.SpaT.obj,                
             celltype_sender = 'Fibroblast',
             ligand = 'IL7',
             celltype_receiver = 'CD4.T.cells',
             receptor = 'IL7R',
             color = c("white","#08A045"))+labs(title = "Downstream targets and TFs")
dev.off()

#############################################################Figure 6I
P240.fov5.id<-as.vector(nano.obj@images$P240_fov5@boundaries$centroids@cells)
P237.fov8.id<-as.vector(nano.obj@images$P237_fov8@boundaries$centroids@cells)
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

svg("Figure6K_2FOV_VCAM1_barplot.svg",8,10)
par(mfrow=c(5,1))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
g.select<-"VCAM1"
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
P240     60         148        362          445    398
P237   1107          26         99          593    296

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

svg("Figure6K_2FOV_VCAM1_barplot_color_by_sample.svg",8,10)
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

