library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(glmGamPoi)
library(dplyr)
library(harmony)

####################################################################################################################################################################################
###################################################################### Seurat QC and clustering #################################################################################### 
####################################################################################################################################################################################

#Load S1
STGG1.obj <- LoadNanostring(
  data.dir = "/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S1",fov = "Run5753_S1"
)
dim(STGG1.obj)
[1]  1008 39587

STGG1.np <- grep("^NegPrb", rownames(STGG1.obj),value = T)

STGG1.obj<-STGG1.obj[setdiff(rownames(STGG1.obj),STGG1.np),]

STGG1.metadata<-read.csv("/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S1/Run5753_S1_metadata_file.csv",header = T,sep = ",")
STGG1.metadata$cell_ID_new<-paste0(STGG1.metadata$cell_ID,"_",STGG1.metadata$fov)
identical(colnames(STGG1.obj),STGG1.metadata$cell_ID_new)

STGG1.obj<-subset(STGG1.obj,cells = STGG1.metadata$cell_ID_new)
dim(STGG1.obj)
[1] 989 39587

STGG1.metadata.f<-subset(STGG1.metadata,STGG1.metadata$cell_ID_new%in%colnames(STGG1.obj))
identical(colnames(STGG1.obj),STGG1.metadata.f$cell_ID_new)
[1] TRUE

STGG1.obj@meta.data[4:24]<-STGG1.metadata.f[,1:21]
STGG1.obj@meta.data[,"Sample"]<-"STGG1"
STGG1.obj@meta.data[,"Sample_name"]<-"Run5753_S1"

#caculate sd of area
STGG1.aera.sd<-sd(STGG1.obj@meta.data$Area)
STGG1.aera.mean<-mean(STGG1.obj@meta.data$Area)
STGG1.area.thresh<-c(STGG1.aera.mean+5*STGG1.aera.sd)

STGG1.double<-subset(STGG1.obj,subset = Area>=as.numeric(STGG1.area.thresh))
STGG1.double.1<-subset(STGG1.obj,subset = nCount_Nanostring>=2000)


pdf("~/Desktop/TLS/CosMx_counts/STGG1_double_QC.pdf",width = 6,height = 4)
VlnPlot(STGG1.obj, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample",pt.size = 0)
VlnPlot(STGG1.double, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
VlnPlot(STGG1.double.1, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
dev.off()

STGG1.f<-subset(STGG1.obj, subset = nCount_Nanostring<2000 & Area<as.numeric(STGG1.area.thresh))
dim(STGG1.f)
[1]   989 39573

save(STGG1.f,file = "STGG1_filtered.RData")

#Load S2
STGU1.obj <- LoadNanostring(
  data.dir = "/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S2",fov = "Run5753_S2"
)

dim(STGU1.obj)
[1]  1008 60254

STGU1.np <- grep("^NegPrb", rownames(STGU1.obj),value = T)

STGU1.obj<-STGU1.obj[setdiff(rownames(STGU1.obj),STGU1.np),]

STGU1.metadata<-read.csv("/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S2/Run5753_S2_metadata_file.csv",header = T,sep = ",")
STGU1.metadata$cell_ID_new<-paste0(STGU1.metadata$cell_ID,"_",STGU1.metadata$fov)
identical(colnames(STGU1.obj),STGU1.metadata$cell_ID_new)

STGU1.obj<-subset(STGU1.obj,cells = STGU1.metadata$cell_ID_new)
dim(STGU1.obj)
[1]   989 60254

STGU1.metadata.f<-subset(STGU1.metadata,STGU1.metadata$cell_ID_new%in%colnames(STGU1.obj))
identical(colnames(STGU1.obj),STGU1.metadata.f$cell_ID_new)
[1] TRUE

STGU1.obj@meta.data[4:24]<-STGU1.metadata.f[,1:21]
STGU1.obj@meta.data[,"Sample"]<-"STGU1"
STGU1.obj@meta.data[,"Sample_name"]<-"Run5753_S2"

#caculate sd of area
STGU1.aera.sd<-sd(STGU1.obj@meta.data$Area)
STGU1.aera.mean<-mean(STGU1.obj@meta.data$Area)
STGU1.area.thresh<-c(STGU1.aera.mean+5*STGU1.aera.sd)

STGU1.double<-subset(STGU1.obj,subset = Area>=as.numeric(STGU1.area.thresh))
STGU1.double.1<-subset(STGU1.obj,subset = nCount_Nanostring>=2000)


pdf("~/Desktop/TLS/CosMx_counts/STGU1_double_QC.pdf",width = 6,height = 4)
VlnPlot(STGU1.obj, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample",pt.size = 0)
VlnPlot(STGU1.double, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
VlnPlot(STGU1.double.1, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
dev.off()

STGU1.f<-subset(STGU1.obj, subset = nCount_Nanostring<2000 & Area<as.numeric(STGU1.area.thresh))
dim(STGU1.f)
[1]   989 60060

save(STGU1.f,file = "STGU1_filtered.RData")

#Load S3
STGU3.obj <- LoadNanostring(
  data.dir = "/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S3",fov = "Run5753_S3"
)
dim(STGU3.obj)
[1]  1008 45056

STGU3.np <- grep("^NegPrb", rownames(STGU3.obj),value = T)

STGU3.obj<-STGU3.obj[setdiff(rownames(STGU3.obj),STGU3.np),]

STGU3.metadata<-read.csv("/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S3/Run5753_S3_metadata_file.csv",header = T,sep = ",")
STGU3.metadata$cell_ID_new<-paste0(STGU3.metadata$cell_ID,"_",STGU3.metadata$fov)
identical(colnames(STGU3.obj),STGU3.metadata$cell_ID_new)

STGU3.obj<-subset(STGU3.obj,cells = STGU3.metadata$cell_ID_new)
dim(STGU3.obj)
[1]   989 45056

STGU3.metadata.f<-subset(STGU3.metadata,STGU3.metadata$cell_ID_new%in%colnames(STGU3.obj))
identical(colnames(STGU3.obj),STGU3.metadata.f$cell_ID_new)
[1] TRUE

STGU3.obj@meta.data[4:24]<-STGU3.metadata.f[,1:21]
STGU3.obj@meta.data[,"Sample"]<-"STGU3"
STGU3.obj@meta.data[,"Sample_name"]<-"Run5753_S3"

#caculate sd of area
STGU3.aera.sd<-sd(STGU3.obj@meta.data$Area)
STGU3.aera.mean<-mean(STGU3.obj@meta.data$Area)
STGU3.area.thresh<-c(STGU3.aera.mean+5*STGU3.aera.sd)

STGU3.double<-subset(STGU3.obj,subset = Area>=as.numeric(STGU3.area.thresh))
STGU3.double.1<-subset(STGU3.obj,subset = nCount_Nanostring>=2000)


pdf("~/Desktop/TLS/CosMx_counts/STGU3_double_QC.pdf",width = 6,height = 4)
VlnPlot(STGU3.obj, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample",pt.size = 0)
VlnPlot(STGU3.double, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
VlnPlot(STGU3.double.1, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
dev.off()

STGU3.f<-subset(STGU3.obj, subset = nCount_Nanostring<2000 & Area<as.numeric(STGU3.area.thresh))
dim(STGU3.f)
[1]   989 45020

save(STGU3.f,file = "STGU3_filtered.RData")

#Load S4
STGU4.obj <- LoadNanostring(
  data.dir = "/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S4",fov = "Run5753_S4"
)
dim(STGU4.obj)
[1]  1008 23458

STGU4.np <- grep("^NegPrb", rownames(STGU4.obj),value = T)

STGU4.obj<-STGU4.obj[setdiff(rownames(STGU4.obj),STGU4.np),]

STGU4.metadata<-read.csv("/Users/fanya562/Desktop/TLS/CosMx/Raw/Run5753_S4/Run5753_S4_metadata_file.csv",header = T,sep = ",")
STGU4.metadata$cell_ID_new<-paste0(STGU4.metadata$cell_ID,"_",STGU4.metadata$fov)
identical(colnames(STGU4.obj),STGU4.metadata$cell_ID_new)

STGU4.obj<-subset(STGU4.obj,cells = STGU4.metadata$cell_ID_new)
dim(STGU4.obj)
[1]   989 23458

STGU4.metadata.f<-subset(STGU4.metadata,STGU4.metadata$cell_ID_new%in%colnames(STGU4.obj))
identical(colnames(STGU4.obj),STGU4.metadata.f$cell_ID_new)
[1] TRUE

STGU4.obj@meta.data[4:24]<-STGU4.metadata.f[,1:21]
STGU4.obj@meta.data[,"Sample"]<-"STGU4"
STGU4.obj@meta.data[,"Sample_name"]<-"Run5753_S4"

#caculate sd of area
STGU4.aera.sd<-sd(STGU4.obj@meta.data$Area)
STGU4.aera.mean<-mean(STGU4.obj@meta.data$Area)
STGU4.area.thresh<-c(STGU4.aera.mean+5*STGU4.aera.sd)

STGU4.double<-subset(STGU4.obj,subset = Area>=as.numeric(STGU4.area.thresh))
STGU4.double.1<-subset(STGU4.obj,subset = nCount_Nanostring>=2000)


pdf("~/Desktop/TLS/CosMx_counts/STGU4_double_QC.pdf",width = 6,height = 4)
VlnPlot(STGU4.obj, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample",pt.size = 0)
VlnPlot(STGU4.double, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
VlnPlot(STGU4.double.1, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample")
dev.off()

STGU4.f<-subset(STGU4.obj, subset = nCount_Nanostring<2000 & Area<as.numeric(STGU4.area.thresh))
dim(STGU4.f)
[1]   989 23446

save(STGU4.f,file = "STGU4_filtered.RData")

#Merge 4 datasets
nano.obj<-merge(x=(STGG1.f),y=c(STGU1.f, STGU3.f, STGU4.f))
dim(nano.obj)
[1]    989 168099

cell.id.f<-read.table("/Users/fanya562/Desktop/TLS/CosMx_counts/cell_id_after_filtering.txt",sep = "\t")
nano.obj<-subset(nano.obj,cells=cell.id.f$V1)
dim(nano.obj)
[1]    989 165572

save(nano.obj,file = "4Sample_after_filtering.RData")


pdf("~/Desktop/TLS/CosMx_counts/4Sample_after_QC.pdf",width = 12,height = 6)
VlnPlot(nano.obj, features = c("nFeature_Nanostring", "nCount_Nanostring"), ncol = 2,group.by = "Sample",pt.size = 0)
dev.off()

write.table(colnames(nano.obj),file = "~/Desktop/TLS/CosMx_counts/Cell_id_after_QC.txt",row.names = F,col.names = F,quote = F)


options(future.globals.maxSize = 8000 * 1024^2)

nano.obj <- SCTransform(nano.obj, assay = "Nanostring", clip.range = c(-10, 10))
nano.obj <- RunPCA(nano.obj, npcs = 30, features = rownames(nano.obj))
ElbowPlot(nano.obj,ndims = 30)

nano.obj <- RunUMAP(nano.obj, dims = 1:30)
nano.obj <- FindNeighbors(nano.obj, reduction = "pca", dims = 1:30)
nano.obj <- FindClusters(nano.obj, resolution = 0.2)

pdf("UMAP_no_aglin.pdf",width = 8,height = 6)
DimPlot(nano.obj,label = T,raster=FALSE,pt.size = 0.1)
DimPlot(nano.obj,split.by = "Sample",ncol = 2,raster=FALSE,pt.size = 0.1)
dev.off()

####################################################################################################################################################################################
###################################################################### Harmony integration ######################################################################################### 
####################################################################################################################################################################################

#Do integration
nano.obj <- nano.obj %>% RunHarmony("Sample")

nano.obj <- FindNeighbors(nano.obj,reduction = "harmony", dims=1:30)
nano.obj <- RunUMAP(nano.obj,reduction = "harmony",  dims = 1:30)

nano.obj <- FindClusters(nano.obj,resolution = 0.6) 

pdf("UMAP_aglined_pc30_r0.6.pdf",width = 8,height = 6)
DimPlot(nano.obj,label=T,raster=FALSE,pt.size = 0.1)
DimPlot(nano.obj,split.by = "Sample",ncol=2,label = T,raster=FALSE,pt.size = 0.1)
dev.off()

pdf("Aglined_pc30_r0.6_QC.pdf",width = 10,height = 4)
VlnPlot(nano.obj,features = c("nCount_Nanostring","nFeature_Nanostring"),pt.size = 0)
dev.off()

DefaultAssay(nano.obj)<-"Nanostring"
all.cluster.markers<-FindAllMarkers(nano.obj,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
library(openxlsx)
write.xlsx(all.cluster.markers,file = "Aligned_pc30_r0.6_cluster_markers.xlsx",rowNames=F,colNames=T)

top5.cluster.markers<-all.cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(nano.obj,features = unique(top5.cluster.markers$gene),cols = "RdYlBu",assay = "SCT")+RotatedAxis()


cluster3.18.markers<-subset(top10.cluster.markers,top10.cluster.markers$cluster==3 | top10.cluster.markers$cluster==18)
DotPlot(nano.obj,features = unique(cluster3.18.markers$gene),cols = "RdYlBu",assay = "SCT")+RotatedAxis()

T.B.markers<-c("IL7R","CD3D","CD3E","CCL5","CD79A","MS4A1","IGKC","IGHM")
DotPlot(nano.obj,features = T.B.markers,cols = "RdYlBu",assay = "SCT")+RotatedAxis()

FeaturePlot(nano.obj,features = c("CD4","CD8A"),cols = c("gray","red"),raster = F,order = T)

cluster16.markerts<- FindConservedMarkers(nano.obj,assay = "Nanostring", ident.1 = "16", grouping.var = "Sample", verbose = FALSE)

save(nano.obj,file = "Aligned_pc30_r0.6_seurat.RData")

cluster.rev<-as.numeric(as.vector(nano.obj@meta.data$SCT_snn_res.0.6))
cluster.rev[cluster.rev %in% c(0,1,3,5,7,8,11,15,17,18,19)]<-"Tumor cell"
cluster.rev[cluster.rev %in% c(2,13,9)]<-"Macrophage"
cluster.rev[cluster.rev %in% c(4)]<-"Fibroblast"
cluster.rev[cluster.rev %in% c(6)]<-"T cell"
cluster.rev[cluster.rev %in% c(10)]<-"Mural cell"
cluster.rev[cluster.rev %in% c(12)]<-"Endothelial"
cluster.rev[cluster.rev %in% c(14)]<-"B cell"
cluster.rev[cluster.rev %in% c(20)]<-"Neutrophil"
cluster.rev[cluster.rev %in% c(16)]<-"Erythrocyte"

nano.obj@meta.data[,"cluster.rev"]<-cluster.rev

save(nano.obj,file = "Aligned_pc30_r0.6_revised.RData")

pdf("UMAP_aglined_pc30_r0.6_revised.pdf",width = 8,height = 6)
DimPlot(nano.obj,group.by = "cluster.rev",raster = F,pt.size = 0.1)+theme(plot.title = element_blank())
dev.off()

metadata.rev<-nano.obj@meta.data
metadata.rev<-metadata.rev[,-c(1,29,31)]
write.table(metadata.rev,file = "metadata_revised.csv",sep = ",",row.names = T,col.names = T )

#Makers of each celltype
Idents(nano.obj)<-"cluster.rev"

all.celltype.markers<-FindAllMarkers(nano.obj,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
write.xlsx(all.celltype.markers,file="all.celltype.markers.xlsx",rowNames=F,colNames=T)

top5.celltype.markers<- all.celltype.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(nano.obj,features = unique(top5.celltype.markers$gene),cols="RdYlBu",assay = "SCT")+RotatedAxis()


##########################################################################################################################################################################
###################################################################### Macrophage and Fibroblast subcluster ##############################################################
##########################################################################################################################################################################
nano.subset<-subset(nano.obj,idents = c("Macrophage","Fibroblast"))

nano.subset <- SCTransform(nano.subset, assay = "Nanostring", clip.range = c(-10, 10))
nano.subset <- RunPCA(nano.subset, npcs = 30, features = rownames(nano.obj))
ElbowPlot(nano.subset,ndims = 30)

#run harmony
nano.subset<- nano.subset %>% RunHarmony("Sample")

nano.subset <- FindNeighbors(nano.subset,reduction = "harmony", dims=1:15)
nano.subset <- RunUMAP(nano.subset,reduction = "harmony",  dims = 1:15)

nano.subset <- FindClusters(nano.subset,resolution = 0.4) 

DimPlot(nano.subset,label = T)
DimPlot(nano.subset,label = T,group.by = "cluster.rev")

#Updata new annotation to all cells seurat object
sub.annotation<-as.numeric(as.vector(nano.subset@meta.data$SCT_snn_res.0.4))
sub.annotation[sub.annotation%in%c(3,8,9)]<-"Fibroblast"
sub.annotation[sub.annotation%in%c(0:2,4:7)]<-"Myeloids"

nano.subset@meta.data[,"sub.annotation"]<-sub.annotation

save(nano.subset,file = "Macrophage_Fibroblast_subcluster.RData")

cell.sub.anno<-data.frame("anno"=nano.subset@meta.data$sub.annotation,
                          "cluster"=nano.subset@meta.data$SCT_snn_res.0.4,
                          row.names = colnames(nano.subset))

cell.anno<-data.frame("anno"=nano.obj@meta.data$cluster.rev,
                      "cluster"=nano.obj@meta.data$SCT_snn_res.0.6,
                      row.names = colnames(nano.obj))

Mac.Fibro.anno<-cell.anno[rownames(cell.sub.anno),]
Mac.Fibro.anno$anno<-cell.sub.anno$anno

other.anno<-cell.anno[setdiff(rownames(cell.anno),rownames(Mac.Fibro.anno)),]

all.anno.new<-rbind(other.anno,Mac.Fibro.anno)
all.anno.new<-all.anno.new[colnames(nano.obj),]

identical(rownames(all.anno.new),colnames(nano.obj))
[1] TRUE

nano.obj@meta.data[,"celltype.new"]<-all.anno.new$anno

Idents(nano.obj)<-"celltype.new"



####################################################################################################################################################################################
###################################################################### Gene signature score ######################################################################################## 
####################################################################################################################################################################################
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



marker_score<-AddModuleScore_UCell(nano.obj,features = markers)


library(stringr)
library(ggplot2)
library(viridis)

a<-colnames(marker_score@meta.data) %>% str_subset("_UCell")
FeaturePlot(marker_score,features = a,order = T,ncol = 2,raster = F,cols = viridis(256))

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

FeaturePlot(marker_score,features = "test_UCell",order = T,ncol = 1,raster = F,cols = viridis(256))+
  theme(plot.title = element_blank())+labs(color=c("Neuron\nSignature\nScore"))



####################################################################################################################################################################################
###################################################################### T cell label transfer ####################################################################################### 
####################################################################################################################################################################################


############################################## Data preparation ################################################

#get negative probe expression
#cd.subset<-subset(cd,cd$is_np=="TRUE"&cd$name=="avg"&cd$value<=0.5)
#rownames(cd.subset)<-paste0(cd.subset$cell_ID,"_",cd.subset$fov,"_",substring(cd.subset$sample_id,10))

#cd.subset<-cd.subset[cell.id.f,]
#rownames(cd.subset)<-paste0(cd.subset$cell_ID,"_",cd.subset$fov,"_",substring(cd.subset$sample_id,10))
#write.xlsx(cd.subset,file = "./CosMx_counts/QC/NP_avg.xlsx",rowNames=T,colNames=T)

#load scRNA reference data
#scRNA.reference<-read.table("./scRNA/exported.txt",header = T,sep = "\t")
#head(scRNA.reference[1:5,1:20])
#rownames(scRNA.reference)<-paste0(scRNA.reference[,2],"_",scRNA.reference[,1])
#scRNA.meta<-scRNA.reference[,c(1:15)]
#scRNA.counts<-scRNA.reference[,c(16:ncol(scRNA.reference))]
#scRNA.counts<-t(scRNA.counts)
#head(scRNA.counts[1:5,1:5])

#scRNA.seurat<-CreateSeuratObject(scRNA.counts,min.cells = 0,min.features = 0)
#dim(scRNA.seurat)

#scRNA.seurat@meta.data[,4:19]<-scRNA.meta

#Idents(scRNA.seurat)<-"T.cell.subtypes"
#T.ref.seurta<-subset(scRNA.seurat,idents = c("CD4 T cells","CD8 T cells"))
#save(T.ref.seurta,file = "./CosMx_counts/InsituType/Tcell_scRNA.RData"


############################################## Load Data #######################################################

#load required data
cd.subset<-read.xlsx("./CosMx_counts/QC/NP_avg.xlsx",rowNames=T,colNames=T)
load(file="./CosMx_counts/Aligned_pc30_r0.6_new_annotation.RData")
load(file = "./CosMx_counts/InsituType/Tcell_scRNA.RData")

T.CosMx<-subset(nano.obj,idents=("T cell"))
dim(T.CosMx)
[1]  989 8661

T.neg.avg<-cd.subset[colnames(T.CosMx),]
rownames(T.neg.avg)<-paste0(T.neg.avg$cell_ID,"_",T.neg.avg$fov,"_",substring(T.neg.avg$sample_id,10))

identical(colnames(T.CosMx),rownames(T.neg.avg))
[1] TRUE

T.neg<-T.neg.avg$value
names(T.neg)<-rownames(T.neg.avg)
head(T.neg)


T.subcluster.ref.exp<-as.matrix(data.frame(AverageExpression(T.ref.seurta)))
dim(T.subcluster.ref.exp)
[1] 27112     2

colnames(T.subcluster.ref.exp)<-c("CD8.T.cells","CD4.T.cells")


###################################### Label Transfer with full gene list ######################################

counts <- t(as.data.frame(T.CosMx@assays$Nanostring@counts))
str(counts)

immunofluordata <- matrix(rpois(n = nrow(counts) * 4, lambda = 100), 
                          nrow(counts))

#Use metadata as reference


# perform automatic cohorting:
cohort <- fastCohorting(immunofluordata,
                        gaussian_transform = TRUE) 
# ("Gaussian_transform = TRUE" maps variables to gaussians in order to 
#  place dramatically different variables on the same scale.)
table(cohort)


sup <- insitutypeML(x =counts,
                    neg=T.neg,
                    cohort = cohort,
                    reference_profiles = T.subcluster.ref.exp)   

str(sup)
round(head(sup$prob), 2)
heatmap(sweep(sup$profiles, 1, pmax(apply(sup$profiles, 1, max), .2), "/"), scale = "none",
        main = "Mean cell type expression profiles")


cols <- colorCellTypes(freqs = table(sup$clust), palette = "brewers")

par(mfrow = c(1, 2))
par(mar = c(0, 0, 3, 0))

cluster<-sup[1]
T.CosMx@meta.data[,"T_sub"]<-cluster[["clust"]]

DimPlot(T.CosMx,group.by = "T_sub")
Idents(T.CosMx)<-"T_sub"

table(T.CosMx@meta.data$T_sub)
CD4.T.cells CD8.T.cells 
4218        4443 

DefaultAssay(T.CosMx)
cd4.vs.cd8<-FindAllMarkers(T.CosMx,min.pct = 0.25,logfc.threshold = 0.25,only.pos = T)


#Update T cell Label transfer result to main seurat object
nano.obj.with.label.transfer<-nano.obj

Idents(T.CosMx)<-"T_sub"
Idents(nano.obj.with.label.transfer, cells = colnames(T.CosMx)) <- Idents(T.CosMx)
nano.obj.with.label.transfer@meta.data$LT_full_gene<-Idents(nano.obj.with.label.transfer)
Idents(nano.obj.with.label.transfer)<-"LT_full_gene"

####################################################################################################################################################################################
###################################################################### CD4 T cell trajectory ####################################################################################### 
####################################################################################################################################################################################

Idents(nano.obj.with.label.transfer)<-"celltype.new"
Tcell.CosMx<-subset(nano.obj.with.label.transfer,idents="T cell")

Tcell.CosMx<-subset_opt(Tcell.CosMx,cells = c(Tcell.CosMx@images$Run5753_S1$centroids@cells,
                                              Tcell.CosMx@images$Run5753_S2$centroids@cells,
                                              Tcell.CosMx@images$Run5753_S3$centroids@cells))


T_TLS.fov.id<-c(Tcell.CosMx@images$STGG1_fov1@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov2@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov5@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov6@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov7@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov1@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov3@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov1@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov2@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov4@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov6@boundaries$centroids@cells)

B_TLS.fov.id<-c(Tcell.CosMx@images$STGU1_fov8@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov9@boundaries$centroids@cells)

M_TLS.fov.id<-c(Tcell.CosMx@images$STGU1_fov26@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov5@boundaries$centroids@cells)

T_PN.fov.id<-c(Tcell.CosMx@images$STGG1_fov3@boundaries$centroids@cells,
               Tcell.CosMx@images$STGG1_fov22@boundaries$centroids@cells,
               Tcell.CosMx@images$STGG1_fov23@boundaries$centroids@cells,
               Tcell.CosMx@images$STGU1_fov6@boundaries$centroids@cells,
               Tcell.CosMx@images$STGU1_fov16@boundaries$centroids@cells,
               Tcell.CosMx@images$STGU1_fov28@boundaries$centroids@cells)

Tumor.fov.id<-c(Tcell.CosMx@images$STGG1_fov9@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov10@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov11@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov12@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov13@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov14@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov15@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov16@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov17@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov18@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov19@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov20@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov21@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov24@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov25@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov26@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov27@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov28@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov29@boundaries$centroids@cells,
                Tcell.CosMx@images$STGG1_fov31@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov2@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov4@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov7@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov10@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov11@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov12@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov13@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov14@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov15@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov17@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov18@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov19@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov20@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov21@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov22@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov23@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov24@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov25@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov27@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov29@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov30@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov31@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU1_fov32@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov3@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov7@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov8@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov9@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov10@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov11@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov12@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov14@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov15@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov16@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov17@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov18@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov19@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov20@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov21@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov22@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov23@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov24@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov25@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov26@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov27@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov28@boundaries$centroids@cells,
                Tcell.CosMx@images$STGU3_fov29@boundaries$centroids@cells)


Tumor.TLS.PN.Tcell.CosMx.<-subset_opt(Tcell.CosMx,cells = c(T_TLS.fov.id,B_TLS.fov.id,M_TLS.fov.id,T_PN.fov.id,Tumor.fov.id))
Tumor.TLS.PN.Tcell.CosMx<-Tumor.TLS.PN.Tcell.CosMx.

save(Tumor.TLS.PN.Tcell.CosMx,file="Tumor.TLS.PN.Tcell.CosMx.RData")


Idents(Tumor.TLS.PN.Tcell.CosMx)<-"celltype.new"
Idents(T.CosMx)<-"T_sub"
Idents(Tumor.TLS.PN.Tcell.CosMx, cells = colnames(T.CosMx)) <- Idents(T.CosMx)

CD4.vs.CD8.DEGs<-FindMarkers(Tumor.TLS.PN.Tcell.CosMx,ident.1 = "CD4.T.cells",ident.2 = "CD8.T.cells",min.pct = 0.25,logfc.threshold = 0.25)

Tcell.subtype<-as.vector(Tumor.TLS.PN.Tcell.CosMx@meta.data$LT_full_gene)
names(Tcell.subtype)<-rownames(Tumor.TLS.PN.Tcell.CosMx@meta.data)

CD4.Tcell.id<-names(Tcell.subtype[Tcell.subtype%in%"CD4.T.cells"])
CD8.Tcell.id<-names(Tcell.subtype[Tcell.subtype%in%"CD8.T.cells"])

Tumor.TLS.PN.CD4.CosMx<-subset_opt(Tumor.TLS.PN.Tcell.CosMx, cells = CD4.Tcell.id)
Tumor.TLS.PN.CD8.CosMx<-subset_opt(Tumor.TLS.PN.Tcell.CosMx, cells = CD8.Tcell.id)

region<-rownames(Tumor.TLS.PN.CD4.CosMx@meta.data)
region[region %in% T_TLS.fov.id ]<-"T_TLS"
region[region %in% M_TLS.fov.id ]<-"M_TLS"
region[region %in% B_TLS.fov.id ]<-"B_TLS"
region[region %in% T_PN.fov.id ]<-"T_PN"
region[region %in% Tumor.fov.id ]<-"Tumor"

Tumor.TLS.PN.CD4.CosMx@meta.data$region<-region

save(Tumor.TLS.PN.CD4.CosMx,file = "./CD4 in 3 regions/Tumor.TLS.PN.CD4.CosMx.RData")
save(Tumor.TLS.PN.CD8.CosMx,file = "./CD8 in 3 regions/Tumor.TLS.PN.CD8.CosMx.RData")
#load CD4 related gene list
#CD4.gene.list<-read.table("CD4_related_genes.txt",sep = "\t",header = F)
#CD4.gene.list<-unique(c(toupper(CD4.gene.list$V1)))
#CD4.gene.list.f<-intersect(rownames(gene_ann),CD4.gene.list)
#write.table(CD4.gene.list.f,file = "CD4_genes_in_CosMx.txt",sep = "\t",quote = F,row.names = F,col.names = F)

CD4.gene.list.f<-read.table("Final list for CD4 trajectory analysis.txt",sep = "\t",header = F)
CD4.gene.list.f<-intersect(rownames(gene_ann),CD4.gene.list.f$V1)

write.table(CD4.gene.list.f,file = 'Gene_list_for_CD4_trajectory.txt',sep = "\t",quote = F)


library(monocle)
CD4_ann<-Tumor.TLS.PN.CD4.CosMx@meta.data
CD4_ann$celltype<-Idents(Tumor.TLS.PN.CD4.CosMx)


gene_ann<-data.frame(gene_short_name=rownames(Tumor.TLS.PN.CD4.CosMx@assays$Nanostring),
                     row.names = rownames(Tumor.TLS.PN.CD4.CosMx@assays$Nanostring))

save(CD4_ann,file = "./CD4 in 3 regions/CD4_ann.RData")
save(gene_ann,file = "./CD4 in 3 regions/gene_ann.RData")


pd<-new("AnnotatedDataFrame",data=CD4_ann)
fd<-new("AnnotatedDataFrame",data=gene_ann)
ce=as.data.frame(Tumor.TLS.PN.CD4.CosMx@assays$Nanostring@counts)

save(ce,file = "CD4_ce.RData")

cds<-newCellDataSet(as.matrix(ce),
                    phenoData = pd,
                    featureData = fd)


cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)


cds<-reduceDimension(cds,max_components = 2,num_dim=6,reduction_method = "tSNE",verbose = T)
cds<-clusterCells(cds,1,2,num_clusters=5)

pData(cds)$Cluster=pData(cds)$celltype

cds<-setOrderingFilter(cds,CD4.gene.list.f)
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,reduction_method="DDRTree",residualModelFormulaStr=c("~Sample"))
cds<-orderCells(cds)

save(cds,file = "CD4_pseudotime_object.RData")
