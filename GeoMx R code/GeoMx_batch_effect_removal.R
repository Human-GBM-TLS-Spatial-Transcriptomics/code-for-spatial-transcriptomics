library(standR)
library(SpatialExperiment)
library(limma)
library(ExperimentHub)
library(openxlsx)

spatial.norm.counts<-read.xlsx("Q3 Norm_All ROIs.xlsx",sheet = 3,rowNames = F,colNames = T)
TargetName<-spatial.norm.counts$TargetName
spatial.norm.counts<-spatial.norm.counts[,-1]
log2.spatial.norm.counts<-log2(spatial.norm.counts+1)
log2.spatial.norm.counts<-cbind(TargetName,log2.spatial.norm.counts)

metadata<-read.xlsx("Q3 Norm_All ROIs.xlsx",sheet = 1,rowNames = F,colNames = T)
metadata$SegmentDisplayName<-gsub(" ",".",metadata$SegmentDisplayName)
metadata$SampleName<-paste0(metadata$SlideName,"_",metadata$ROILabel)
metadata$SampleName<-gsub(" ","",metadata$SampleName)
sample.name<-paste0(metadata$SlideName,"_",metadata$ROILabel)
sample.name<-gsub(" ","",sample.name)
  

annotation<-read.xlsx("221011 Annotation template file.xlsx",rowNames = F,colNames = T)

annotation$sample.name<-paste0(annotation$Patient.ID,"_",annotation$`ROI.(label)`)
annotation$sample.name<-gsub(" ","",annotation$sample.name)

####extract TLS and PN samples
annotation.f<-subset(annotation,annotation$ClusterType==c("TLS")|annotation$ClusterType==c("PN"))
annotation.f<-subset(annotation.f,annotation.f$Composition=="T"|annotation.f$Composition=="B"|annotation.f$Composition=="M")
annotation.f$group<-paste0(annotation.f$ClusterType,"_",annotation.f$Composition)

metadata.f<-subset(metadata,metadata$SampleName%in%annotation.f$sample.name)

spatial.norm.counts.f<-log2.spatial.norm.counts[,c("TargetName",metadata.f$SegmentDisplayName)]

featureAnnoFile<-read.xlsx("Q3 Norm_All ROIs.xlsx",sheet = 4,colNames = T,rowNames = F)

spe.norm <- readGeoMx(spatial.raw.counts.f, metadata.f,rmNegProbe = F)


SummarizedExperiment::assayNames(spe.norm)

library(ggalluvial)

colData(spe.norm)$regions <- paste0(colData(spe.norm)$ClusterType,"_",colData(spe.norm)$Composition)


spe.norm <- addPerROIQC(spe.norm, rm_genes = F)
metadata(spe.norm) |> names()

pdf("CellCounts_with_Reads_cor.pdf",width = 6,height = 6)
plotROIQC(spe.norm, y_threshold = 50000, col = SlideName)+labs(title = "Q3 Norm Counts")
dev.off()

plotRLExpr(spe.norm)

pdf("QC.pdf",width = 6,height = 6)

plotRLExpr(spe.norm, ordannots = "SlideName", assay = 1, col = SlideName)+labs(title = "Q3 Norm Counts")

dev.off()

pdf("PCA_before_BatchRemoval.pdf",width = 6,height = 4)
drawPCA(spe.norm, assay = 1, col = SlideName, shape = regions)+labs(title = "Q3 Norm Counts")
dev.off()


plotPairPCA(spe.norm, col = regions,n_dimension = 3, 
            shape = SlideName, assay = 1)


spe.norm <- findNCGs(spe.norm,n_assay = 1, batch_name = "SlideName", top_n = 500)
metadata(spe.norm) |> names()

findBestK(spe.norm, maxK = 10, factor_of_int = "regions", NCGs = metadata(spe.norm)$NCGs, factor_batch = "SlideName")

spe.norm.ruv <- geomxBatchCorrection(spe.norm, n_assay = 1,factors = "regions", 
                                    NCGs = metadata(spe.norm)$NCGs, k = 6)

pdf("Log2_Q3_PCA_after_RUV4.pdf",width = 10,height = 8)

plotPairPCA(spe.norm.ruv, assay = 2, color =regions , shape = SlideName, title = "RUV4 with Log2 Q3 Norm Counts")

dev.off()



Q3.norm.counts.BR<-data.frame(spe.norm.ruv@assays@data@listData[["logcounts"]])
write.xlsx(Q3.norm.counts.BR,file = "Q3_Norm_counts_after_BR.xlsx",rowNames=T,colNames=T)


spe.norm.list <- list(spe.norm, spe.norm.ruv)

pdf("BatchRemoval_evaluation.pdf",width = 12,height = 8)

plotClusterEvalStats(spe_list = spe.norm.list,
                     bio_feature_name = "regions",
                     batch_feature_name = "SlideName",
                     data_names = c("Raw","RUV4"))

plotRLExpr(spe.norm.ruv, assay = 2, color = SlideName) + ggtitle("RUV4 with Q3 Norm")


dev.off()

