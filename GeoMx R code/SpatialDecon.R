spatial.norm.counts<-read.xlsx("Q3 Norm_All ROIs.xlsx",sheet = 3,rowNames = F,colNames = T)
spatial.raw.counts<-read.xlsx("Filtered_counts.xlsx",sheet = 3,rowNames = F,colNames = T)


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

spatial.norm.counts.f<-spatial.norm.counts[,c("TargetName",metadata.f$SegmentDisplayName)]
spatial.raw.counts.f<-spatial.raw.counts[,c("TargetName",metadata.f$SegmentDisplayName)]

featureAnnoFile<-read.xlsx("Q3 Norm_All ROIs.xlsx",sheet = 4,colNames = T,rowNames = F)

spe.norm <- readGeoMx(spatial.norm.counts.f, metadata.f,rmNegProbe = F)

####deconvolution
spatialdecon.norm<-prepareSpatialDecon(spe.norm,assay2use = "counts")

library(SpatialDecon)
safe_TME<-read.table("safeTME-for-tumor-immune.csv",sep = ",",header = T)
rownames(safe_TME)<-safe_TME$X
safe_TME<-safe_TME[,-1]

raw.counts<-spatial.raw.counts.f
rownames(raw.counts)<-raw.counts$TargetName
raw.counts<-raw.counts[,-1]

res.norm = spatialdecon(norm = as.matrix(spatialdecon.norm$normCount),
                        raw=as.matrix(raw.counts),
                        bg = as.matrix(spatialdecon.norm$backGround),
                        X = as.matrix(safe_TME),
                        cell_counts=spe.norm@colData$AOINucleiCount,
                        align_genes = TRUE)
str(res.norm)

restils.norm = spatialdecon(norm = as.matrix(spatialdecon.norm$normCount),                     # normalized data                    
                            bg = as.matrix(spatialdecon.norm$backGround),   # expected background counts for every data point in norm
                            raw= as.matrix(raw.counts), # raw data, used to down-weight low-count observations
                            X = as.matrix(safeTME), # safeTME matrix, used by default
                            cell_counts=spe.norm@colData$AOINucleiCount, 
                            cellmerges = safeTME.matches)
str(restils.norm)

layout(mat = (matrix(c(1, 2), 1)), widths = c(7, 3))
TIL_barplot(restils.norm$prop_of_nontumor, 
            draw_legend = T, cex.names = 0.45)

spatialdecon.result<-data.frame(restils.norm[["prop_of_nontumor"]])
library(openxlsx)
write.xlsx(spatialdecon.result,"SpatialDecon_with_Q3_norm.xlsx",rowNames=T,colNames=T)



