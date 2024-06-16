library(BiocParallel)
library(data.table)
library(magick)
library(SingleCellExperiment)

names(sid) <- sid <- list.dirs("~/Desktop/TLS/CosMx_Counts", FALSE, FALSE) # sample IDs
names(iid) <- iid <- c("morphology", "segmentation") # image IDs

#load Expression data
dir <- list.dirs("~/Desktop/TLS/CosMx_Counts", recursive = FALSE)
sce <- lapply(sid, \(s) {
  y <- fread(file.path("~/Desktop/TLS/CosMx_Counts", s, paste0(s,"_","exprMat_file.csv")))
  cd <- fread(file.path("~/Desktop/TLS/CosMx_Counts", s, paste0(s,"_","metadata_file.csv")))
  # subset & sparsify
  y <- y[y$cell_ID %in% cd$cell_ID, ]
  y <- as.matrix(y[, -seq(2)])
  y <- as(t(y), "dgCMatrix")
  # construct SCE
  gs <- rownames(y)
  rd <- DataFrame(gene = gs, sample_id = s)
  cs <- paste0("cell", seq(ncol(y)))
  cd$cell <- colnames(y) <- cs
  g <- gsub("[0-9]+", "", s)
  cd$sample_id <- s
  cd$group_id <- g
  sce <- SingleCellExperiment(
    assays = list(counts = y),
    rowData = rd, colData = cd)
  # separate negative probes
  np <- grepl("^NegPrb", gs)
  altExp(sce, "np") <- sce[np, ]
  sce <- sce[!np, ]
})


library(dplyr)
library(tidyr)
library(tidytext)
library(SingleCellExperiment)

sce <- bplapply(sce, BPPARAM = bpparam(), \(sce) {
  lys <- list(sce, altExp(sce))
  lys <- lapply(lys, \(sce) {
    z <- (y <- assay(sce)) != 0 
    # gene-level
    rd <- data.frame(
      sum = rowSums(y),
      avg = rowMeans(y),
      det = rowSums(z),
      frq = rowMeans(z))
    # cell-level
    cd <- data.frame(
      sum = colSums(y),
      avg = colMeans(y),
      det = colSums(z),
      frq = colMeans(z))
    rowData(sce) <- cbind(rowData(sce), rd)
    colData(sce) <- cbind(colData(sce), cd)
    return(sce)
  })
  sce <- lys[[1]]
  altExp(sce) <- lys[[2]]
  return(sce)
})

lt <- c("sum", "det")
qc <- c(lt, c("avg", "frq"))
# extract gene/cell metadata,
# log-transform & reformat
.md <- \(dat, dim) {
  .md <- list(rowData, colData)[[dim]]
  lys <- lapply(dat, \(.) data.frame(.md(.)))
  do.call(rbind, lys) %>% 
    mutate_at(all_of(lt), log1p) %>% 
    pivot_longer(all_of(qc))
}
# join & tidy data from 
# reporters & negative controls
.df <- \(dat, dim) {
  list(
    "FALSE" = .md(dat[[1]], dim),
    "TRUE" = .md(dat[[2]], dim)) %>% 
    bind_rows(.id = "is_np") %>% 
    mutate_at("is_np", as.logical)
}
# overview of cell & reporter 
# / negative control counts
.tbl <- \(df, val = "sum")
filter(df, name == val) %>% 
  group_by(is_np, sample_id) %>% 
  summarise(
    .groups = "drop", n = n(),
    mean = round(mean(value), 2)) %>% 
  pivot_wider(
    names_from = "is_np", 
    values_from = "mean") %>% 
  dplyr::rename(
    control = "TRUE",
    reporter = "FALSE")


# Quality control
# separate negative controls
dat <- list(sce, lapply(sce, altExp))
# prettify gene/cell metadata
rd <- .df(dat, 1)
cd <- .df(dat, 2)

library(ggplot2)
## Gene-level
pdf("~/Desktop/TLS/CosMx_Counts/CosMx_QC.pdf",width = 10,height = 6)
ggplot(rd, aes(y = value, fill = is_np,
               reorder_within(sample_id, value, name))) +
  facet_wrap(~ name, scales = "free") +
  geom_boxplot(size = 0.2, 
               outlier.size = 0.2, outlier.shape = 16) +
  scale_x_reordered() + theme_bw() + theme(
    axis.title = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1))

## Cell-level
gg <- cd %>% filter(
  is_np, name == "sum") %>% 
  group_by(sample_id, fov) %>% 
  summarise_at("value", mean)
ggplot(gg, aes(factor(fov), sample_id, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(
    "mean negative\ncontrol count",
    option = "A", limits = c(0, NA)) +
  coord_equal(expand = FALSE) +
  theme_bw() + theme(
    axis.title = element_blank(),
    panel.grid = element_blank())

gg <- cd %>% filter(
  !is_np, name == "frq") %>% 
  group_by(sample_id, fov) %>% 
  summarise_at("value", mean)
ggplot(gg, aes(factor(fov), sample_id, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(
    "mean RNA\ndetection\nfrequency",
    option = "A", limits = c(0, NA)) +
  coord_equal(expand = FALSE) +
  theme_bw() + theme(
    axis.title = element_blank(),
    panel.grid = element_blank())

ggplot(cd, aes(y = value, fill = is_np,
               reorder_within(sample_id, value, name))) +
  facet_wrap(~ name, scales = "free") +
  geom_boxplot(size = 0.2, 
               outlier.size = 0.2, outlier.shape = 16) +
  scale_x_reordered() + theme_bw() + theme(
    axis.title = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.key.size = unit(0.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

sub <- BiocParallel::bplapply(sce, BPPARAM = bpparam(), \(sce) {
  # keep cells with average negative control count below 0.5
  sce <- sce[, colData(altExp(sce))$avg <= 0.5]
  # keep cells with at least 20 detected features
  sce <- sce[, colSums(assay(sce) > 0) >= 20]
  return(sce)
})

list(
  before = .tbl(.df(list(sce, lapply(sce, altExp)), 2)),
  after = .tbl(.df(list(sub, lapply(sub, altExp)), 2))) %>% 
  bind_rows(.id = ".") %>% 
  pivot_wider(
    names_from = ".", names_sep = ".", 
    values_from = c("n", "reporter", "control")) %>% 
  mutate(
    .after = "n.after", 
    "%" = round(100 * n.after / n.before, 2)) %>% 
  knitr::kable()

|sample_id  | n.before| n.after|     %| reporter.before| reporter.after| control.before| control.after|
  |:----------|--------:|-------:|-----:|---------------:|--------------:|--------------:|-------------:|
  |P240 |    39590|   39241| 99.12|            5.77|           5.79|           0.65|          0.65|
  |P237 |    60259|   59612| 98.93|            5.93|           5.95|           0.83|          0.83|
  |P238 |    45062|   43937| 97.50|            5.28|           5.34|           0.61|          0.62|
  |P239 |    23462|   22985| 97.97|            5.53|           5.58|           0.73|          0.74|


sub <- lapply(sub, \(sce) {
  altExp(sce) <- NULL
  rowData(sce) <- NULL
  return(sce)
})

(sce <- do.call(cbind, sub))

cell.id.f<-paste0(sce@colData$cell_ID,"_",sce@colData$fov,"_",substring(sce@colData$sample_id,10))
library(openxlsx)
write.table(cell.id.f,file = "cell_id_after_filtering.txt",quote=F,row.names = F,col.names = F)

rm(list = setdiff(ls(),"cell.id.f"))


#######################################################################################################
###########################################  Session Info  ############################################
#######################################################################################################
> sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] ggplot2_3.4.2               tidytext_0.4.1              tidyr_1.3.0                
[4] dplyr_1.1.2                 SeuratObject_4.1.3          Seurat_4.3.0               
[7] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0             
[10] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0             
[13] S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
[16] matrixStats_0.63.0          magick_2.7.4                data.table_1.14.8          
[19] BiocParallel_1.32.6        

loaded via a namespace (and not attached):
  [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6          
[4] ellipsis_0.3.2         ggridges_0.5.4         XVector_0.38.0        
[7] spatstat.data_3.0-1    rstudioapi_0.14        farver_2.1.1          
[10] leiden_0.4.3           listenv_0.9.0          SnowballC_0.7.1       
[13] ggrepel_0.9.3          fansi_1.0.4            codetools_0.2-18      
[16] splines_4.2.0          knitr_1.42             polyclip_1.10-4       
[19] jsonlite_1.8.4         ica_1.0-3              cluster_2.1.3         
[22] png_0.1-8              uwot_0.1.14            spatstat.sparse_3.0-1 
[25] shiny_1.7.4            sctransform_0.3.5      BiocManager_1.30.20   
[28] compiler_4.2.0         httr_1.4.5             Matrix_1.5-4          
[31] fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.1             
[34] later_1.3.1            htmltools_0.5.5        tools_4.2.0           
[37] igraph_1.4.2           gtable_0.3.3           glue_1.6.2            
[40] GenomeInfoDbData_1.2.9 RANN_2.6.1             reshape2_1.4.4        
[43] Rcpp_1.0.10            scattermore_0.8        vctrs_0.6.2           
[46] nlme_3.1-157           spatstat.explore_3.2-1 progressr_0.13.0      
[49] lmtest_0.9-40          spatstat.random_3.1-5  xfun_0.39             
[52] stringr_1.5.0          globals_0.16.2         mime_0.12             
[55] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1         
[58] goftest_1.2-3          future_1.32.0          zlibbioc_1.44.0       
[61] MASS_7.3-56            zoo_1.8-12             scales_1.2.1          
[64] promises_1.2.0.1       spatstat.utils_3.0-3   parallel_4.2.0        
[67] RColorBrewer_1.1-3     reticulate_1.28        pbapply_1.7-0         
[70] gridExtra_2.3          stringi_1.7.12         tokenizers_0.3.0      
[73] rlang_1.1.0            pkgconfig_2.0.3        bitops_1.0-7          
[76] lattice_0.20-45        tensor_1.5             ROCR_1.0-11           
[79] purrr_1.0.1            labeling_0.4.2         patchwork_1.1.2       
[82] htmlwidgets_1.6.2      cowplot_1.1.1          tidyselect_1.2.0      
[85] parallelly_1.35.0      RcppAnnoy_0.0.20       plyr_1.8.8            
[88] magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
[91] DBI_1.1.3              DelayedArray_0.24.0    withr_2.5.0           
[94] pillar_1.9.0           fitdistrplus_1.1-11    abind_1.4-5           
[97] survival_3.3-1         RCurl_1.98-1.12        sp_1.6-0              
[100] tibble_3.2.1           future.apply_1.10.0    crayon_1.5.2          
[103] janeaustenr_1.0.0      KernSmooth_2.23-20     utf8_1.2.3            
[106] spatstat.geom_3.2-1    plotly_4.10.1          grid_4.2.0            
[109] digest_0.6.31          xtable_1.8-4           httpuv_1.6.10         
[112] munsell_0.5.0          viridisLite_0.4.1    


