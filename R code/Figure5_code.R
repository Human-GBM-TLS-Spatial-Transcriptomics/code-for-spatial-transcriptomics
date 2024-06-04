####################################################################################################################################################################################
###################################################################### code for figure5 ####################################################################################### 
####################################################################################################################################################################################

#############################################################Figure 5B
p13<-ImageDimPlot(nano.obj.with.label.transfer, fov = "STGU1_fov16",size = 0.5, coord.fixed = FALSE,axes = T,
                  cells = WhichCells(nano.obj.with.label.transfer,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                  cols=c("red","#AC844B","cyan","magenta","white"),border.color = "black")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))

STGU1.fov16.cell.frac<-data.frame(table(STGU1.fov16.obj@meta.data$celltype.new))
STGU1.fov16.cell.frac$FOV<-"M-TLS"
STGU1.fov16.cell.frac<-STGU1.fov16.cell.frac[-c(4,7),]

#
STGU1.fov26.id<-as.vector(nano.obj.with.label.transfer@images$STGU1_fov26@boundaries$centroids@cells)
STGU1.fov26.obj<-subset_opt(nano.obj.with.label.transfer, cells =STGU1.fov26.id)

table(STGU1.fov26.obj@meta.data$celltype.new)
B cell  Endothelial  Erythrocyte   Fibroblast   Mural cell Myeloid cell   Neutrophil 
247           57           63           47           41          445           17 
T cell   Tumor cell 
244         1160 

colnames(STGU1.fov16.cell.frac) = c('Celltype','Number of cell','Group')
STGU1.fov16.cell.frac<-STGU1.fov16.cell.frac[order(STGU1.fov16.cell.frac$`Number of cell`,decreasing = F),]
STGU1.fov16.cell.frac$Celltype<-factor(STGU1.fov16.cell.frac$Celltype,levels = STGU1.fov16.cell.frac$Celltype)


p14<-ggplot(STGU1.fov16.cell.frac)+geom_col(aes(x=Celltype,y=STGU1.fov16.cell.frac$`Number of cell`,fill=Celltype),width = 0.8)+
  mytheme+coord_flip()+ylab("Number of cells")+NoLegend()+theme(plot.margin = margin(0,0,0,-1))+
  scale_fill_manual(values = c("magenta","red","white","cyan","#AC844B"))

svg("Figure5B_STGU1_fov16_with_cell_fraction.svg",width = 8,height = 6)
p13+p14+plot_layout(widths = c(5, 1))
dev.off()

#############################################################Figure 5C
p15<-ImageDimPlot(nano.obj.with.label.transfer, fov = "STGG1_fov22",size = 0.5, coord.fixed = FALSE,axes = T,
                  cells = WhichCells(nano.obj.with.label.transfer,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                  cols=c("red","#AC844B","cyan","magenta","white"),border.color = "black")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))

STGG1.fov22.id<-as.vector(nano.obj@images$STGG1_fov22@boundaries$centroids@cells)
STGG1.fov22.obj<-subset_opt(nano.obj, cells =STGG1.fov22.id)

table(STGG1.fov22.obj@meta.data$celltype.new)
B cell  Endothelial  Erythrocyte   Fibroblast   Mural cell Myeloid cell     T cell   Tumor cell 
48           37            1            7           22          322          249          802 

STGG1.fov22.cell.frac<-data.frame(table(STGG1.fov22.obj@meta.data$celltype.new))
STGG1.fov22.cell.frac<-STGG1.fov22.cell.frac[-c(3,5,8),]
STGG1.fov22.cell.frac$FOV<-"T-PN"

p16<-ggplot(STGG1.fov22.cell.frac)+geom_col(aes(x=Celltype,y=STGG1.fov22.cell.frac$`Number of cell`,fill=Celltype),width = 0.8)+
  mytheme+coord_flip()+ylab("Number of cells")+NoLegend()+theme(plot.margin = margin(0,0,0,-1))+
  scale_fill_manual(values = c("white","magenta","red","cyan","#AC844B"))

svg("Figure5C_STGG1_fov22_with_cell_fraction.svg",width = 8,height = 6)
p13+p14+plot_layout(widths = c(5, 1))
dev.off()


#############################################################Figure 5E
STGG1_fov5.id<-nano.obj.with.label.transfer@images$STGG1_fov5@boundaries$centroids@cells
STGG1_fov5.obj<-subset_opt(nano.obj.with.label.transfer,cells = STGG1_fov5.id)

dotplot.obj<-subset_opt(nano.obj.with.label.transfer,cells = c(STGG1_fov5.id,STGG1_fov22.id,STGU1_fov16.id))

DefaultAssay(dotplot.obj)<-"SCT"

Idents(dotplot.obj)<-"celltype.new"
type<-as.character(as.vector(dotplot.obj@meta.data$fov))
type[type%in%"5"]<-"T-TLS"
type[type%in%"22"]<-"T-PN (Tumor)"
type[type%in%"16"]<-"T-PN (Meninges)"
dotplot.obj@meta.data$type<-type


###MMP2 dotplot size by number
heatmap_legend_params <- list(title = "Scaled exp",
                              legend_width = unit(3, "cm"),
                              title_position = "topcenter")

HeatmapDotPlot1 = function(colour, size,colour.label="Scaled exp", size.label="Fraction expressing", cell.size = 0.5, scale=TRUE, cluster_columns=FALSE, cluster_rows=FALSE, col=c("blue","red"), show_heatmap_legend=TRUE, row_title_rot = 0, ...){
  
  if(scale){
    colour = colour %>% t() %>% scale() %>% t()
  }
  if(is.character(col)){
    if(scale){
      col_fun = colorRamp2(seq(from=-3, to=3, length.out=length(col)), col)
    } else {
      col_fun = colorRamp2(seq(from=min(colour), to=max(colour), length.out=length(col)), col)
    }
  } else {
    col_fun = col
  }
  hm = Heatmap(matrix = colour,
               row_gap = unit(2, "mm"), 
               column_gap = unit(2, "mm"), 
               cell_fun = function(j, i, x, y, width, height, fill){
                 grid.rect(x = x, y = y, width = width, height = height, 
                           gp = gpar(col = "grey", fill = NA))
                 grid.circle(x = x, y = y, r = 0.5*(size[i,j])*max(unit.c(width, height)), 
                             gp = gpar(fill = col_fun(colour[i, j]), col = NA))
               },
               rect_gp = gpar(type="none"), height = unit(nrow(colour)*cell.size, "cm"), width = unit(ncol(colour)*cell.size, "cm"),
               cluster_columns=cluster_columns, cluster_rows=cluster_rows, show_heatmap_legend = TRUE, row_title_rot = row_title_rot, ...)
  hm.legend = list()
  if(!is.null(size.label)){
    hm.legend = c(hm.legend, list(Legend(title = size.label1,
                                         labels = c(0.25, 0.50, 0.75, 1.00) %>% as.character,
                                         size=unit.c(unit(sqrt(0.25)*cell.size, "cm"),
                                                     unit(sqrt(0.5)*cell.size, "cm"),
                                                     unit(sqrt(0.75)*cell.size , "cm"),
                                                     unit(sqrt(1.0)*cell.size, "cm")),
                                         type = "points",
                                         grid_height = unit(cell.size,"cm"),
                                         grid_width=unit(cell.size,"cm"),
                                         legend_height=unit(4*cell.size*2, "cm"),
                                         background = NA)))
  }
  if(!is.null(colour.label)){
    hm.legend = c(hm.legend, list(Legend(title=colour.label,
                                         col_fun = col_fun)))
  }
  return(hm) 
}


aggr.by<-c("celltype.new","type")
aggr.fun<-"mean"
features<-c("MMP2")

data = GetAssayData(MMP2.plot.obj, slot="data", assay="SCT")
data = t(data.frame(data[intersect(features, rownames(data)), ]))
rownames(data)<-"MMP2"
aggr = MMP2.plot.obj@meta.data[,aggr.by,drop=F]
aggr = apply(aggr, 1, paste, collapse="_____")
aggr.levels = levels(factor(MMP2.plot.obj@meta.data[,aggr.by[1]]))
if(length(aggr.by) > 1){
  for(a in aggr.by[-1]){
    aggr.levels = unlist(lapply(aggr.levels, function(aggr.level){paste(aggr.level, levels(factor(MMP2.plot.obj@meta.data[,a])), sep='_____')}))
  }
}
aggr = factor(aggr, levels = aggr.levels)
meta.data = MMP2.plot.obj@meta.data

aggr = factor(aggr, levels=intersect(levels(aggr), unique(as.character(aggr))))
data = lapply(levels(aggr), function(x){data[,aggr==x,drop=F]}) %>% setNames(levels(aggr))
colour = do.call(cbind, lapply(data, function(data){apply(data,1,aggr.fun)})) # average expression
size = do.call(cbind, lapply(data, function(data){apply(as.matrix(data) > min(data),1, sum)})) # number of cells expressed, assuming min(data) corresponds to 0 counts
size = as.matrix(size)/max(size)
size = sqrt(size)
meta.data = meta.data[match(names(data), aggr),,drop=F]
colnames(colour) = gsub("_____.*$", "", colnames(colour))
colnames(size) = gsub("_____.*$", "", colnames(size))

split.by="type"
if(length(split.by) == 1){
  grouping = meta.data[,split.by]
}

gene_grouping = NULL
annot.columns=NULL
annot.colours = NULL
show_legend=TRUE
show_annotation_name=FALSE
annotation_labels = NULL
legend_title = NULL
size.label="Number of cells expressed"
cell.size = 1



if(length(annot.columns) > 0){
  if(length(annot.colours) > 0){
    annot = HeatmapAnnotation(df = meta.data[,annot.columns,drop=F], col = annot.colours, which="column", show_legend = show_legend, show_annotation_name = show_annotation_name)
  } else {
    annot = HeatmapAnnotation(df = meta.data[,annot.columns,drop=F], which="column", show_legend=show_legend, show_annotation_name = show_annotation_name)
  }
  hm = HeatmapDotPlot(colour=colour,size=size,col=c("blue", "grey", "red"), cell.size=cell.size, top_annotation = annot, show_heatmap_legend=show_legend, column_split=grouping, row_split = gene_grouping
  )
} else {
  hm = HeatmapDotPlot(colour=colour, size=size, col=c("blue", "grey", "red"), cell.size=cell.size, column_split=grouping, row_split = gene_grouping
  )
}


lgd_list = list(Legend(title = size.label,
                       labels = c(43, 86, 129, 172) %>% as.character,
                       size=unit.c(unit(sqrt(0.25)*1.33, "cm"),
                                   unit(sqrt(0.5)*1.33, "cm"),
                                   unit(sqrt(0.75)*1.33 , "cm"),
                                   unit(sqrt(1.0)*1.33, "cm")),
                       type = "points",
                       grid_height = unit(cell.size,"cm"),
                       grid_width=unit(2,"cm"),
                       legend_height=unit(4*cell.size*2, "cm"),
                       background = NA,
                       legend_width = 4))



svg("Figure5_MMP2_expression_dotplot_size_by_number.svg",width = 14,height=6)
draw(hm, annotation_legend_list = lgd_list)
dev.off()



########COL4A1
feature1<-c("COL4A1")

data1 = GetAssayData(MMP2.plot.obj, slot="data", assay="SCT")
data1 = t(data.frame(data1[intersect(feature1, rownames(data1)), ]))
rownames(data1)<-"COL4A1"

meta.data1 = MMP2.plot.obj@meta.data

data1 = lapply(levels(aggr), function(x){data1[,aggr==x,drop=F]}) %>% setNames(levels(aggr))

colour1 = do.call(cbind, lapply(data1, function(data1){apply(data1,1,aggr.fun)})) # average expression

size1 = do.call(cbind, lapply(data1, function(data1){apply(as.matrix(data1) > min(data1),1, sum)})) # number of cells expressed, assuming min(data) corresponds to 0 counts
size1 = as.matrix(size1)/max(size1)
size1 = sqrt(size1)

meta.data1 = meta.data1[match(names(data1), aggr),,drop=F]
colnames(colour1) = gsub("_____.*$", "", colnames(colour1))
colnames(size1) = gsub("_____.*$", "", colnames(size1))

if(length(split.by) == 1){
  grouping1 = meta.data1[,split.by]
}

if(length(annot.columns) > 0){
  if(length(annot.colours) > 0){
    annot1 = HeatmapAnnotation(df = meta.data1[,annot.columns,drop=F], col = annot.colours, which="column", show_legend = show_legend, show_annotation_name = show_annotation_name)
  } else {
    annot1 = HeatmapAnnotation(df = meta.data1[,annot.columns,drop=F], which="column", show_legend=show_legend, show_annotation_name = show_annotation_name)
  }
  hm1 = HeatmapDotPlot(colour=colour1,size=size1,col=c("blue", "grey", "red"), cell.size=cell.size, top_annotation = annot1, show_heatmap_legend=show_legend, column_split=grouping, row_split = gene_grouping
  )
} else {
  hm1 = HeatmapDotPlot(colour=colour1, size=size1, col=c("blue", "grey", "red"), cell.size=cell.size, column_split=grouping, row_split = gene_grouping
  )
}

hm1

lgd_list1 = list(Legend(title = size.label,
                       labels = c(42, 83, 125, 166) %>% as.character,
                       size=unit.c(unit(sqrt(0.25)*1.33, "cm"),
                                   unit(sqrt(0.5)*1.33, "cm"),
                                   unit(sqrt(0.75)*1.33 , "cm"),
                                   unit(sqrt(1.0)*1.33, "cm")),
                       type = "points",
                       grid_height = unit(cell.size,"cm"),
                       grid_width=unit(2,"cm"),
                       legend_height=unit(4*cell.size*2, "cm"),
                       background = NA,
                       legend_width = 4))

svg("Figure5_COL4A1_expression_dotplot_size_by_number.svg",width = 14,height=6)
draw(hm1, annotation_legend_list = lgd_list1)
dev.off()


#############################################################Figure 5F
STGU1_fov16.obj<-subset_opt(nano.obj.with.label.transfer, cells =STGU1_fov16.id)

Idents(STGU1.fov16.obj)<-"LT_full_gene"

STGU1.fov16.meta.data<-data.frame(STGU1_fov16.obj@meta.data)
STGU1.fov16.meta.data.f<-STGU1.fov16.meta.data[,c(10,11)]
STGU1.fov16.meta.data.f$celltype<-as.vector(Idents(STGU1_fov16.obj))
STGU1.fov16.meta.data.f<-cbind(paste0("X",rownames(STGU1.fov16.meta.data.f)),STGU1.fov16.meta.data.f)
colnames(STGU1.fov16.meta.data.f)<-c("cell","x","y","celltype")

STGU1.fov16.counts<-data.frame(STGU1_fov16.obj@assays$SCT@data)

# create SpaTalk data
STGU1.fov16.SpaT.obj <- createSpaTalk(st_data = as.matrix(STGU1.fov16.counts),
                     st_meta = STGU1.fov16.meta.data.f[,-4],
                     species = "Human",
                     if_st_is_sc = T,
                     spot_max_cell = 1,
                     celltype = STGU1.fov16.meta.data.f$celltype)

# Filter LRIs with downstream targets
STGU1.fov16.SpaT.obj <- find_lr_path(object = STGU1.fov16.SpaT.obj, lrpairs = lrpairs, pathways = pathways)

# Infer cell-cell communications from all celltype pairs
STGU1.fov16.SpaT.obj <- dec_cci_all(object = STGU1.fov16.SpaT.obj,use_n_cores = 10)
save(obj,file = "STGU1_fov16_cell-cell_interatcion.RData")

STGU1.fov16.LR_pair<-STGU1.fov16.SpaT.obj@lrpair
write.xlsx(STGU1.fov16.LR_pair,file = "STGU1_fov16_cell-cell_interaction.xlsx",rowNames=T,colNames=T)


svg("Figure5F_STGU1_fov16_MMP2_interaction.svg",width = 7,height = 4)
plot_lrpair(object = obj,
            celltype_sender = 'Fibroblast',
            ligand = 'MMP2',
            celltype_receiver = 'Endothelial',
            receptor = 'PECAM1',
            if_plot_density = F,
            size = 1.5,
            color = c("#619CFF","red","lightgray"),
            arrow_length = 0.05)
dev.off()


#############################################################Figure 5G
STGG1.fov5.obj<-subset_opt(nano.obj.with.label.transfer, cells =STGG1.fov5.id)

Idents(STGG1.fov5.obj)<-"LT_full_gene"

STGG1.fov5.meta.data<-data.frame(STGG1.fov5.obj@meta.data)
STGG1.fov5.meta.data.f<-STGG1.fov5.meta.data[,c(10,11)]
STGG1.fov5.meta.data.f$celltype<-as.vector(Idents(STGG1.fov5.obj))
STGG1.fov5.meta.data.f<-cbind(paste0("X",rownames(STGG1.fov5.meta.data.f)),STGG1.fov5.meta.data.f)
colnames(STGG1.fov5.meta.data.f)<-c("cell","x","y","celltype")

STGG1.fov5.counts<-data.frame(STGG1.fov5.obj@assays$SCT@data)

# create SpaTalk data
STGG1.fov5.SpaT.obj <- createSpaTalk(st_data = as.matrix(STGG1.fov5.counts),
                     st_meta =STGG1.fov5.meta.data.f[,-4],
                     species = "Human",
                     if_st_is_sc = T,
                     spot_max_cell = 1,
                     celltype = STGG1.fov5.meta.data.f$celltype)

# Filter LRIs with downstream targets
STGG1.fov5.SpaT.obj <- find_lr_path(object = STGG1.fov5.SpaT.obj, lrpairs = lrpairs, pathways = pathways)

# Infer cell-cell communications from all celltype pairs
STGG1.fov5.SpaT.obj <- dec_cci_all(object = STGG1.fov5.SpaT.obj,use_n_cores = 10)
save(STGG1.fov5.SpaT.obj,file = "STGG1_fov5_cell-cell_interatcion.RData")

STGG1.fov5.LR_pair<-STGG1.fov5.SpaT.obj@lrpair
write.xlsx(STGG1.fov5.LR_pair,file = "STGG1_fov5_cell-cell_interaction.xlsx",rowNames=T,colNames=T)

svg("Figure5G_STGG1_fov5_MMP2_interaction.svg",width = 7,height = 4)
plot_lrpair(object = obj,
            celltype_sender = 'Fibroblast',
            ligand = 'MMP2',
            celltype_receiver = 'Endothelial',
            receptor = 'PECAM1',
            if_plot_density = F,
            size = 1.5,
            color = c("#619CFF","red","lightgray"),
            arrow_length = 0.05)
dev.off()
