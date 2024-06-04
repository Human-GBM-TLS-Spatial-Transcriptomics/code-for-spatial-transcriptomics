####################################################################################################################################################################################
###################################################################### code for figure7S ########################################################################################### 
####################################################################################################################################################################################

#Figure S7A
STGU3.fov2.id<-as.vector(nano.obj.with.label.transfer@images$STGU3_fov2@boundaries$centroids@cells)
STGU3.fov2.obj<-subset_opt(nano.obj.with.label.transfer, cells =STGU3.fov2.id)

table(STGU3.fov2.obj@meta.data$celltype.new)
B cell  Endothelial  Erythrocyte   Fibroblast   Mural cell Myeloid cell   Neutrophil 
69           43            7          238           42          247            1 
T cell   Tumor cell 
257          742 


STGU3.fov2.cell.frac<-data.frame(table(STGU3.fov2.obj@meta.data$celltype.new))
STGU3.fov2.cell.frac$FOV<-"T-TLS"
STGU3.fov2.cell.frac<-STGU3.fov2.cell.frac[-c(3,5,7,9),]

p19<-ImageDimPlot(nano.obj.with.label.transfer, fov = "STGU3_fov2",size = 0.5, coord.fixed = FALSE,axes = T,
                  cells = WhichCells(nano.obj.with.label.transfer,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                  cols=c("red","#AC844B","cyan","magenta","white"),border.color = "black")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))

colnames(STGU3.fov2.cell.frac) = c('Celltype','Number of cell','Group')
STGU3.fov2.cell.frac<-STGU3.fov2.cell.frac[order(STGU3.fov2.cell.frac$`Number of cell`,decreasing = F),]
STGU3.fov2.cell.frac$Celltype<-factor(STGU3.fov2.cell.frac$Celltype,levels = STGU3.fov2.cell.frac$Celltype)

p20<-ggplot(STGU3.fov2.cell.frac)+geom_col(aes(x=Celltype,y=STGU3.fov2.cell.frac$`Number of cell`,fill=Celltype),width = 0.8)+
  mytheme+coord_flip()+ylab("Number of cells")+NoLegend()+theme(plot.margin = margin(0,0,0,-1))+
  scale_fill_manual(values = c("red","magenta","white","#AC844B","cyan"))

svg("FigureS7A_STGU3_fov2_with_cell_fraction.svg",width = 8,height = 6)
p19+p20+plot_layout(widths = c(5, 1))
dev.off()


STGU1.fov26.id<-as.vector(nano.obj.with.label.transfer@images$STGU1_fov26@boundaries$centroids@cells)
STGU1.fov26.obj<-subset_opt(nano.obj.with.label.transfer, cells =STGU1.fov26.id)

table(STGU1.fov26.obj@meta.data$celltype.new)
B cell  Endothelial  Erythrocyte   Fibroblast   Mural cell Myeloid cell   Neutrophil 
247           57           63           47           41          445           17 
T cell   Tumor cell 
244         1160 


STGU1.fov26.cell.frac<-data.frame(table(STGU1.fov26.obj@meta.data$celltype.new))
STGU1.fov26.cell.frac$FOV<-"M-TLS"
STGU1.fov26.cell.frac<-STGU1.fov26.cell.frac[-c(3,5,7,9),]

p15<-ImageDimPlot(nano.obj.with.label.transfer, fov = "STGU1_fov26",size = 0.5, coord.fixed = FALSE,axes = T,
                  cells = WhichCells(nano.obj.with.label.transfer,idents = c("T cell","Endothelial","B cell","Myeloid cell","Fibroblast")),
                  cols=c("red","#AC844B","cyan","magenta","white"),border.color = "black")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))

colnames(STGU1.fov26.cell.frac) = c('Celltype','Number of cell','Group')
STGU1.fov26.cell.frac<-STGU1.fov26.cell.frac[order(STGU1.fov26.cell.frac$`Number of cell`,decreasing = F),]
STGU1.fov26.cell.frac$Celltype<-factor(STGU1.fov26.cell.frac$Celltype,levels = STGU1.fov26.cell.frac$Celltype)


p16<-ggplot(STGU1.fov26.cell.frac)+geom_col(aes(x=Celltype,y=STGU1.fov26.cell.frac$`Number of cell`,fill=Celltype),width = 0.8)+
  mytheme+coord_flip()+ylab("Number of cells")+NoLegend()+theme(plot.margin = margin(0,0,0,-1))+
  scale_fill_manual(values = c("white","red","cyan","magenta","#AC844B"))

svg("FigureS7A_STGU1_fov26_with_cell_fraction.svg",width = 8,height = 6)
p15+p16+plot_layout(widths = c(5, 1))
dev.off()

#Figure S7B
svg("FigureS7B_STGG1_fov12_CD4_trajectory.svg",width=6,height = 4)
ImageDimPlot(nano.obj.with.label.transfer,fov = "STGG1_fov12",
             cols = c("#8A817C","#80ED99","#8A817C","#8A817C","#8A817C","#008000","red","black",
                      "black","black","black","black","black","black","black","black"),
             border.color = "white")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))
dev.off()

svg("FigureS7B_STGU1_fov18_CD4_trajectory.svg",width=6,height = 4)
ImageDimPlot(nano.obj.with.label.transfer,fov = "STGU1_fov18",
             cols = c("#8A817C","#8A817C","#8A817C","red","black",
                      "black","black","black","black","black","black","black","black"),
             border.color = "white")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))
dev.off()
