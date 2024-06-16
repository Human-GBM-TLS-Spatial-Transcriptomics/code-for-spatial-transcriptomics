####################################################################################################################################################################################
###################################################################### code for figure7 ####################################################################################### 
####################################################################################################################################################################################
library(ggplot2)
library(viridisLite)
library(Seurat)
library(monocle)
library(reshape2)
library(pheatmap)

#Load CD4 trajectory result (Analysis code can be found at 'Analysis code.R')
cds<-load(file = "CD4_pseudotime_object.RData")

#Figure 7D
pdf("Figure7D_CD4_trajectory_pseudotime.pdf",width = 6,height = 5)
plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.7)
dev.off()

#Figure 7E
pdf("Figure7E_CD4_trajectory_State.pdf",width = 6,height = 5)
plot_cell_trajectory(cds, color_by = "State",cell_size = 0.7)
dev.off()


#Figure 7H
Tumor.TLS.PN.CD4.CosMx@meta.data$CD4_trajctory<-paste0("CD4","_State",Tumor.TLS.PN.CD4.CosMx@meta.data$trajectory.state)

Idents(Tumor.TLS.PN.CD4.CosMx)<-"CD4_trajctory"
Idents(nano.obj.with.label.transfer)<-"LT_full_gene"

Idents(nano.obj.with.label.transfer, cells = colnames(Tumor.TLS.PN.CD4.CosMx)) <- Idents(Tumor.TLS.PN.CD4.CosMx)

table(Idents(nano.obj.with.label.transfer))
CD4_State1   CD4_State7   CD4_State3   CD4_State5   CD4_State4   CD4_State6   CD4_State2  Endothelial   Tumor cell 
1090          650          695          186          763          360            7         5790        97999 
Myeloid cell   Neutrophil   Mural cell       T cell       B cell  Erythrocyte   Fibroblast 
33017          613         6206         4910         4300         3098         5888 

svg("Figure7G_P240_fov22_CD4_trajectory.svg",width=6,height = 4)
ImageDimPlot(nano.obj.with.label.transfer,fov = "P240_fov22",
             cols = c("#8A817C","#80ED99","#8A817C","#8A817C","#8A817C","#008000","red","black",
                      "black","black","black","black","black","black","black","black"),
             border.color = "white")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))
dev.off()

svg("Figure7G_P238_fov2_CD4_trajectory.svg",width=6,height = 4)
ImageDimPlot(nano.obj.with.label.transfer,fov = "P238_fov2",
             cols = c("#8A817C","#80ED99","#8A817C","#8A817C","#8A817C","#008000","red","black",
                      "black","black","black","black","black","black","black","black"),
             border.color = "white")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))
dev.off()

svg("Figure7G_P237_fov26_CD4_trajectory.svg",width=6,height = 4)
ImageDimPlot(nano.obj.with.label.transfer,fov = "P237_fov26",
             cols = c("#8A817C","#80ED99","#8A817C","#8A817C","#8A817C","#008000","red","black",
                      "black","black","black","black","black","black","black","black"),
             border.color = "white")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))
dev.off()

svg("Figure7G_P237_fov8_CD4_trajectory.svg",width=6,height = 4)
ImageDimPlot(nano.obj.with.label.transfer,fov = "P237_fov8",
             cols = c("#8A817C","#80ED99","#8A817C","#8A817C","#8A817C","#008000","red","black",
                      "black","black","black","black","black","black","black","black"),
             border.color = "white")+NoLegend()+
  theme(plot.margin = margin(0,0,0,0),legend.position = "bottom",axis.text=element_text(face="bold", size=7,colour = 'white'))
dev.off()

#Figure 7G
state.fraction<-data.frame(table(cds@phenoData@data$State,cds@phenoData@data$region))
state.fraction$Region<-factor(state.fraction$Region,levels = c("Tumor","T_PN","T_TLS","M_TLS","B_TLS"))
colnames(state.fraction)<-c("State","Region","Number")

state.fraction.table<-dcast(state.fraction,Region~State)
write.xlsx(state.fraction.table,file = "CD4_trajectory_State_fraction.xlsx",rowNames=F,colNames=T)

mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 axis.line.y = element_blank(),  
                 axis.ticks.y = element_blank(), 
                 panel.background = element_rect(fill="white"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank()) 

svg("CD4_Trajectory_State_fraction.svg",width = 6,height=4)
ggplot(state.fraction, aes( x = Region, weight = Number, fill = State))+
  geom_bar( position = "fill")+ylab("Percent")+coord_flip()+ mytheme+
  scale_y_continuous(labels = c("0%","25%","50%","75%","100%"),expand = c(0,0))
dev.off()

#Figure 7K
Tumor.TLS.PN.CD4.CosMx@meta.data$trajectory.state<-cds@phenoData@data$State

State.average.expression<-AverageExpression(Tumor.TLS.PN.CD4.CosMx,features =CD4.gene.list.f,assays = "SCT")
State.average.expression<-data.frame(State.average.expression$SCT)
colnames(State.average.expression)<-c(paste0("State","_",c(1:7)))

svg("CD4_genes_average_expression_by_state(state2 exclued).svg",width = 8,height = 12)
pheatmap(State.average.expression[,-2],cluster_rows = T,cluster_cols = F,scale = "row",fontsize_row = 6,clustering_method = "ward.D2")
dev.off()

write.xlsx(State.average.expression,file = "CD4_Average_expression_of_each_state.xlsx",rowNames=T,colNames=T)


