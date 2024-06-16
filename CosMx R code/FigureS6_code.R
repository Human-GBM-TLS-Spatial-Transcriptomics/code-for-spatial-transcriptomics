####################################################################################################################################################################################
###################################################################### code for figureS6 ########################################################################################### 
####################################################################################################################################################################################

#Figure S6A
svg("FigureS6A_P240_fov22_CCL19(EC)_CCR7(CD4)_interaction.svg",width = 7,height = 4)
plot_lrpair(object = P240.fov22.SpaT.obj,
            celltype_sender = 'Myeloid_cell',
            ligand = 'CCL21',
            celltype_receiver = 'CD4.T.cells',
            receptor = 'CCR7',
            if_plot_density = F,
            size = 1.5,
            color = c("#AC844B","#08A045","lightgray"),
            arrow_length = 0.05)
dev.off()

#FigureS6B
svg("FigureS6B_P240_fov22_CCL19(EC)_CCR7(CD4)_downstream.svg",width = 7,height = 4)
plot_lr_path(object = P240.fov5.SpaT.obj,                
             celltype_sender = 'Endothelial',
             ligand = 'CCL19',
             celltype_receiver = 'CD4.T.cells',
             receptor = 'CCR7',
             color = c("red","#08A045"))+labs(title = "Downstream targets and TFs")
dev.off()

#FigureS6C
svg("FigureS6C_P240_fov22_CCL21(Mye)_CCR7(CD4)_downstream.svg",width = 7,height = 4)
plot_lr_path(obj12,
             celltype_sender = c('Myeloid_cell'),
             ligand = 'CCL21',
             celltype_receiver = 'CD4.T.cells',
             receptor = 'CCR7',
             color = c("#AC844B","#08A045"))+labs(title = "Downstream targets and TFs")
dev.off()
