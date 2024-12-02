# Load libraries 
library(future)  
library(Seurat) 
library(SingleR)
library(ggplot2)
library(dplyr)

# Increase memory limits
memory.limit(size = 900000000 )
options(future.globals.maxSize = 50 * 1024^3)  

# Load the Seurat object 
saple_obj_har_after <- readRDS("../results_data/saple_obj_anno_har_after.rds")

# Add stage 
saple_obj_har_after$stage <- "Normal"
saple_obj_har_after$stage[which(saple_obj_har_after$Stage == 'AIS'|saple_obj_har_after$Stage == 'MIA'| 
                                    saple_obj_har_after$Stage == 'IAC'| saple_obj_har_after$Stage == 'IB' ) ] = "early"
saple_obj_har_after$stage[which(saple_obj_har_after$Stage == 'IIA'|saple_obj_har_after$Stage == 'IIB'|
                                      saple_obj_har_after$Stage == 'IIIA')] = "advance"

# Marker genes
vaildmaker  <- list( 'Epithelial cells' = c("EPCAM","MUC1", "KRT19","KRT18","CDH1","KRT8")  , 
                     'T cells'  = c("CD3D" ,"CD3E", "CD3G" ,"CD8A"),  
                     'NK cells'  = c("NKG7","KLRD1","GNLY","KLRB1"),  
                     'B cells' = c("CD79A" ,"CD19","JCHAIN","MS4A1","IGHG1","IGHA1"),  
                     'Myeloid cells' = c("CD68", "LYZ","MARCO",  "AIF1","CTSD") ,  
                     'Mast cells' = c("TPSB2","MS4A2" ,  "TPSAB1","CPA3"),
                     'Endothelial cells'= c("PECAM1","VWF","CDH5" , "CLDN5", "RAMP2"),
                     'Fibroblasts cells'= c("COL1A1","PDGFRA" , "COL1A2" ,"ACTA2","DCN","LUM"))   

# Save the Dot Plot
pdf("../results_figures/marker_exp_in_clusters.pdf", width =12 , height=6 )
DotPlot(object =saple_obj_har_after , features = vaildmaker, group.by = "celltype",
        assay = "RNA" ) +theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("lightgrey", "#ff0000"))
dev.off()

# Assign cell types to clusters based on predefined markers
celltype <- data.frame(clusterID = 0:29,celltype =c("T cells","NK cells","Myeloid cells","Myeloid cells","Epithelial cells","Myeloid cells", "B cells", "Mast cells","Endothelial cells","Myeloid cells","B cells", 
                                                      "Epithelial cells", "Fibroblasts cells","Myeloid cells", "Epithelial cells","Myeloid cells","Myeloid cells", "Epithelial cells","Epithelial cells","Epithelial cells","Fibroblasts cells",
                                                      "Epithelial cells", "Mast cells","Epithelial cells","Myeloid cells","B cells", "Endothelial cells","T cells","Epithelial cells", "Myeloid cells"))
for (i in 1:dim(saple_obj_har_after)[2]) {
  index = saple_obj_har_after@meta.data[["seurat_clusters"]][i]
  saple_obj_har_after@meta.data$celltype[i] = celltype[index,2]
}
saple_obj_har_after$celltype = factor(saple_obj_har_after$celltype,levels =c("Epithelial cells",'T cells',"NK cells",'B cells', 
                                                                             'Myeloid cells','Mast cells','Endothelial cells','Fibroblasts cells'  ) )

# Save Seurat object
saveRDS(saple_obj_har_after,file = "../results_data/saple_obj_anno_har_after.rds")


# Save figures as PDF
pdf("../results_figures/marker_exp_in_celltype.pdf", width =12 , height=10)
 DotPlot(object =saple_obj_har_after , features = vaildmaker, group.by = "celltype",
         assay = "RNA" ) +theme_bw()+
   theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
   labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
   scale_color_gradientn(values = seq(0,1,0.2),colours = c("lightgrey", "#ff0000")) 
dev.off()

pdf("../results_figures/marker_exp_in_celltype_umap.pdf", width =12 , height=8)
FeaturePlot(saple_obj_har_after,features =  c("MUC1","CD3D","NKG7","CD79A", "LYZ","TPSAB1","CLDN5","DCN"),
            reduction = "umap",ncol = 4 , cols = c("lightgrey", "#ff0000" ) )    
dev.off() 

pdf("../results_figures/celltype_umap.pdf")
DimPlot(saple_obj_har_after, reduction = "umap",group.by = "celltype",label = T,repel = T )+scale_color_npg()+scale_fill_npg()+ ggplot2::theme(legend.position = "none")
dev.off() 
