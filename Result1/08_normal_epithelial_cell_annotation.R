# Load libraries 
library(future)  
library(Seurat)
library(SingleR)
library(dplyr) 
library(ggplot2)
library(scater) 

# Increase memory limits 
options(future.globals.maxSize = 50 * 1024^3)  
memory.limit(size = 900000000 ) 

# Normal epithelial cell annotation
saple_EP <- readRDS(file = "../results_data/saple_EP.rds")
saple_EP_N <- subset(saple_EP, T_N_cell == "Normal")
saple_EP_N$seurat_clusters <- as.character(saple_EP_N$seurat_clusters)
vaildmaker  <- list( 'NE'= c("CALCA","CHGA", "ASCL1"), 
                     'AT2' = c("SFTPB" ,"SFTPC","SFTPD", "SFTPA1")  , 
                     'AT1'  = c("AGER", "PDPN","CLIC5","CAV1"),  
                     'Club' = c("SCGB3A2", "SCGB3A1","SCGB1A1","MUC5B"),  
                     'Ciliated'  = c("TUBA1A","FOXJ1", "CAPS", "CCDC78") ,
                     'Tumor'  = c("MDK","TIMP1")
                      )   

# Figure 1G
pdf("../results_figures/marker_exp_in_clusters.pdf")
DotPlot(object =saple_EP_N , features = vaildmaker,  assay = "RNA" ) +theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("lightgrey", "#ff0000"))
dev.off()


# singleR cell annotation
Ref_sce   <- readRDS("../../results_1/results_data/lung_58celltype.rds")
saple06_annot  <-  saple_EP_N@assays[["RNA"]]@counts 
dim(saple06_annot)
common_gene <- intersect(rownames(Ref_sce), rownames(saple06_annot))
Ref_sce <- Ref_sce[common_gene,]
Ref_sce <- Ref_sce[common_gene,c(which(Ref_sce@colData@listData[["free_annotation"]]=="Club"),which(Ref_sce@colData@listData[["free_annotation"]]=="Ciliated"),
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Proximal Ciliated"),which(Ref_sce@colData@listData[["free_annotation"]]=="Basal"),
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Proximal Basal"), 
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Proliferating Basal"),which(Ref_sce@colData@listData[["free_annotation"]]=="Goblet"),
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Mucous"),which(Ref_sce@colData@listData[["free_annotation"]]=="Serous"),
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Ionocyte"),which(Ref_sce@colData@listData[["free_annotation"]]=="Neuroendocrine"),
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Alveolar Epithelial Type 1"),which(Ref_sce@colData@listData[["free_annotation"]]=="Alveolar Epithelial Type 2"),
                                 which(Ref_sce@colData@listData[["free_annotation"]]=="Signaling Alveolar Epithelial Type 2")
)]
saple06_annot <- saple06_annot[common_gene,]
saple06_annot <- SummarizedExperiment(assays=list(counts=saple06_annot)) %>% logNormCounts()
pred <- SingleR(test=saple06_annot, ref=Ref_sce,
                labels=Ref_sce$free_annotation,
                clusters = saple_EP_N@meta.data[["seurat_clusters"]])
pdf("../results_figures/singleR_celltype.pdf" )
plotScoreHeatmap(pred,show_colnames =T,cluster_cols = T ,show.labels = F) 
dev.off()

# Add normal celltype
saple_EP@meta.data$T_N_cells <- saple_EP@meta.data$T_N_cell
saple_EP@meta.data$T_N_cells <- as.character(saple_EP@meta.data$T_N_cells)
saple_EP$T_N_cells[which(saple_EP$seurat_clusters == '2')] =  "AT2"
saple_EP$T_N_cells[which(saple_EP$seurat_clusters == '13')] =  "AT2"
saple_EP$T_N_cells[which(saple_EP$seurat_clusters == '4')] =  "Ciliated"
saple_EP$T_N_cells[which(saple_EP$seurat_clusters == '8')] =  "Club"
saple_EP$T_N_cells[which(saple_EP$seurat_clusters == '10')] =  "AT1"
saple_EP$T_N_cells[which(saple_EP$seurat_clusters == '22')] =  "NE"
saple_EP@meta.data$T_N_cells <- factor(saple_EP@meta.data$T_N_cells,levels = c("Tumor","NE","AT2","AT1","Club", "Ciliated"))

# Save Seurat object
saveRDS(saple_EP,file = "../results_data/saple_EP_anno.rds")

# Figure 1F
pdf("../results_figures/umap_celltype.pdf")
DimPlot(saple_EP, reduction = "umap",group.by = "T_N_cells",label = T,repel = T,raster = T )+scale_color_npg()+scale_fill_npg()+NoLegend() 
dev.off()




