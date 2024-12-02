if(1){
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(dplyr)
  library(KernSmooth)
  library(RColorBrewer)
  library(plotly)
  library(BiocParallel)
  library(grid)
  library(ComplexHeatmap)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(ggforce)
}


saple_EP_T <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_T_subcluster.Rds")   
sce_SCENIC <- open_loom("E:/LUAD_CTS/Result3/data/SCENIC/output/out_SCENIC.loom") 

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
 
# Calculate cell-specific transcription factors
cellinfo <- saple_EP_T@meta.data[,c('cell_state','stage')]#细胞meta信息
cellinfo$cell_state<-factor(cellinfo$cell_state,
                             levels =c("NEF","Immune","Cellcycle","EMT"))

cellinfo <-cellinfo[order(cellinfo$cell_state),] 
cell_state <-  as.data.frame(subset(cellinfo,select = 'cell_state'))
selectedResolution <- "cell_state"
sub_regulonAUC <- regulonAUC[,rownames(cell_state)]

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cell_state[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)
rssPlot <-  plotRSS(
    rss, 
    zThreshold = 1.5 ,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.01,
    varName = "cellType" )

pdf(file = "E:/LUAD_CTS/Result3/Figure4D.pdf")
rssPlot
dev.off()
#  Calculate AUCell activity
auc_active = getAUC(sub_regulonAUC)

auc_active <- auc_active[rev(unique(c("HMGA1(+)",rssPlot$rowOrder ))),rownames(cell_state) ]

ann_colors = list( cell_state = c("NEF" = "#82B0D2","Immune" = "#8ECFC9", 
                                  "Cellcycle" = "#FA7F6F",  "EMT"= "#FFBE7A"   ))
pdf(file = "../results_figures/AUC_score.pdf")
pheatmap( auc_active  , show_colnames=F, cluster_cols = F,cluster_rows = F, annotation_colors = ann_colors,
                    annotation_col=cell_state , treeheight_row=10, scale = "row", breaks=seq(-4, 4, length.out = 100),
                    treeheight_col=10,colorRampPalette(colors = c("lightblue","white","red"))(100) )
dev.off()

TF_name =c("NFIX","NKX2-1","FOXA2","ELF5","CEBPA","NFATC2","ELF4","ETS1","IRF4",
           "RUNX3","STAT4","ZNF816","IKZF1","E2F6","PRDM1","IKZF3","STAT5A","NFYB", 
           "TFDP1","TEAD4","E2F1","ZNF467","FOSL1","E2F8","MYBL1","BARX1","GABPB1",
           "MXD3","FOXM1","HNF4A","SMAD3","HMGA1" )
TF_exp <- saple_EP_T@assays[["RNA"]]@data[TF_name,]
TF_exp <-TF_exp[,rownames(cell_state)]
pheatmap::pheatmap( TF_exp  , show_colnames=F, cluster_cols = F,cluster_rows = F,
                    annotation_col=cell_state , treeheight_row=10, scale = "row", breaks=seq(-4, 4, length.out = 100),
                    treeheight_col=10,colorRampPalette(colors = c("lightblue","white","red"))(100) )
