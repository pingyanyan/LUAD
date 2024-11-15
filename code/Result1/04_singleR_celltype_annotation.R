# Load libraries 
library(future)  
library(scater) 
library(SingleR) 
library(dplyr)

# Increase memory limits
memory.limit(size = 900000000 )
options(future.globals.maxSize = 50 * 1024^3)  
 
# Load the Seurat object
saple_obj_har_after  <- readRDS("../results_data/saple_obj_anno_har_after.rds")
Ref_sce <- readRDS("../results_data/lung_58celltype.rds")

# Create a SummarizedExperiment object 
saple_counts  <-  saple_obj_har_after@assays[["RNA"]]@counts  
common_gene <- intersect(rownames(Ref_sce), rownames(saple_counts))
Ref_sce <- Ref_sce[common_gene,]
saple_counts <- saple_counts[common_gene,]
saple_counts <- SummarizedExperiment(assays=list(counts=saple_counts)) %>% 
  logNormCounts()

# Perform cell type annotation using SingleR
pred <- SingleR(test=saple_counts, ref=Ref_sce,
                labels=Ref_sce$celltype,
                clusters = saple_obj_har_after@meta.data[["seurat_clusters"]])

# Save the heatmap plot
pdf("../results_figures/singleR_celltype.pdf" )
plotScoreHeatmap(pred,show_colnames =T,cluster_cols = T,show.labels = F)
dev.off() 


