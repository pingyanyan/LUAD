# Load libraries
library(future)  
library(Seurat)
library(dplyr)
library(harmony)
cols <- c( "#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF" )

# Increase memory limits 
memory.limit(size = 900000000 )
options(future.globals.maxSize = 50 * 1024^3)  

# Load the Seurat object
saple_obj <- readRDS("E:/LUAD_CTS/Result1/data/saple_obj_har_before.rds") 

# Run Harmony to remove batch effects, grouping by the 'samples' variable
saple_obj <- RunHarmony(saple_obj, group.by.vars = "samples")
saple_obj <- FindNeighbors(saple_obj,reduction = "harmony", dims=1:25) %>% 
  FindClusters(resolution = 0.5) 
saple_obj <- saple_obj %>% RunUMAP(reduction = "harmony",dims=1:25)  


# Save the Seurat object after Harmony integration
saveRDS(saple_obj,file = "E:/LUAD_CTS/Result1/data/saple_obj_har_after.rds")

#save figures
pdf("E:/LUAD_CTS/Result1/figure/seurat_clusters_har_after.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "seurat_clusters",label = T,repel = T ,cols =  colp)+ ggplot2::theme(legend.position = "none")
dev.off()

# Save UMAP plots as PDF figures
pdf("E:/LUAD_CTS/Result1/figure/samples_har_after.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "samples" ,cols =  cols)
dev.off()

pdf("E:/LUAD_CTS/Result1/figure/Phase_har_after.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "Phase" )
dev.off()

pdf("E:/LUAD_CTS/Result1/figure/orig.ident_har_after.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "orig.ident" )+scale_color_jco()+scale_fill_jco()
dev.off()

 

