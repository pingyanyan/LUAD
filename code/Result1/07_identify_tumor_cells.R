# Load libraries 
library(future)  
library(Seurat)
library(SingleR)
library(dplyr) 
library(ggsci)
library(scales)

# Increase memory limits 
options(future.globals.maxSize = 50 * 1024^3)  
memory.limit(size = 900000000 ) 

# Extract and re-cluster epithelial cells
saple_obj_har_after  <- readRDS("../../results_1/results_data/saple_obj_anno_har_after.rds")
saple_EP <- subset(saple_obj_har_after,subset = celltype == "Epithelial cells")
saple_EP$tiss_origin <- factor(saple_EP$tiss_origin,levels = c("Tumor_origin","Normal_origin"))
saple_EP <- NormalizeData(saple_EP) %>%  FindVariableFeatures() %>% ScaleData()  %>%  RunPCA()
saple_EP <- FindNeighbors(saple_EP,  dims=1:25) %>% FindClusters(resolution = 0.5) 
saple_EP <- saple_EP %>% RunUMAP(dims = 1:25)

# Save UMAP plots 
pdf("../results_figures/umap_seurat_clusters_before.pdf" )
DimPlot(saple_EP, reduction = "umap",label = T,repel = T) 
dev.off()

pdf("../results_figures/umap_doublet_info_before.pdf" )
DimPlot(saple_EP, reduction = "umap",group.by = "doublet_info",label = T,repel = T)+scale_color_lancet()+scale_fill_lancet()
dev.off()

# Filter Doublet
saple_EP <- saple_EP[ ,!(saple_EP$seurat_clusters %in% c("6", "12", "4","14","18", "24" ))]
saple_EP$seurat_clusters <- as.character(saple_EP$seurat_clusters)
saple_EP <- NormalizeData(saple_EP) %>%  FindVariableFeatures() %>% ScaleData()  %>%  RunPCA()
saple_EP <- FindNeighbors(saple_EP,  dims=1:25) %>% FindClusters(resolution = 0.5) 
table(saple_EP@active.ident)   
saple_EP <- saple_EP %>%RunUMAP(dims  = 1:25)

# Figure 1C 
pdf("../results_figures/umap_seurat_clusters_after.pdf" )
DimPlot(saple_EP, reduction = "umap",label = T,repel = T) 
dev.off()

# Figure 1D
pdf("../results_figures/ratio_in_seurat_clusters.pdf" )
Cellratio <-as.data.frame( prop.table(table(saple_EP$tiss_origin,saple_EP$seurat_clusters), margin = 2)) 
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic()+scale_color_jco()+scale_fill_jco() +coord_flip()+
  labs(x='Sample',y = 'Ratio')+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off() 

pdf("../results_figures/ratio_in_seurat_samples.pdf" )
Cellratio <-as.data.frame( prop.table(table(saple_EP$samples,saple_EP$seurat_clusters), margin = 2)) 
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +coord_flip()+
  labs(x='Sample',y = 'Ratio')+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off() 

# Add tumor or normal label
saple_EP$T_N_cell <- 'Tumor'
saple_EP$T_N_cell[which(saple_EP$seurat_clusters == '2'|saple_EP$seurat_clusters=='4'|saple_EP$seurat_clusters=='8'|
                          saple_EP$seurat_clusters=='10'|saple_EP$seurat_clusters=='13'|saple_EP$seurat_clusters=='22')] =  "Normal" 
saple_EP$T_N_cell <- factor(saple_EP$T_N_cell,levels = c("Tumor","Normal"))
saveRDS(saple_EP,file = "../results_data/saple_EP.rds")

