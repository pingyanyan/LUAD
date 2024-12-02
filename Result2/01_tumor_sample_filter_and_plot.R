library(Seurat)
library(dplyr)
saple_EP <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_stage.Rds")
saple_EP <-  subset(saple_EP,subset = T_N_cell  == "Normal" |T_N_cell  == "Tumor" &  tiss_origin  == "Tumor_origin" )   

# Extract patients with more than 200 tumor cells
saple_EP_T <- subset(saple_EP,subset = T_N_cells  == "Tumor")   
saple_EP_T <- subset(saple_EP_T,subset = samples  != "GSM79"& samples  != "GSM80" 
                     &samples  != "GSM81"&samples  != "GSM83"  &samples  != "GSM84" )
saple_EP_T <- NormalizeData(saple_EP_T) %>%  FindVariableFeatures() %>% ScaleData()  %>%  RunPCA()
saple_EP_T <- FindNeighbors(saple_EP_T,  dims=1:25) %>% FindClusters(resolution = 0.5)  
saple_EP_T <- saple_EP_T %>%RunUMAP(dims  = 1:25)

pdf(file = "E:/LUAD_CTS/Result2/Figure2A.pdf")
DimPlot(saple_EP_T, reduction = "umap", group.by = "samples", label = TRUE ,repel = T  )
dev.off()

saveRDS(saple_EP_T,file = "E:/LUAD_CTS/Result2/data/saple_EP_T.Rds")
