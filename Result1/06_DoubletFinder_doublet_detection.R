# Load libraries 
  library(future)  
  library(Seurat) 
  library(DoubletFinder)
  library(dplyr)
  
# Increase memory limits
  memory.limit(size = 900000000 )
  options(future.globals.maxSize = 50 * 1024^3)  

# Load the Seurat object
saple_obj_anno_har_after  <- readRDS("../results_data/saple_obj_anno_har_after.rds")

# Define the number of principal components
dim.usage = "30"

# Function to perform doublet detection using DoubletFinder
Find_doublet <- function(data){
  # Perform parameter sweep for pK values using DoubletFinder
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Estimate homotypic doublet proportion based on cell types
  homotypic.prop <- modelHomotypic(data$celltype)
  
  # Estimate expected doublet rate
  DoubletRate = ncol(data)*8*1e-6
  nExp_poi <- round(DoubletRate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # Select the best pK value based on the highest mean BC score
  p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  
  # Run DoubletFinder using the selected pK and adjusted expected doublet count
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}

# Split the Seurat object by "samples" for separate doublet detection
scRNAlist <- SplitObject( object = saple_obj_anno_har_after  ,split.by ="samples" ) 

# Apply normalization, scaling, PCA, and doublet detection for each sample
for(i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData( scRNAlist[[i]]) %>% FindVariableFeatures( selection.method = "vst") %>%
    ScaleData() %>% RunPCA() %>%  Find_doublet() 
}

saple_obj_har_after_exclud_doule <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])  
saple_obj_anno_har_after$doublet_info  = saple_obj_har_after_exclud_doule$doublet_info

saveRDS(saple_obj_har_after_exclud_doule,file = "../results_data/saple_obj_har_after_exclud_doule.rds")
saveRDS(saple_obj_anno_har_after,file = "../results_data/saple_obj_anno_har_after.rds")
 