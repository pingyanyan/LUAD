# Load libraries 
library(future)  
library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)

# Increase memory limits 
memory.limit(size = 900000000 )
options(future.globals.maxSize = 50 * 1024^3)  
 
# Set working directory for the LUAD single-cell data
setwd("F:/LUAD/data/single_cell/")

# GSE189357
assays <- dir("GSE189357_RAW/")
dir <- paste0("GSE189357_RAW/", assays)
samples <- c('GSM77','GSM78','GSM79','GSM80','GSM81','GSM82','GSM83','GSM84',"GSM85") 
Stage <- c("IAC","IAC","MIA","MIA","AIS","MIA","AIS","AIS","IAC")
tiss_origin<- c("Tumor_origin","Tumor_origin","Tumor_origin","Tumor_origin",
                "Tumor_origin","Tumor_origin","Tumor_origin","Tumor_origin","Tumor_origin")
tiss_site <- c(  "Primary",  "Primary",    "Primary",    "Primary",   "Primary",   "Primary", 
                 "Primary",   "Primary",   "Primary")

# Initialize an empty list to store Seurat objects
scRNAlist <- list()

# Loop through each sample, process the data and create Seurat objects
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i]) 
  
  scRNAlist[[i]] <- CreateSeuratObject(counts, project= "GSE189357",min.cells=3 )     
  # Add metadata (sample info)
  scRNAlist[[i]]@meta.data[["samples"]] <- rep(samples[i], dim(counts)[2]) 
  scRNAlist[[i]]@meta.data[["Stage"]] <- rep(Stage[i], dim(counts)[2]) 
  scRNAlist[[i]]@meta.data[["tiss_origin"]] <- rep(tiss_origin[i], dim(counts)[2])  
  scRNAlist[[i]]@meta.data[["tiss_site"]] <- rep(tiss_site[i], dim(counts)[2])  
  
  # Rename cells by adding sample identifiers
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples[i])   

  # Calculate the percentage of mitochondrial genes (MT) and ribosomal genes (RP)    
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")

  
  # Calculate the percentage of hemoglobin-related genes (HB)  
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
  scRNAlist[[i]][["percent.HB"]] <- PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 

  # Quality control 
  scRNAlist[[i]]  <- subset(scRNAlist[[i]],
                            subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & 
                            nCount_RNA > 500 & percent.mt < 20 & percent.HB < 10 )  
  
} 
scRNAlist_1 =  scRNAlist 

# CO24433  
assays <- dir("CO24433/")
dir <- paste0("CO24433/", assays)
samples = c(  "p018n", "p018", "p019n", "p019" ,"p023" ,"p024", "p027n" ,"p027" ,"p028n", "p029n" ,"p030n", "p030" ,"p031n",
              "p031",  "p032n" ,"p032", "p033n", "p033",   "p034n", "p034") 
Stage <- c(  "Normal", "IAC", "Normal",  "IIIA", "IAC", "IAC",  "Normal", "IB", "Normal","Normal","Normal","IIIA", "Normal",
             "IB",    "Normal", "IB",  "Normal", "IIB",   "Normal", "IIA")
tiss_origin<- c("Normal_origin","Tumor_origin","Normal_origin","Tumor_origin","Tumor_origin","Tumor_origin","Normal_origin",
                "Tumor_origin","Normal_origin","Normal_origin","Normal_origin","Tumor_origin","Normal_origin","Tumor_origin",
                "Normal_origin",  "Tumor_origin","Normal_origin", "Tumor_origin","Normal_origin", "Tumor_origin")
tiss_site <-  c(  "Normal", "Primary", "Normal",  "Primary",  "Primary",   "Primary",  "Normal",   "Primary",   "Normal","Normal","Normal","Primary", "Normal",
                  "Primary",    "Normal", "Primary",  "Normal", "Primary",   "Normal",  "Primary")

# Initialize an empty list to store Seurat objects
scRNAlist <- list()

# Loop through each sample, process the data and create Seurat objects
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  
  scRNAlist[[i]] <- CreateSeuratObject(counts,project = "CO24433" , min.cells=3 ) 
  # Add metadata (sample info)
  scRNAlist[[i]]@meta.data[["samples"]] <- rep( samples[i] ,dim(counts)[2]) 
  scRNAlist[[i]]@meta.data[["Stage"]] <- rep( Stage[i] ,dim(counts)[2]) 
  scRNAlist[[i]]@meta.data[["tiss_origin"]] <- rep( tiss_origin[i] ,dim(counts)[2])  
  scRNAlist[[i]]@meta.data[["tiss_site"]] <- rep( tiss_site[i] ,dim(counts)[2])  
  
  # Rename cells by adding sample identifiers 
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples[i])   
 
  # Calculate the percentage of mitochondrial genes (MT) and ribosomal genes (RP)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  
 
  # Calculate the percentage of hemoglobin-related genes (HB)  
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
  scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 

  scRNAlist[[i]]  <- subset(scRNAlist[[i]] ,subset = nFeature_RNA > 200 &nFeature_RNA < 10000 &nCount_RNA > 500& percent.mt < 20 & percent.HB < 10 )  
  
} 
scRNAlist_2  <- scRNAlist

scRNAlist <- c(scRNAlist_1,scRNAlist_2)
rm(scRNAlist_1,scRNAlist_2)

# Merge all Seurat objects into a single Seurat object for further analysis
saple_obj <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])  
saveRDS(saple_obj,file = "../../results_1/results_data/saple_obj.rds")


# Normalize the merged Seurat object
saple_obj <- readRDS("../../results_1/results_data/saple_obj.rds")
saple_obj <- NormalizeData(saple_obj) 

# Perform cell cycle scoring based on predefined gene sets
g2m_genes = cc.genes$g2m.genes %>%  CaseMatch(match = rownames(saple_obj))
s_genes   = cc.genes$s.genes  %>%  CaseMatch(match = rownames(saple_obj))
saple_obj <- CellCycleScoring(object=saple_obj, g2m.features=g2m_genes, s.features=s_genes)

# Identify variable features, scale the data and run PCA for dimensionality reduction
saple_obj <- FindVariableFeatures(saple_obj) %>% ScaleData() 
saple_obj <- RunPCA(saple_obj, features = VariableFeatures(object = saple_obj))

# Find neighbors and clusters, then perform UMAP for visualization
ElbowPlot(saple_obj, ndims = 50)
saple_obj <- FindNeighbors(saple_obj, dims = 1:25) %>% 
  FindClusters(resolution = 0.5)
saple_obj <- RunUMAP(saple_obj,dims=1:25) 

# Save the Seurat object after preprocessing (before harmony integration)
saveRDS(saple_obj,file = "../../results_1/results_data/saple_obj_har_before.rds")

# Save UMAP plots 
pdf("../../results_1/results_figures/samples_har_befor.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "samples" ,cols =  colp)
dev.off()

pdf("../../results_1/results_figures/Phase_har_befor.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "Phase" )
dev.off()

pdf("../../results_1/results_figures/orig.ident_har_befor.pdf" )
DimPlot(saple_obj, reduction = "umap",group.by = "orig.ident" )+scale_color_jco()+scale_fill_jco()
dev.off()
