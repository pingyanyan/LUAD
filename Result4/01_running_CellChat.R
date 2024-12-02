#### Select four tumor states and tumor-related cells for cell-cell communication analysis


## Step 0: Set up the working environment
setwd("E:/LUAD_CTS/Result4/data/input/")
library(CellChat); library(patchwork) 
library(mindr); library(Seurat); library(ggalluvial)
options(stringsAsFactors = FALSE)


## Step 1: Preparing the input data for CellChat
# 1) Load Data
files <- list.files()
seu_list <- lapply(files, readRDS)

# 2) Select four tumor states and tumor-related cells
names(seu_list) <- gsub("\\.rds$", "", files, ignore.case = TRUE)
seu_list[[1]] <- subset(seu_list[[1]], subset = M_Cell_type %in%  c("IgA+ PC", "IgG+ PC", "Igm+ PC(stress)", "memory_B", "memory_B(stress)", "naive B"))
seu_list[[2]] <- subset(seu_list[[2]], subset = M_Cell_type == "tumor ECs")
seu_list[[3]] <- subset(seu_list[[3]], subset = M_Cell_type == "THBS2+ Myofibroblasts")
seu_list[[4]] <- subset(seu_list[[4]], subset = M_Cell_type == "mo_mac")
seu_list[[5]] <- subset(seu_list[[5]], subset = M_Cell_type %in% c("CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex", "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex"))
seu_list[[6]] <- subset(seu_list[[6]], subset = cell_state %in% c("EMT", "Immune", "Cellcycle", "NEF"))
seu_list[[6]]$M_Cell_type <- seu_list[[6]]$cell_state

# 3) Integrate Seurat objects and normalize, set idents as cell types
seu_obj <- merge(seu_list[[1]], y = seu_list[-1])
seu_obj <- NormalizeData(seu_obj)
Idents(seu_obj) <- seu_obj$M_Cell_type


## Step 2: Run CellChat
# 1) Create a CellChat object: including the normalized expression matrix and cell types
cellchat <- createCellChat(seu_obj) 

# 2) Set the ligand-receptor database
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# 3) Filter genes not in the ligand-receptor database to reduce the data size
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

# 4) Run the main CellChat functions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Visualize the number of interactions in a circular plot
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, 
                 weight.scale = TRUE,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 vertex.weight = 20,  # 从groupSize变成固定大小
                 weight.scale = TRUE,
                 label.edge= F, title.name = "Interaction weights/strength")

# Save the CellChat results
saveRDS(cellchat, "../cellchat.rds")
cellchat <- readRDS("../cellchat.rds")
