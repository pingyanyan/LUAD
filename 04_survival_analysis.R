## Step 0: Set up the working environment ----
setwd("E:/LUAD_CTS/Result4/data/")
library(CellChat); library(patchwork); options(stringsAsFactors = FALSE)
library(mindr); library(Seurat); library(ggalluvial)
library("survminer"); library("survival")

## Step 1: data preparation----
# survival data
load("input/TCGA_GEO7.rdata") 
# gene set
cellchat <- readRDS("cellchat.rds")
cells.level <- levels(cellchat@idents)
d1 <- subsetCommunication(cellchat, slot.name = "net", thresh = 0.05,
                          sources.use = cells.level[14:17], targets.use = cells.level[1:17]) 
d2 <- subsetCommunication(cellchat, slot.name = "net", thresh = 0.05,
                          sources.use = cells.level[1:17], targets.use = cells.level[14:17]) 
df.net <- rbind(d1, d2)
pairs <- unique(df.net$interaction_name_2)
pairs <- gsub("\\(", "", pairs)
pairs <- gsub("\\)", "", pairs)
pairs <- gsub(" ", "", pairs)

# KM
p <- list()
for(pair in pairs){
  genes <- unlist(strsplit(pair, "-"))
  lig <- genes[1]
  Rec <- genes[2]
  Rec <- unlist(strsplit(Rec, "\\+"))
  
  genes <- c(lig, Rec)
  if(!all(genes %in% rownames(expr_list[[1]]))){
    cat(genes, "Pair not found\n")
    next
  }
  clin_list[[1]]$l <- t(ifelse(expr_list[[1]][lig,] > median(t(expr_list[[1]][lig,])), "High", "Low"))
  clin_list[[1]]$r <- ifelse(colSums(expr_list[[1]][Rec,]) > median(colSums(expr_list[[1]][Rec,])), "High", "Low")
  clin_list[[1]]$l_r <- paste0(clin_list[[1]]$l,"_",clin_list[[1]]$r)
  fit <- survfit(Surv(survival_time, vital_status) ~ l_r, data = clin_list[[1]])
  p[[pair]] <- ggsurvplot(fit, data = clin_list[[1]],
                           title = pair,
             font.x=14,font.y =14,
             font.tickslab = 12,
             pval = TRUE, 
             pval.size=4.5,
             size = 1.2, 
             palette = c("#E41A1C", "#a1a1a1", "#dfdfdf", "#377EB8" ),                    
             legend = c(0.9, 0.8),
             legend.title = "group", 
             tables.theme=theme_cleantable()) + 
    xlab("Time")+ylab("Survival probability")
}


