library(CytoTRACE)
library(ggpubr)
library(dplyr)
library(monocle)
library(Seurat)
memory.limit(size = 900000000 )
#tumor monocle
mono <- readRDS("E:/LUAD_CTS/Result5/data/mono.RDS")
saple_EP_T <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_T.Rds")
#monocle 
cd <- mono(saple_EP_T)
cd$cell_state  <- factor(cd$cell_state, levels = c("NEF","Immune","Cellcycle","EMT"))

pdf(file = "E:/LUAD_CTS/Result5/Figure6A_1.pdf")
plot_cell_trajectory(cd,color_by="cell_state", size=1,show_backbone=TRUE  ) +  
                    scale_color_manual(breaks = c("NEF","Immune","Cellcycle",  "EMT"), 
                                       values=c("#82B0D2","#8ECFC9","#FA7F6F",  "#FFBE7A"))  
dev.off()

pdf(file = "E:/LUAD_CTS/Result5/Figure6A_2.pdf")
plot_cell_trajectory(cd,color_by="Pseudotime", size=1,show_backbone=TRUE  ) +  scale_color_gradient(low = "grey",high = "#E41A1C")
dev.off()


data <- cd@phenoData@data
mpg <- data[which(data$State == "4"),]
pdf(file = "E:/LUAD_CTS/Result5/Figure6A_3.pdf" )
group_by(mpg, cell_state) %>%
  summarise(percent = n() / nrow(mpg)) %>%
  ggplot(aes(x = factor(1), y = percent, fill = cell_state)) +  theme_bw() +
  geom_col(colour = "white") + 
  coord_polar(theta = "y", start = 1.65)   +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+scale_fill_manual(values = c("Cellcycle" = "#FA7F6F", "Immune" = "#8ECFC9", "NEF" = "#82B0D2",
                                 "EMT"= "#FFBE7A"))
dev.off()

# BEAM 分支
signature = readRDS(file = "E:/LUAD_CTS/Result2/data/signature.Rds")
express_genes <- VariableFeatures(FindVariableFeatures(saple_EP_T,selection.method = "vst",nfeatures = 2000)) 


BEAM_res <- BEAM(cd[express_genes,], branch_point = 1, cores = 10) 
TF_name[which(TF_name %in% BEAM_res$gene_short_name == "TRUE")]
BEAM_res <- BEAM_res[order(BEAM_res$qval),] 
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf(file = "E:/LUAD_CTS/Result5/Figure6A_2.pdf" ,width = 5, height = 8)
 
p = plot_genes_branched_heatmap(cd[row.names(subset(BEAM_res, qval < 0.01 & pval<0.01  )),],
                            branch_point = 1, 
                            num_clusters = 4, 
                            cores = 10, 
                            branch_labels = c("branch 2", "branch 3"),
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T) 
                            
dev.off()





# AT2  Club origin 
pdf(file = "../results_figures/mono_AT2_Club.pdf" ,width = 6, height = 5)
plot_genes_in_pseudotime(cd[c("SFTPC","SCGB1A1"),], color_by = "cell_state") +  
  scale_color_manual(breaks = c("NEF","Immune","Cellcycle",  "EMT"), 
                     values=c("#82B0D2","#8ECFC9","#FA7F6F",  "#FFBE7A"))
dev.off()
#EMT in different state
exp <- data.frame( EMT = cd$EMT1, State  = cd$State, cell_state=cd$cell_state ) 
exp  <- exp[which(exp$cell_state == "EMT"),]
pdf(file = "../results_figures/EMT_in_state.pdf"  )
 ggboxplot(exp, 'State' ,'EMT' ,  fill = "State", outlier.shape = NA) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'right') + 
  stat_compare_means(comparisons = list(c("1","2"),c("1","3") ,c("2","3") ),
                      method = 'wilcox.test') +
  scale_fill_manual(values=c("1" = "#82B0D2", '3' ="#FFBE7A",'2' = "#FA7F6F" ))
dev.off()
#stem score cytoTrace   
setwd("../results_data/cyto_stem/")
mat_3k <- as.matrix(saple_EP_T@assays$RNA@counts) 
results <- CytoTRACE(mat = mat_3k)
cd$StemScore <- results[["CytoTRACE"]] 
saveRDS(cd,file = "../../results_data/mono_tumor.Rds")

pdf(file = "E:/LUAD_CTS/Result5/Figure6_4")
plot_cell_trajectory(cd,color_by="StemScore", size=1,show_backbone=TRUE  ) +  scale_color_gradient(low = "grey",high = "#E41A1C")  
dev.off()
