if(1){
  library(future)  
  library(Seurat)  
  library(ggplot2) 
  library(dplyr) # A Grammar of Data Manipulation
  library(RColorBrewer) # ColorBrewer Palettes
  library(grid) # The Grid Graphics Package
  library(scales) # Scale Functions for Visualizationlibrary(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)  
  library(forcats)
  library(ggrepel)
  memory.limit(size = 900000000 )
  options(future.globals.maxSize = 50 * 1024^3)  
} 

saple_EP_T <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_T_subcluster.Rds")
celltype = c("NEF","Immune","Cellcycle","EMT" )
DEGss <- list()
for (i in 1:length(celltype)) { 
  DEGs <-  FindMarkers(object = saple_EP_T, ident.1 = celltype[i], only.pos = T, 
                       group.by = "cell_state" ,min.pct = 0.25,logfc.threshold = 0.25 )
  DEGs$gene <- rownames(DEGs)
  DEGs$tumor_state <- celltype[i]
  DEGs$group <-  "adjust Pvalue > 0.05"
  DEGs <- DEGs[order(DEGs$avg_log2FC,decreasing = T),]
  DEGs$group[which(DEGs$p_val_adj <= 0.05 ) ] <- "adjust Pvalue <= 0.05" 
  DEGs$label <- NA
  DEGs$label[1:30] <- rownames(DEGs)[1:30] 
  DEGss[[i]]  <- DEGs
}

DEGs = c()
for (i in 1:length(celltype)) {  
  DEGs = rbind(DEGs, DEGss[[i]]   )  
}

DEGs$group <- factor(DEGs$group, levels = c("adjust Pvalue <= 0.05", "adjust Pvalue > 0.05"))
DEGs$tumor_state <- factor(DEGs$tumor_state, levels = c("NEF","Immune","Cellcycle","EMT") )

df_bg <- DEGs %>%
  group_by(tumor_state) %>%
  summarize(max_log2FC = max(avg_log2FC),min_log2FC = min(avg_log2FC))

pdf(file = "E:/LUAD_CTS/Result2/Figure3A.pdf")
ggplot()+
  geom_col(data = df_bg, 
           mapping = aes(tumor_state,max_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5) +
  geom_col(data = df_bg, 
           mapping = aes(tumor_state, min_log2FC),
           fill = "grey85", width = 0.8, alpha = 0.5)+
  geom_jitter(data = DEGs,
              mapping = aes(x = tumor_state, y = avg_log2FC, color = group),
              size= 1.5,width = 0.4, alpha = 0.7)+
  geom_col(data = df_bg,
           mapping = aes(x= tumor_state, y = 0.2, fill = tumor_state),
           width = 0.8)+ 
  geom_text(data=df_bg,
            mapping = aes(x=tumor_state, y=0.125, label=tumor_state),
            size = 4, color ="black",fontface = "bold")+
  scale_color_manual(values = c("#e42313", "#0061d5" ))+
  scale_fill_manual(values = c("#82B0D2","#8ECFC9","#FA7F6F","#FFBE7A"))+
  geom_text_repel(data = DEGs,
                  mapping = aes(x = tumor_state, y = avg_log2FC, label = label),
                  max.overlaps = 10000,
                  size=3, 
                  segment.color='black',
                  show.legend=FALSE)+
  theme_classic()+
  theme_minimal() + 
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 13))+
  labs(x = "group", y = "Log2FoldChange", fill= NULL, color = NULL)+
  guides(color=guide_legend(override.aes = list(size=6,alpha=1)))
 dev.off()

 
# GO enrichment

markers <- DEGs |> group_by(tumor_state) |>
 filter(p_val_adj <= 0.05) |>
 ungroup()

gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')

markers <- full_join(markers, gid, by=c('gene' = 'SYMBOL'))

gene_up.go.BP <- enrichGO(gene = markers$ENTREZID[which(markers$tumor_state == "EMT")],
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = T)

NEF <- as.data.frame(gene_up.go.BP)
NEF<- NEF[NEF$Description %in% c("respiratory gaseous exchange by respiratory system",
                                 "surfactant homeostasis","epithelial tube morphogenesis",
                                 "epithelial tube branching involved in lung morphogenesis",
                                 "antigen processing and presentation of exogenous peptide antigen",
                                 "regulation of lipid metabolic process",
                                 "negative regulation of apoptotic signaling pathway")  ,]
NEF$tumor_state <- "NEF"
  
Immune <-  as.data.frame(gene_up.go.BP)
Immune<- Immune[Immune$Description %in% c("positive regulation of cytokine production",
                                          "immune response-regulating signaling pathway",
                                          "response to reactive oxygen species",
                                          "regulation of T cell activation",
                                          "negative regulation of immune system process",
                                          "regulation of adaptive immune response",
                                          "regulation of inflammatory response")  ,]
Immune$tumor_state <- "Immune"

Cellcycle <-  as.data.frame(gene_up.go.BP)
Cellcycle<- Cellcycle[Cellcycle$Description %in% c("chromosome segregation",
                                                   "oxidative phosphorylation",
                                                   "spindle organization","protein folding",
                                                   "cytoplasmic translation",
                                          "positive regulation of cell cycle",
                                          "ATP metabolic process")  ,]
Cellcycle$tumor_state <- "Cellcycle"

EMT <- as.data.frame(gene_up.go.BP)
EMT<- EMT[EMT$Description %in% c("epithelial cell migration",
                                 "oxidative phosphorylation",
                                 "endothelial cell differentiation",
                                 "ATP metabolic process",
                                 "actin filament organization",
                                 "cell junction assembly",
                                 "cell-matrix adhesion")  ,]
EMT$tumor_state <- "EMT"


NEF$p.adjust <- -log(NEF$p.adjust)
NEF <- NEF[order(NEF$Count),]
NEF$Description <- factor(NEF$Description, levels = NEF$Description)
colors = c("#418ebb", "#c5e7b5", "#ffc99a","#d94d57")

pdf(file = "E:/LUAD_CTS/Result2/Figure3B.pdf",width = 10, height = 5)
ggplot(data=NEF, aes(x=Count,y= Description, fill=p.adjust)) +
  geom_bar(stat="identity", width=0.8) +scale_fill_gradientn(colors  = colors ) +   theme_classic() +
  xlab("Num of Genes") + ylab("GO term") + labs(title = "NEF state")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) +
  geom_text(aes(x=0, label = Description), vjust = 0, hjust = 0, size = 3)  +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
dev.off()

