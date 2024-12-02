if(1){
  library(future)  
  library(ggplot2)
  library(igraph)
  library(Seurat)
  library("Scissor")
  memory.limit(size = 900000000 )
  options(future.globals.maxSize = 50 * 1024^3)  
}   
 
# scRNA data 
saple_EP_T <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_T_subcluster.Rds") 
sc_dataset <- Seurat_preprocessing(saple_EP_T@assays$RNA@counts, verbose = F)
sc_dataset@meta.data = saple_EP_T@meta.data

pdf(file = "E:/LUAD_CTS/Result6/Figure7A_1.pdf")
DimPlot(sc_dataset, reduction = 'umap', label = T, repel = T,label.size = 5,group.by = "cell_state",
        cols = c("Cellcycle" = "#FA7F6F", "Immune" = "#8ECFC9", 
                 "EMT"= "#FFBE7A","NEF" = "#82B0D2"))
dev.off()

# bulk : expr and phenotype
load("E:/LUAD_CTS/data/LUAD_bulk_sur_8.rdata")

phenotype <- clin_list[[1]][,7:8]
colnames(phenotype) <- c("time", "status")
head(phenotype) 


# scissor
infos1 <- Scissor(expr_list[[1]], sc_dataset, phenotype, alpha = 0.05, 
                  family = "cox", Save_file = 'E:/LUAD_CTS/Result6/data/Scissor_LUAD_survival.RData') 
Scissor_select <- rep("Background cells", ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor + cells"
Scissor_select[infos1$Scissor_neg] <- "Scissor - cells"
table(Scissor_select)
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "Scissor")

pdf(file = "E:/LUAD_CTS/Result6/Figure7A_2.pdf", width = 9)
DimPlot(sc_dataset, reduction = 'umap', group.by = 'Scissor', 
        cols = c("Scissor - cells" = 'royalblue' ,"Scissor + cells" = 'indianred1',
               "Background cells" ='grey'), pt.size = 1.2,  
        order = c("Scissor - cells","Scissor + cells"))
dev.off()

# cell ratio
Cellratio <-as.data.frame( prop.table(table(sc_dataset$cell_state,sc_dataset$Scissor), margin = 2)) 
Cellratio$Var1 <- factor(Cellratio$Var1,levels = c('NEF', 'Immune','Cellcycle','EMT'))
pdf("E:/LUAD_CTS/Result6/Figure7B.pdf.pdf" )
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x=' ',y = 'Ratio')+
  scale_fill_manual (values=c("Cellcycle" = "#FA7F6F", "Immune" = "#8ECFC9", 
                              "EMT"= "#FFBE7A","NEF" = "#82B0D2")) +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off() 
