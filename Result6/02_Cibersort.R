memory.limit(size = 900000000)
library(Seurat)
library(CIBERSORT)  
library(reshape2)
setwd("E:/LUAD_CTS/Result6/data")
clinical <- read.table("clinical.txt", header = T,sep = "\t" )

LM22.file <- "tumor_subcluster_maker_exp_tmp.txt"
TCGA_exp.file <- ".mRNA_exp.txt"
TCGA_TME.results <- cibersort(LM22.file ,TCGA_exp.file, perm = 1000, QN = F) 
 
TME_data <- as.data.frame(TCGA_TME.results[,1:4])  
TME_data$stage <- clinical$stage  


library(tidyverse)
library(rstatix)
library(ggpubr)   
testdata <- TME_data %>%
  gather(key = "tumor_cluster", value = "score", Cellcycle,Immune, EMT, NEF) %>%
  convert_as_factor(stage, tumor_cluster)
testdata$tumor_cluster <-  factor(testdata$tumor_cluster , levels = c("NEF","Immune","Cellcycle","EMT" ))
testdata$stage <-  factor(testdata$stage , levels = c("normal", "early","advance"))  
stat.test <- testdata %>%
  group_by(tumor_cluster) %>%
  pairwise_t_test(
    score ~ stage, paired = F, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-p) # Remove details
stat.test

# 绘制图形
bxp <- ggboxplot(
  testdata, x = "tumor_cluster", y = "score", 
  fill  = "stage", palette = "jco",  outlier.shape = NA
)+ scale_fill_manual(values=c("normal"="#377EB8",   
                               "early"="#fce38a", 
                               "advance"="#E41A1C"))
# 添加显著性标记
stat.test <- stat.test %>% add_xy_position(x = "tumor_cluster")

pdf(file = "E:/LUAD_CTS/Result6/cybersort_cellratio.pdf", width = 10, height = 6)
  bxp + stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    step.increase = 0.001
  )  + theme(legend.position="right",axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

 