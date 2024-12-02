library(dplyr)
library(ggplot2)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(GSVA)

# Data preparation for cNMF
saple_EP_T <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_T.Rds") 
for (si in as.character(unique(saple_EP_T$samples))) {
  metadata = saple_EP_T@meta.data %>% filter(samples == si)
  count = as.data.frame(saple_EP_T[["RNA"]]@counts[,rownames(metadata)])
  count = count[rowSums(count) >0,]
  #count = count[!str_detect(rownames(count),"^MT-"),]
  count = data.frame(t(count))
  write.table(count, file = paste(si,"count.txt", sep = "_"),quote = F,sep = "\t",row.names = T,col.names = T)
}

# cNMF Results Loading
# progroam score 
setwd("E:/LUAD_CTS/Result2/data/cnmf_usage")
assays <- dir(".")
sample = c("GSM77",'GSM78','GSM82' , 'GSM85' ,'p018','p019','p023', 
           'p024','p027', 'p030','p031', 'p032','p033' ,'p034') 

usage_1 <- list()
mean_usage <- c()
for (i in 1:length(assays)) { 
  usage_1[[i]] <- read.csv(file = assays[i] ,header = TRUE,sep=",",row.names = 1)
  names(usage_1[[i]]) <- paste(sample[i],sep = "_",gsub('[X]', '',names(usage_1[[i]])))
  usage_1[[i]] <- sweep(usage_1[[i]],1,Matrix::rowSums(usage_1[[i]]),FUN = "/")
  me_age <- colMeans(usage_1[[i]])
  names(me_age) <- names(usage_1[[i]])
  mean_usage <- c(mean_usage,me_age) 
}
mean_usage  <-  as.data.frame(mean_usage)
mean_usage$sample <-  rownames(mean_usage)
mean_usage <- mean_usage[order(mean_usage$mean_usage),]
mean_usage$sample<-factor(mean_usage$sample, 
                  levels =c(  mean_usage$sample))
mean_usage$group = "1"

# gene spectra for the programs
mean_usage <- mean_usage[mean_usage$mean_usage>=0.01,]
mean_usage$sample = as.character(mean_usage$sample)

setwd("E:/LUAD_CTS/Result2/data/cnmf_gene_spectra/") 
assays <- dir(".")
for (i in 1:length(assays)) { 
  gene_score <- as.data.frame(t(read.table(file = assays[i] ,header = TRUE,sep=",",row.names = 1)))
  colnames(gene_score) <- paste(sample[i],sep = "_", colnames(gene_score)) 
  gene_score$gene <- rownames(gene_score)
  if (sample[i] =="GSM77") {
    all_gene_score = gene_score
  }else{
    all_gene_score = all_gene_score %>% full_join(gene_score,by = "gene")
  }
}
 
# saveRDS(all_gene_score,file = "E:/LUAD_CTS/Result2/data/all_gene_score.Rds")
all_gene_score <- readRDS(file = "E:/LUAD_CTS/Result2/data/all_gene_score.Rds")
gene_score = all_gene_score[, -6]
rownames(gene_score) = all_gene_score$gene 
gene_score = gene_score[,mean_usage$sample]

# Compute correlation
res <- rcorr(as.matrix(gene_score)) 
df <- data.frame(program = rownames(res$r),programs = 'Other')

si_program <- list(NEF = c("p018_2", "p023_4","p030_4","p033_3" ,"p019_2"),
                          Immune = c("p031_5","p033_4","p027_5","p023_5","GSM85_4", "p032_3","p018_3","p030_5","GSM77_3","p023_3"),
                          Cellcycle = c("p018_5","p027_4","p034_5","p033_5" , "p024_5", "p019_4","p030_3", "p023_2","GSM77_2" ),
                          EMT = c("p031_3","p032_4","p018_4", "GSM82_5", "p019_3","p033_1"),
                          MT =  c("GSM85_2","GSM78_1","GSM77_1","GSM77_4", "p034_2", "p027_2", "p018_1","p023_1") )
 for (i in 1:length(si_program)) {
   df$programs[df$program %in% si_program[[i]]] =  names(si_program)[i]
 } 

df <- data.frame(programs = df$programs )
df$programs <- factor(df$programs,levels = c('Cellcycle','Immune','EMT','NEF','MT','Other' ))
Anno = HeatmapAnnotation(df = df, col = list(programs = c("Cellcycle" = "#FA7F6F", "Immune" = "#8ECFC9", "EMT"= "#FFBE7A", 
                                                      "NEF" = "#82B0D2","Other" = "lightgrey", "MT" = "grey")))

pdf(file = "E:/LUAD_CTS/Result2/Figure2B.pdf",width = 10, height = 8 ) 
Heatmap(res$r, col = colorRamp2(c(0,0.2,0.4, 0.6 ), c( "white","yellow",  "red",  "#67001F")),border = "black", 
        top_annotation = Anno,show_row_dend = F,show_column_dend = F,show_row_names = FALSE ,show_column_names = FALSE )
dev.off()

# four programs
signature_program = c("p018_5","p027_4","p034_5","p033_5" , 
                      "p024_5", "p019_4","p030_3", "p023_2","GSM77_2",
                      "p031_5","p033_4","p027_5","p023_5","GSM85_4",
                      "p032_3","p018_3","p030_5","GSM77_3","p023_3",
                      "p031_3","p032_4","p018_4","GSM82_5", 
                      "p019_3","p033_1","p018_2", "p023_4","p030_4","p033_3" ,"p019_2") 
# program Hallmarker 
df <- df[match(signature_program,df$program),]
df <- data.frame(programs = df$programs )
exp <- as.matrix(gene_score[,signature_program])
row.names(df) <- colnames(exp)
df$programs <- factor(df$programs,levels = c('NEF','Immune','Cellcycle','EMT'))
 
load("E:/R/Driver_SCNAs-master/Driver_SCNAs/example/data/input/m_hallmark.RData")
gsea  <- gsva(exp, m_hallmark, method="gsva",kcdf = "Gaussian") 
 
pdf(file = "E:/LUAD_CTS/Result2/Figure2C.pdf",width = 10, height = 8  )
ann_colors = list( programs = c("Cellcycle" = "#FA7F6F", "Immune" = "#8ECFC9", "EMT"= "#FFBE7A", 
                                "NEF" = "#82B0D2"))
pheatmap(gsea,    annotation_col = df,annotation_colors = ann_colors,
         color =colorRampPalette(c("#091A8C", "white","#B30000"))(100), 
         show_colnames = T, cluster_cols = F, 
         show_rownames = T, cluster_rows = T) 
dev.off()

# signature
signature <- list()
for (i in 1:length(si_program)) {
  signature_loading = all_gene_score[,c("gene",si_program[[i]])]
  used.gene = c()
  for (pi in si_program[[i]]) {
    tmp.df = signature_loading[,c("gene",pi)]
    tmp.loading = tmp.df[,2]
    names(tmp.loading) = tmp.df[,1]
    tmp.loading = tmp.loading[!is.na(tmp.loading)]
    used.gene = append(used.gene,names(tail(sort(tmp.loading),100)))
  }
  used.gene = unique(used.gene) 
  
  signature_loading = signature_loading[signature_loading$gene %in% used.gene ,]
  rownames(signature_loading) = signature_loading$gene
  signature_loading$gene = NULL
  signature_loading$total_loading = rowSums(signature_loading)
  signature_loading$average_loading = signature_loading$total_loading/length(signature_program)
  
  signature_loading = signature_loading %>% arrange(desc(average_loading))
  signature_gene = head(rownames(signature_loading),30)
  signature[[names(si_program)[i]]] <- signature_gene
}

saveRDS(signature, file = "E:/LUAD_CTS/Result2/data/signature.Rds")




 