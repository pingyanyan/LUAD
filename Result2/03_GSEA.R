library(clusterProfiler)
library(myenrichplot)
library(enrichplot)
library(org.Hs.eg.db)

# Load preprocessed data
load(file = "E:/LUAD_CTS/data/H1.Rdata")
all_gene_score <- readRDS("E:/LUAD_CTS/Result2/data/all_gene_score.Rds")

gene_score = all_gene_score[, -6]
rownames(gene_score) = all_gene_score$gene 
signature_program = c("p018_5","p027_4","p034_5","p033_5" , 
                      "p024_5", "p019_4","p030_3", "p023_2","GSM77_2" )
signature_program = c("p031_5","p033_4","p027_5","p023_5","GSM85_4",
                      "p032_3","p018_3","p030_5","GSM77_3","p023_3")
signature_program = c("p031_3","p032_4","p018_4","GSM82_5", 
                      "p019_3","p033_1", "p019_1", "p024_3")
signature_program = c( "p018_2", "p023_4","p030_4","p033_3" ,"p019_2"  )
signature_loading = gene_score[,signature_program]

# GSEA for each selected signature program
sigle_program <-list()
gsea_cellcyle = paste("gsea_NET",1:10,".pdf",sep = "_")
for (i in 1:dim(signature_loading)[2]) {
  print(i)
  geneList = signature_loading[,i]
  names(geneList) = rownames(signature_loading)
  geneList = sort(geneList, decreasing = T)
  head(geneList)
  egmt <- GSEA(geneList, TERM2GENE = H1, verbose=FALSE, pvalueCutoff = 0.05)
    p =   enrichplot::gseaplot2(egmt, c( 1:3  ),  pvalue_table = T, #显示p值
                           color = ggsci::pal_lancet()(3),  title = colnames(signature_loading)[i])
  sigle_program[[colnames(signature_loading)[i]]] <- p
  pdf(file =  gsea_cellcyle[i],width = 10, height = 10)
  p 
  print(p)
  dev.off()
}   


for (i in 1:dim(signature_loading)[2]) {
  print(i)
  geneList = signature_loading[,i]
  names(geneList) = rownames(signature_loading)
  geneList = sort(geneList, decreasing = T)
  
  entrezIDs <- mget(names(geneList), org.Hs.egSYMBOL2EG, ifnotfound=NA)   
  entrezIDs <- as.character(entrezIDs)
  geneList <- geneList[-which(entrezIDs=="NA")]
  entrezIDs <- entrezIDs[-which(entrezIDs=="NA")]
  names(geneList) = entrezIDs
  gseGO <- gseGO(geneList =  geneList,OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05,  pAdjustMethod = "BH")
  p =   enrichplot::gseaplot2(gseGO, c( a = c( which(gseGO@result[["Description"]] == "respiratory gaseous exchange by respiratory system" ) 
                                               ,which(gseGO@result[["Description"]] == "surfactant homeostasis" ),
                                               which(gseGO@result[["Description"]] == "lung development" ),
                                               which(gseGO@result[["Description"]] == "respiratory system development" ),
                                               which(gseGO@result[["Description"]] == "lung alveolus development" ))   ),  pvalue_table = T, #显示p值
                              color = ggsci::pal_lancet()(3),  title = colnames(signature_loading)[i])
  sigle_program[[colnames(signature_loading)[i]]] <- p
  pdf(file =  gsea_cellcyle[i],width = 10, height = 10)
  p 
  print(p)
  dev.off()
}   

 
 
 

 