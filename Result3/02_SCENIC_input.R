library(SCENIC)
memory.limit(size = 900000000)
saple_EP_T <- readRDS("E:/LUAD_CTS/Result2/data/saple_EP_T_subcluster.Rds")  

setwd("E:/LUAD_CTS/Result3/data/SCENIC/int")  

exprMat <-as.matrix(saple_EP_T@assays[["RNA"]]@counts  )   
exprMat <- exprMat[which(rowSums(exprMat)>0),]
 

hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org="hgnc", dbDir="int", dbs = hg38_dbs, nCores=1)
 
# Gene filter/selection  (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)   
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
 
write.csv(t(as.matrix(exprMat_filtered)),file = "tumor_counts.csv")

