# Load libraries 
library(future)   
library(dplyr) 
library(infercnv)
library(copykat) 
library(ggVennDiagram)

# Increase memory limits 
options(future.globals.maxSize = 50 * 1024^3)  
memory.limit(size = 900000000 ) 

# Load the Seurat object
options(stringsAsFactors = F) 
saple_obj_anno_har_after  <- readRDS("../../results_1/results_data/saple_obj_anno_har_after.rds")
saple_EP <- readRDS(file = "../results_data/saple_EP_anno.rds")

## InferCNV input data
# Add counts
saple_EP_N <- subset(saple_EP, T_N_cell == "Normal")
saple_EP_N$seurat_clusters <- as.character(saple_EP_N$seurat_clusters)
saple_EP_N$seurat_clusters[which(saple_EP_N$tiss_origin=="Normal_origin")] = paste(saple_EP_N$seurat_clusters[which(saple_EP_N$tiss_origin=="Normal_origin")],"_","N")
saple_EP_T <- subset(saple_EP, T_N_cell == "Tumor")
saple_EP_T$seurat_clusters <- as.character(saple_EP_T$seurat_clusters)
dat_N <-  as.matrix(saple_EP_N@assays[["RNA"]]@counts)    
dat_T <-  as.matrix(saple_EP_T@assays[["RNA"]]@counts)   
dat=cbind(dat_N,dat_T) 
dat <- dat[which(rowSums(dat) > 0),]
# Add group info
groupinfo=data.frame(cellname= c(colnames(dat))  ,
                     cluster= c(saple_EP_N$seurat_clusters, saple_EP_T$seurat_clusters) ) #seurat_clusters
table(groupinfo$cluster)

# Add gene info
geneInfor  <- read.table("../results_data/infercnv/gencode_v32_gene_pos_gene_name.txt"  )
colnames(geneInfor) <- c("gene_name","chr","start","end")
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
geneInfor <- geneInfor[geneInfor$chr %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                            "chr20","chr21","chr22"),]
geneInfor$chr <-as.numeric(sub("^...","",geneInfor$chr)) 
geneInfor <- geneInfor[match(rownames(dat),geneInfor$gene_name),]
geneInfor = na.omit(geneInfor)
geneInfor <-  geneInfor[with(geneInfor, order(chr, start)),] 

dat=dat[match(geneInfor$gene_name, rownames(dat) ),groupinfo$cellname]
dim(dat)

setwd("../results_data/infercnv/") 
expFile='expFile.txt';groupFiles='groupFiles.txt';geneFile='geneFile.txt' 
write.table(dat,file = expFile,sep = '\t',quote = F)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
 
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("2 _ N","4 _ N","8 _ N","10 _ N",
                                                      "13 _ N","22 _ N"))   
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                              cluster_by_groups=TRUE,   # cluster
                              hclust_method="ward.D2",  
                              plot_steps=F,
                              HMM = F, 
                              denoise = T,
                              scale_data = T)
plot_cnv(infercnv_obj,out_dir = "F:/LUAD/results_1_identify_tumor/results_figures/infercnv/",output_format = "pdf")
saveRDS(infercnv_obj2,file = "../infercnv/infercnv_obj.rds")

##  copykat 
setwd("../copycat/")  
exp.File <- saple_EP@assays[["RNA"]]@counts
res <- copykat(rawmat=exp.File,   
               id.type="S",  
               ngene.chr=5, 
               win.size=25, 
               KS.cut=0.1, 
               sam.name="test",   
               distance="euclidean",  
               norm.cell.names="",  
               n.cores=1,
               output.seg="FLASE")
Sys.time() 


# Venn
copykat <- read.table(file = "E:/R/lung_adenocarcinoma/2个数据集/copycat/test_copykat_prediction.txt",header = T)
  
mydata <- list("two criterion" = colnames(saple_EP_T),
                 "copycat" =  copykat$cell.names[which(copykat$copykat.pred=="aneuploid")]    )
pdf(file = "../results_figures/infercnv_cat.pdf")
ggVennDiagram(mydata)
dev.off()
