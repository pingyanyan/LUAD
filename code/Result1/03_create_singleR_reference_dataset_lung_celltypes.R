library(scater) 
# Read the metadata of the reference dataset
Ref_annot <- read.csv("../../data/single_cell/singleR_ref/krasnow_hlca_10x_metadata.csv")
Ref_annot <- Ref_annot[which(Ref_annot$tissue == "lung"), ]

# Lower the resolution of cell labels
a  <-  as.data.frame.array(table(Ref_annot$free_annotation)) ;
a$celltype  =  c("Fibroblast","Fibroblast","Epithelial","Epithelial","Fibroblast","Endothelial","B",
                 "Epithelial","Mast","Mast","Endothelial","Endothelial","Endothelial","Endothelial",  
                 "Endothelial","Endothelial" ,"T","T", "T","T","Epithelial" ,"Myeloid","Epithelial","Epithelial","Myeloid",
                 "Fibroblast","Epithelial","Myeloid" , "Myeloid","Epithelial","Fibroblast","Endothelial","Myeloid","Fibroblast",
                 "Epithelial" ,"Myeloid","Myeloid","Fibroblast","NK","T","Epithelial","Myeloid", "Myeloid","Fibroblast","B","Myeloid","Myeloid",
                 "Epithelial","Myeloid","T","Epithelial" ,"Epithelial","Epithelial","Epithelial","Myeloid","Fibroblast" ,"Endothelial")
Ref_annot$celltype = a$celltype[match(Ref_annot$free_annotation, rownames(a))]

# Read the gene expression dataRef <- read.csv("../../data/single_cell/singleR_ref/krasnow_hlca_10x_UMIs.csv")
Ref <- textshape::column_to_rownames(Ref, loc = 1)
Ref <- Ref[,match(Ref_annot$X,colnames(Ref))]

# Filter the reference data to keep genes expressed in at least 3 cells and cells with at least 200 expressed genes
Ref <- Ref[which(rowSums(Ref!=0)>3) , which(colSums(Ref!=0)>200)]
Ref_annot <- Ref_annot[match(colnames(Ref), Ref_annot$X), c("X", "free_annotation", "celltype")]

# Create a SummarizedExperiment object with the counts matrix and the sample metadata
Ref_sce <- SummarizedExperiment(assays=list(counts=Ref), colData = Ref_annot) 
Ref_sce <- logNormCounts(Ref_sce) 
 
saveRDS(Ref_sce,file = "../results_data/lung_58celltype.rds")


