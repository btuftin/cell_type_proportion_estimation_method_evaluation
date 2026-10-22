library(Biobase)
library(BisqueRNA)
library(data.table)
library(tidyverse)
library(Seurat)

log_this <- function(x){
  cat(paste0(x,"\n"))
}

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/JHS_bulk_data")

reference_file_name <- "refD.rds"
n_bulk <- 100

reference_name <- commandArgs(trailingOnly = TRUE)[1]
reference_file_name <- paste0(reference_name, ".rds")

n_bulk <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
print(paste0(n_bulk, " bulk samples"))

bulkRNA <- fread(paste0("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/",
                        n_bulk, "_JHS_samples_from_reiner.txt"),
                 header = TRUE)

# This is 17GB all on its own, so size your virtual machine accordingly
# Running with a substantial fraction of the genes available requires
# tons more memory anyway
scRNA <- readRDS(file.path("../datasets_for_comparison", reference_file_name))

# We're extracting the data from the Seurat object later, and need
# to save the donor_ids first
scRNA.individual_ids <- scRNA$donor_id

scRNA.celltypes <- scRNA$cell_type


  
# Getting the counts, this frame also has row and column names with
# genes and bar codes
scRNA.data <- scRNA@assays$RNA@data

gc(verbose=FALSE)
rm(scRNA)  

bulk_matrix <- as.matrix(bulkRNA %>% select(-c(Name)))

# filtering out null genes
bulk_rowsums <- rowSums(bulk_matrix)
not_null <- bulk_rowsums > 0

bulkRNA <- bulkRNA[not_null,]

bulk_matrix <- as.matrix(bulkRNA %>% ungroup %>% select(-c(Name)))

rownames(bulk_matrix) <- bulkRNA$Name

ref_rowsums <- rowSums(as.matrix(scRNA.data))
ref_rowsums <- ref_rowsums[ref_rowsums > 0]

gc()

common_ENSGID <- intersect(names(ref_rowsums),bulkRNA$Name)
rm(bulkRNA)

log_this(paste(length(common_ENSGID), "genes in common"))

log_this("Filtering scRNA")
#seurat
scRNA.data <- scRNA.data[which(rownames(scRNA.data) %in% common_ENSGID),]
log_this("scRNA basic filter complete")
log_this(paste(nrow(scRNA.data), "rows remaining"))

scRNA.genes <- rownames(scRNA.data)
scRNA.genes <- unique(scRNA.genes)

bulk_matrix <- bulk_matrix[rownames(bulk_matrix) %in% common_ENSGID,]
bulkRNA.genes <- rownames(bulk_matrix)


  log_this("scRNA basic filter complete")
  log_this(paste(nrow(scRNA.data), "rows remaining"))

  bulk.eset <- Biobase::ExpressionSet(assayData = bulk_matrix)
  
  scRNA.sample_ids <- colnames(scRNA.data)
  
  sc.pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors = F,
                         row.names = scRNA.sample_ids,
                         SubjectName = scRNA.individual_ids,
                         cellType=scRNA.celltypes)
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName","cellType"))
  sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno, varMetadata=sc.meta)
  
  sc.eset<-Biobase::ExpressionSet(assayData = as.matrix(scRNA.data),phenoData=sc.pdata)
  
  log_this("analyzing")
  res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers= NULL,
                                                use.overlap=FALSE)  
  
  log_this("saving results")
  #  save(Est.prop, file = file.path(current_wd, paste0(test_prefix, "JHS_RNAseq_MuSiC_new_reference.RData")))
  write.table(res$bulk.props,
              file = paste0("results/Bisque_results_", reference_name, "_",
                            n_bulk, "_JHS_reiner.tsv"),
              quote = F,sep="\t",row.names = T,col.names = T)
  
