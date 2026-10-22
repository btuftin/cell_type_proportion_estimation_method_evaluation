library(data.table)
library(tidyverse)
library(Biobase)
library(Seurat)
library(SingleCellExperiment)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/")

# Converting to tab separated files for CibersortX
# Conversion of the pseudobulk files and of JHS bulk below
ref_name <- "refA"

scRNA <- readRDS(paste0(ref_name, ".rds"))

length(unique(scRNA@meta.data$donor_id))

# drop ~5k genes with expression = 0

ref_rowsums <- rowSums(scRNA@assays$RNA@data)
ref_rowsums <- ref_rowsums[ref_rowsums > 0]

scRNA <- subset(scRNA, features = names(ref_rowsums))

sc_df <- as.matrix(scRNA@assays$RNA@data)
gene_names <- unlist(scRNA@assays$RNA@data@Dimnames[1])
sc_df <- as_tibble(sc_df) %>% mutate(GeneSymbol = gene_names, .before = 1)
# add names after conversion, so tibble doesn't rename duplicate column names
colnames(sc_df) <- c("GeneSymbol", as.character(scRNA@meta.data$cell_type))
write_tsv(sc_df, paste0(ref_name, "_for_cibersortX.txt"))

rm(list = ls())

#Subset even further to investigate CIBERSORTx lack of meaningful p value
#with these larger sets
samples_cells <- tibble(donor = scRNA@meta.data$donor_id, cell = scRNA@assays$RNA@data@Dimnames[[2]])

samples_cells <- samples_cells %>% group_by(donor) %>% slice_sample(n = 3)
samples_cells <- samples_cells[1:9,]

scRNA <- subset(scRNA, cells = samples_cells$cell)

ref_rowsums <- rowSums(scRNA@assays$RNA@data)
ref_rowsums <- ref_rowsums[ref_rowsums > 0]

# Downsampling the number of genes as well to see how low we need
# to go to get bad p-values

scRNA <- subset(scRNA, features = sample(scRNA@assays$RNA@data@Dimnames[[1]], 200))

sc_df <- as.matrix(scRNA@assays$RNA@data)
gene_names <- unlist(scRNA@assays$RNA@data@Dimnames[1])
sc_df <- as_tibble(sc_df) %>% mutate(GeneSymbol = gene_names, .before = 1)
# add names after conversion, so tibble doesn't rename duplicate column names
colnames(sc_df) <- c("GeneSymbol", as.character(scRNA@meta.data$cell_type))
write_tsv(sc_df, paste0("refMINI_for_cibersortX.txt"))


# Conversion of pseudobulk files

pseudoName <- "pseudoB"

pseudobulk <- readRDS(paste0(pseudoName, "_samples.rds"))

gene_names <- unlist(pseudobulk$RNA@Dimnames[1])
pb_df <- as_tibble(pseudobulk$RNA) %>% mutate(GeneSymbol = gene_names, .before = 1)
write_tsv(pb_df, paste0(pseudoName, "_samples_for_cibersortX.txt"))


# conversion of JHS bulk
bulkRNA <- fread("/proj/yunligrp/users/jwen/NWGC_RNASeq/NWGC_RNAseq_02112022/aggregate_data/reiner_grc_rnaseq_1_gene_reads.gct.gz",
                 skip = 2, header = TRUE)

n_bulk <- 400

  samples <- colnames(bulkRNA)[2:ncol(bulkRNA)]
  set.seed(7702)
  samples <- sample(samples, n_bulk)
  bulkRNA <- bulkRNA %>% select(Name, all_of(samples))

colnames(bulkRNA) <- paste0("g", colnames(bulkRNA))

# bulk data gene ids have version numbers, and without version numbers
# they are not unique. Removing genes with 0 expression, removing transcript
# id and keeping transcript with highest version number for duplicate ids

bulkRNA <- bulkRNA %>% mutate(gene_id_clean = str_split_i(gName,"\\.",1),
                              .before = gName)

bulk_matrix <- as.matrix(bulkRNA %>% select(-c(gene_id_clean, gName, gDescription)))

# filtering out null genes
bulk_rowsums <- rowSums(bulk_matrix)
not_null <- bulk_rowsums > 0

bulkRNA <- bulkRNA[not_null,]

bulkRNA <- bulkRNA %>% arrange(gName) %>% group_by(gene_id_clean) %>%
  slice_head()

colnames(bulkRNA)[1] <- "Name"
bulkRNA <- bulkRNA %>% select(-gName)

write_tsv(bulkRNA, paste0(n_bulk, "_JHS_samples_for_cibersortX.txt"))
