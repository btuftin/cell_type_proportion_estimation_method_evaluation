library(tidyverse)
library(Biobase)
library(Seurat)
library(SingleCellExperiment)

# Loading reference. 4GB on disk, 17GB object in R, requires 
# at minimum 20GB memory for RStudio
scRNA <- readRDS("/proj/yunligrp/users/jwen/snRNA_collection/Yazar_2022/local.rds")

# exploring some questions about the reference
# scRNA <- read_csv("/proj/yunligrp/users/chenweit/data/Yazar2022cts_by_id.csv")
# 
# scRNA <- read_csv("/proj/yunligrp/users/chenweit/data/AIDA2023/cts_by_id.csv")

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/")

# The Yazar data has cellxgene cell types, while we want the types from the
# original Yazar paper. We've received that data, and will replace the cell types
# rather than use the new data files, to not have to redo other code since
# the new data isn't in Seurat format

# The overlap isn't complete, but the new file either doesn't include the counts
# is they're malformed somehow, so remapping the cell types in the other file is
# what we'll do.

celltype_mapping <- read_csv("../../div/origin_cell_type_mapping.csv")

# The original paper dropped platelets and erythrocytes, so we will too

celltype_mapping <- celltype_mapping %>% filter(!types %in% c("Platelets", "Erythrocytes"))

cell_overlap <- colnames(scRNA@assays$RNA@data)

cell_overlap <- cell_overlap[cell_overlap %in% celltype_mapping$fingerprint]

cell_order <- tibble(fingerprint = cell_overlap)

celltype_mapping <- cell_order %>% inner_join(celltype_mapping %>% select(-c(`...1`)))

sc_Ref <- subset(scRNA, cells = cell_overlap)

# pulling old cell types to compare as sanity check

# celltype_mapping <- celltype_mapping %>% mutate(old_celltypes = sc_Ref@meta.data$cell_type)
# 
# comparison <- celltype_mapping %>% group_by(types, old_celltypes) %>%
#   summarize(n = n())

# The mappings appear sane, even if there is substantial difference in
# labeling. I don't know the cell types well enough to tell how surprising
# some of the mismatches are.

# relabel cells
rm(scRNA)

sc_Ref$cell_type <- factor(str_replace(celltype_mapping$types, " ", "_"))

# Preserve original file to loop through for
# pseudobulk creation

scRNA <- sc_Ref

# Subset samples
# set seed to enable exact replication of random sampling!
set.seed(7702)
#sample_subjects <- 1:201 # test line
sample_subjects <- unique(scRNA@meta.data[["donor_id"]])

refA <- sample(sample_subjects, 200)
refB <- sample(refA, 100)
refC <- sample(refB, 50)
refD <- sample(refC, 25)
refE <- sample(refD, 10)

for(ref in c("refA", "refB", "refC", "refD", "refE")) {
#for(ref in c("refD", "refE")) {
    
  sc_Ref = subset(sc_Ref,donor_id %in% get(ref))
  saveRDS(sc_Ref, file = paste0(ref, ".rds"))
  if(ref %in% c("refD", "refE")) {
    #save cibersort matrix
    ref_rowsums <- rowSums(sc_Ref@assays$RNA@data)
    ref_rowsums <- ref_rowsums[ref_rowsums > 0]
    
    sc_Ref <- subset(sc_Ref, features = names(ref_rowsums))
    
    sc_df <- as.matrix(sc_Ref@assays$RNA@data)
    gene_names <- sc_Ref@assays$RNA@data@Dimnames[1]
    colnames(sc_df) <- sc_Ref@meta.data$cell_type 
    sc_df <- as_tibble(sc_df) %>% mutate(GeneSymbol = gene_names, .before = 1)
    # add names after conversion, so tibble doesn't rename duplicate column names
    colnames(sc_df) <- c("GeneSymbol", as.character(scRNA@meta.data$author_cell_type))
    write_tsv(sc_df, paste0(ref, "_for_cibersort.txt"))
  }
}

rm(sc_Ref)

# create pseudo_bulk datasets from subset of remaining subjects
set.seed(7702)
pseudo_samples <- sample_subjects[!(sample_subjects %in% refA)]

pseudoA_samples <- sample(pseudo_samples, 400)
pseudoB_samples <- sample(pseudoA_samples, 200)
pseudoC_samples <- sample(pseudoB_samples, 100)
pseudoD_samples <- sample(pseudoC_samples, 50)
pseudoE_samples <- sample(pseudoD_samples, 10)

# We end up doing this twice for the biggest sample set,
# but that's still some saving in time
scRNA = subset(scRNA,donor_id %in% pseudoA_samples)
scRNA@assays$RNA@counts <- scRNA@assays$RNA@data

for(ref in c("pseudoA_samples", "pseudoB_samples", "pseudoC_samples", 
             "pseudoD_samples", "pseudoE_samples")) {
  scRNA = subset(scRNA,donor_id %in% get(ref))
  pseudo_bulk <- AggregateExpression(scRNA, group.by = "donor_id")
  cell_type_truth <- tibble(donor = scRNA$donor_id, cell_type = scRNA$cell_type)
  saveRDS(pseudo_bulk, file = paste0(ref, ".rds"))
  write_csv(cell_type_truth, paste0("cell_type_truth_", ref, ".csv"))
}
