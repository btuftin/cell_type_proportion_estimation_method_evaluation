library(data.table)
library(tidyverse)
library(Biobase)
library(Seurat)
library(SingleCellExperiment)
library(AnnotationDbi)

library(org.Hs.eg.db)


setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/")


# Conversion of pseudobulk files

pseudoName <- "pseudoD"

pseudobulk <- readRDS(paste0(pseudoName, "_samples.rds"))

genes <- unlist(pseudobulk$RNA@Dimnames[1])

columns(org.Hs.eg.db)
uniKeys <- keys(org.Hs.eg.db, keytype="UNIPROT")
cols <- c("SYMBOL", "ENSEMBL")
mapping <- select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT")

# Cibersort signature matrix has 547 genes
LM22 <- read_tsv("../Cibersort/LM22.txt")

LM22 <- inner_join(mapping, LM22, by = c("SYMBOL" = "Gene symbol"))

genes <- genes[genes %in% LM22$ENSEMBL]

LM22 <- LM22 %>% filter(ENSEMBL %in% genes)

LM22 <- LM22 %>% dplyr::select(c(SYMBOL, ENSEMBL)) %>% distinct()

length(unique(LM22$SYMBOL))
LM22$SYMBOL[duplicated(LM22$SYMBOL)]
# 496 out of 497. One is duplicated. Dropping the one that is named HLA-DQA2
# in ensemble.org
LM22 <- LM22 %>% filter(ENSEMBL != "ENSG00000237541")

pseudobulk$RNA <- pseudobulk$RNA[LM22$ENSEMBL,]

pb_df <- as_tibble(pseudobulk$RNA) %>% mutate(GeneSymbol = LM22$SYMBOL, .before = 1)
write_tsv(pb_df, paste0("../Cibersort/", pseudoName, "_samples_for_LM22cibersortX.txt"))


# conversion of JHS bulk
bulkRNA <- fread("/proj/yunligrp/users/jwen/NWGC_RNASeq/NWGC_RNAseq_02112022/aggregate_data/reiner_grc_rnaseq_1_gene_reads.gct.gz",
                 skip = 2, header = TRUE)

n_bulk <- 400

  samples <- colnames(bulkRNA)[2:ncol(bulkRNA)]
  set.seed(7702)
  samples <- sample(samples, n_bulk)
  bulkRNA <- bulkRNA %>% dplyr::select(Name, all_of(samples))

colnames(bulkRNA) <- paste0("g", colnames(bulkRNA))

# bulk data gene ids have version numbers, and without version numbers
# they are not unique. Removing genes with 0 expression, removing transcript
# id and keeping transcript with highest version number for duplicate ids

bulkRNA <- bulkRNA %>% mutate(gene_id_clean = str_split_i(gName,"\\.",1),
                              .before = gName)

bulk_matrix <- as.matrix(bulkRNA %>% dplyr::select(-c(gene_id_clean, gName, gDescription)))

# filtering out null genes
bulk_rowsums <- rowSums(bulk_matrix)
not_null <- bulk_rowsums > 0

bulkRNA <- bulkRNA[not_null,]

bulkRNA <- bulkRNA %>% arrange(gName) %>% group_by(gene_id_clean) %>%
  slice_head()

colnames(bulkRNA)[1] <- "Name"
bulkRNA <- bulkRNA %>% dplyr::select(-gName)

bulkRNA <- LM22 %>% inner_join(bulkRNA, by = c("ENSEMBL" = "Name"))
length(unique(bulkRNA$SYMBOL))

colnames(bulkRNA)[1] <- "Name"
bulkRNA <- bulkRNA %>% dplyr::select(-ENSEMBL)
bulkRNA <- bulkRNA %>% dplyr::select(-gDescription)

write_tsv(bulkRNA, paste0("../Cibersort/", n_bulk, "_JHS_samples_for_LM22cibersortX.txt"))
