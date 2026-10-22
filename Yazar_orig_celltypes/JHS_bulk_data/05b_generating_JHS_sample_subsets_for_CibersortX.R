library(data.table)
library(tidyverse)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/")

# Conversion of bulk files
# manually select one of the two options

bulk_n <- "100"
bulk_n <- "400"

bulkRNA <- fread(paste0(bulk_n, "_JHS_samples_from_reiner.tsv"))

colnames(bulkRNA) <- unlist(bulkRNA[1,])

bulkRNA <- bulkRNA[-1,]
# bulk data gene ids have version numbers, and without version numbers
# they are not unique. Removing genes with 0 expression, removing transcript
# id and keeping transcript with highest version number for duplicate ids

bulkRNA <- bulkRNA %>% mutate(gene_id_clean = str_split_i(Name,"\\.",1),
                              .before = Name)

bulk_matrix <- as.matrix(bulkRNA %>% select(-c(gene_id_clean, Name)))

# # filtering out null genes (turns out there are none in this set)
# bulk_rowsums <- rowSums(bulk_matrix)
# not_null <- bulk_rowsums > 0
# 
# bulkRNA <- bulkRNA[not_null,]

bulkRNA <- bulkRNA %>% arrange(Name) %>% group_by(gene_id_clean) %>%
  slice_head()

bulkRNA <- bulkRNA %>% select(Name)
colnames(bulkRNA)[1] <- "Name"

write_tsv(bulkRNA, paste0(n_bulk, "_JHS_samples_for_cibersortX.txt"))
