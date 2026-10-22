library(data.table)
library(tidyverse)

bulkRNA <- fread("/proj/yunligrp/users/jwen/NWGC_RNASeq/NWGC_RNAseq_02112022/aggregate_data/reiner_grc_rnaseq_1_gene_reads.gct.gz",
                 skip = 2, header = TRUE)
samples <- colnames(bulkRNA)[3:ncol(bulkRNA)]

# bulk data gene ids have version numbers, and without version numbers
# they are not unique. Removing genes with 0 expression, removing transcript
# id and keeping transcript with highest version number for duplicate ids

bulkRNA <- bulkRNA %>% mutate(gene_id_clean = str_split_i(Name,"\\.",1),
                              .before = Name)

bulk_matrix <- as.matrix(bulkRNA %>% select(-c(gene_id_clean, Name, Description)))

# filtering out null genes
bulk_rowsums <- rowSums(bulk_matrix)
not_null <- bulk_rowsums > 0

bulkRNA <- bulkRNA[not_null,]

bulkRNA <- bulkRNA %>% arrange(Name) %>% group_by(gene_id_clean) %>%
  slice_head()


# select samples in common with flow_cytometry

flow_cyt_cells <- read_tsv("/proj/raffieldlab/projects/JHS/flow_jhs/JHS_FLOW_data_transmit_20240408_TOPMed_ID.txt")

ID_map<- fread("/proj/yunligrp/users/jwen/NWGC_RNASeq/NWGC_RNAseq_02112022/aggregate_data/rnaseq_1_gene_reads.gct_IDsconversion.map")

ID_map <- ID_map  %>% mutate(subjid = as.character(JHSID)) %>% select(c(subjid, NWDID))

flow_cyt_cells <- ID_map %>% inner_join(flow_cyt_cells) %>% select(-NWDID)

samples <- samples[samples %in% flow_cyt_cells$subjid]

bulkRNA <- bulkRNA %>% ungroup()

bulkRNA <- bulkRNA %>% select(-"Name")
colnames(bulkRNA)[1] <- "Name"


for (n_bulk in c(100, 400)){
  set.seed(7702)
  loop_samples <- sample(samples, n_bulk)
  loop_bulkRNA <- bulkRNA %>% select(Name, all_of(loop_samples))
  
  write_tsv(loop_bulkRNA, paste0("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/",
                            n_bulk, "_JHS_samples_from_reiner.txt"))
}
