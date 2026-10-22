library(data.table)
library(tidyverse)

bulkRNA <- fread("/proj/yunligrp/users/jwen/NWGC_RNASeq/NWGC_RNAseq_02112022/aggregate_data/reiner_grc_rnaseq_1_gene_reads.gct.gz",
                 skip = 2, header = TRUE)
samples <- colnames(bulkRNA)[3:ncol(bulkRNA)]


rm(bulkRNA)

# select samples in common with flow_cytometry

flow_cyt_cells <- read_tsv("/proj/raffieldlab/projects/JHS/flow_jhs/JHS_FLOW_data_transmit_20240408_TOPMed_ID.txt")

ID_map<- fread("/proj/yunligrp/users/jwen/NWGC_RNASeq/NWGC_RNAseq_02112022/aggregate_data/rnaseq_1_gene_reads.gct_IDsconversion.map")

ID_map <- ID_map  %>% mutate(subjid = as.character(JHSID)) %>% select(c(subjid, NWDID))

flow_cyt_cells <- ID_map %>% inner_join(flow_cyt_cells) %>% select(-NWDID)

samples <- samples[samples %in% flow_cyt_cells$subjid]

