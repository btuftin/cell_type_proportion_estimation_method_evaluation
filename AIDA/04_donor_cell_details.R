library(data.table)
library(tidyverse)

setwd("/proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/AIDA/datasets_for_comparison/")

metadata_ref <- fread("reference_all/metadata.csv", 
                      select =c("donor_id", "Country", "sex", "development_stage"))

nrow(metadata_ref) # I.e. number of cells
# 394139

metadata_ref <- metadata_ref %>% unique()
nrow(metadata_ref) # I.e. number of donors
#200

metadata_test <- fread("test_set_all/metadata.csv", 
                       select =c("donor_id", "Country", "sex", "development_stage"))
nrow(metadata_test) # I.e. number of cells
#664770

metadata_test <- metadata_test %>% unique()
nrow(metadata_test) # I.e. number of donors
#317


gc()

metadata <- rbind(metadata_ref, metadata_test)

metadata <- metadata %>% mutate(age = as.numeric(str_split_i(development_stage, "-", 1)))

avg_cells <- (394139 + 664770) / 517
# 2048

mean_age <- mean(metadata$age)
# 42
range_age <- range(metadata$age)
# 19 75
pct_female <- 100*sum(metadata$sex == "female")/nrow(metadata)
#55.5%
table(metadata$Country)
