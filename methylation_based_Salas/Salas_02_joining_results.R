library(tidyverse)

setwd("/proj/raffieldlab/users/btuftin/Various/methylation_cell_abundance_correlation/")

result_list <- list()

result_files <- list.files("results", ".RData")

for (i in c(1:length(result_files))) {
  load(file.path("results", result_files[i]))
  result_list[[i]] <- prop_ext$prop
}

results_joined <- NULL
for (i in c(1:length(result_files))) {
  if (is.null(results_joined)) {
    results_joined <- tibble(id = rownames(result_list[[i]]))
    results_joined <- cbind(results_joined, result_list[[i]])
  } else {
    results_new <- tibble(id = rownames(result_list[[i]]))
    results_new <- cbind(results_new, result_list[[i]])
    results_new <- results_new %>% anti_join(results_joined)
    results_joined <- rbind(results_joined, results_new)
  }
}

write_csv(results_joined, "results/estimated_cell_proportions_JHS_using_Salas.csv")
