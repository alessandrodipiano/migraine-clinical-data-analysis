
path_longitudinal <- "/Users/alessandro.distefano/Desktop/AI/3rd year/Medical/migraine-clinical-data-analysis/data/imputed/imputed_longitudinal"
baseline <- "/Users/alessandro.distefano/Desktop/AI/3rd year/Medical/migraine-clinical-data-analysis/data/imputed/imputed_data_baseline"

library(miceadds)

folder <- path_longitudinal

files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
files <- sort(files)  # ensure imp_01, imp_02, ... order

completed_list <- lapply(files, read.csv, stringsAsFactors = FALSE)

imp_long_reloaded <- miceadds::datlist2mids(completed_list)

fit <- with(imp_long_reloaded, lm(MMDs ~ CYCLE + MONTH))
pool(fit)