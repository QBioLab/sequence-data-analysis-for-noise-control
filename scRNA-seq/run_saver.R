library(SAVER)

transform_by_saver <- function(mtx){
  mtx <- t(mtx)
  mtx.saver <- saver(mtx, ncores = 20,estimates.only = TRUE,size.factor = 1)
  return(mtx.saver)
}

write_data <- function(mtx,out){
  write.csv(t(mtx),out,row.names = FALSE)
}

run_main <- function(prefix){
  sample_hela_select <- read.csv(paste0('../result/scRNA_seq_2020_11_1-7/',prefix,'_hela_select.csv'))
  sample_hela_select_saver <- transform_by_saver(sample_hela_select)
  write_data(sample_hela_select_saver,paste0('../result/scRNA_seq_2020_11_1-7/',prefix,'_hela_select_saver.csv'))

  sample_f9_select <- read.csv(paste0('../result/scRNA_seq_2020_11_1-7/',prefix,'_f9_select.csv'))
  sample_f9_select_saver <- transform_by_saver(sample_f9_select)
  write_data(sample_f9_select_saver,paste0('../result/scRNA_seq_2020_11_1-7/',prefix,'_f9_select_saver.csv'))
}

# there is a side effect for saver in R, which is the '-' in gene symbol will be replaced with '.'
names = c('A1A2','B1A3','C1C2')
lapply(names,run_main)