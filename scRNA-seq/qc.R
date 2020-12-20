library(DropletUtils)
library(scater)

second_analysis <- function(mtx){
  rownames(mtx) <- rowData(mtx)$Symbol
  colnames(mtx) <- mtx$Barcode
  mtx_stat <- perCellQCMetrics(mtx,subsets=list(Mito=grep("mt-", rowData(mtx)$Symbol)))
  mtx_filter <- quickPerCellQC(mtx_stat,percent_subsets="subsets_Mito_percent")
  mtx <- mtx[,!mtx_filter$discard]
  tmtx <- t(counts(mtx))
  return(tmtx)
}

write_data <- function(mtx,out){
  write.csv(as.matrix(mtx),out)
}

# Read the 10x genomics data
A1A2_1 <- read10xCounts("../metadata/scRNA_seq_2020_11_1/A1A2/filtered_feature_bc_matrix")
B1A3_1 <- read10xCounts("../metadata/scRNA_seq_2020_11_1/B1A3/filtered_feature_bc_matrix")
C1C2_1 <- read10xCounts("../metadata/scRNA_seq_2020_11_1/C1C2/filtered_feature_bc_matrix")
A1A2_2 <- read10xCounts("../metadata/scRNA_seq_2020_11_7/A1A2/filtered_feature_bc_matrix")
B1A3_2 <- read10xCounts("../metadata/scRNA_seq_2020_11_7/B1A3/filtered_feature_bc_matrix")
C1C2_2 <- read10xCounts("../metadata/scRNA_seq_2020_11_7/C1C2/filtered_feature_bc_matrix")

xA1A2_1 <- second_analysis(A1A2_1)
xA1A2_2 <- second_analysis(A1A2_2)
xB1A3_1 <- second_analysis(B1A3_1)
xB1A3_2 <- second_analysis(B1A3_2)
xC1C2_1 <- second_analysis(C1C2_1)
xC1C2_2 <- second_analysis(C1C2_2)

write_data(xA1A2_1,"../result/scRNA_seq_2020_11_1-7/A1A2_1.csv")
write_data(xB1A3_1,"../result/scRNA_seq_2020_11_1-7/B1A3_1.csv")
write_data(xC1C2_1,"../result/scRNA_seq_2020_11_1-7/C1C2_1.csv")

write_data(xA1A2_2,"../result/scRNA_seq_2020_11_1-7/A1A2_2.csv")
write_data(xB1A3_2,"../result/scRNA_seq_2020_11_1-7/B1A3_2.csv")
write_data(xC1C2_2,"../result/scRNA_seq_2020_11_1-7/C1C2_2.csv")
