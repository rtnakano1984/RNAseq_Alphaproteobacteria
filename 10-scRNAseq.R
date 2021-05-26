#!/netsctach/dep_psl/grp_psl/ThomasN/tools/R403/bin/Rscript
# Note!! A different R binary than usual
# R script to analyze scRNAseq data from 

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(readr,  quietly=T, warn.conflicts=F)
library(Seurat,  quietly=T, warn.conflicts=F)
library(dplyr,  quietly=T, warn.conflicts=F)


ref_data       <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/ref_data/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"

df     <-   read_csv(paste(ref_data, "GSE123818_Root_single_cell_wt_datamatrix.csv.gz", sep=""))
RhiDEG <- read.table(paste(stat, "Rhi_specific_DEGs.txt", sep=""), sep="\t", header=T, row.names=NULL, stringsAsFactors=F)

mat <- as.matrix(df[, -1])
rownames(mat) <- df$X1


dat <- CreateSeuratObject(counts=mat, project="Denyer_Ma", min.cells=3, min.features=200)
norm <- NormalizeData(dat) %>% FindVariableFeatures() %>% ScaleData()

norm_mat <- as.matrix(GetAssayData(norm))


# norm_mat <- norm_mat[1:1500, 1:1500]

cor <- cor(t(norm_mat))

idx <- colSums(is.na(cor)) == nrow(cor)-1
cor <- cor[!idx, !idx]

rank <- sapply(1:nrow(cor), function(x){
	rank(-cor[,x])
	})

mr <- data.frame(ID=rownames(rank), sqrt(rank*t(rank)))
colnames(mr) <- c("ID", rownames(rank))

write_delim(mr, paste(processed_data, "GSE123818_MR.text.gz", sep=""), delim="\t")


idx <- colnames(mr) %in% c("ID", RhiDEG$Gene)
mr_rhi <- mr[, idx]
write.table(mr_rhi, paste(processed_data, "GSE123818_MR.RhiDEG.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")


