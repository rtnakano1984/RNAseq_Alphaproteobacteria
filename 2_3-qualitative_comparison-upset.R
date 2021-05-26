
# R script for RNAseq projects 3362 and 3598
# k mean clustering based on log2cpm; Fig S1
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,      quietly=T, warn.conflicts=F)
# library(dplyr,        quietly=T, warn.conflicts=F)
# library(eulerr,      quietly=T, warn.conflicts=F)
# library(reshape2,     quietly=T, warn.conflicts=F)
# library(ggplot2,      quietly=T, warn.conflicts=F)
# library(patchwork,    quietly=T, warn.conflicts=F)
# library(ComplexHeatmap,    quietly=T, warn.conflicts=F)
library(UpSetR,    quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))
# theme_RTN <- theme_RTN + theme(panel.spacing=unit(3, "pt"))

sig <- read.table(paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F) %>% as.matrix

# rename
idx <- match(colnames(sig), treatment_1$names)
colnames(sig) <- treatment_1$ID[idx]

# up and down
DE_up <- sig
idx <- DE_up != 1
DE_up[idx] <- 0
colnames(DE_up) <- paste(colnames(sig), "_up", sep="")

DE_down <- -sig
idx <- DE_down != 1
DE_down[idx] <- 0
colnames(DE_down) <- paste(colnames(sig), "_down", sep="")

DE_df <- data.frame(ID=rownames(sig), DE_up, DE_down, row.names=NULL)


#
pdf(paste(fig, "upset_diagrams.pdf", sep=""), width=15, height=5)
	upset(DE_df, order.by="freq", sets=colnames(DE_df)[-1], mb.ratio=c(.35, .65), nintersects=150, keep.order=T)
dev.off()

pdf(paste(fig, "upset_diagrams.up.pdf", sep=""), width=10, height=4)
	upset(DE_df, order.by="freq", sets=paste(treatment_1$ID[-1], "_up", sep=""), mb.ratio=c(.35, .65), nintersects=150, keep.order=T, sets.bar.color=c_dark_green)
dev.off()

pdf(paste(fig, "upset_diagrams.down.pdf", sep=""), width=10, height=4)
	upset(DE_df, order.by="freq", sets=paste(treatment_1$ID[-1], "_down", sep=""), mb.ratio=c(.35, .65), nintersects=150, keep.order=T, sets.bar.color=c_cudo_magenta)
dev.off()




# 







