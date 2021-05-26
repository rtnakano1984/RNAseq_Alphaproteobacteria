#!/netscratch/dep_psl/grp_psl/ThomasN/tools/R403/bin/Rscript

# R script for RNAseq projects 3362 and 3598
# Comparing responses between treatments; Fig 1
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(reshape2,      quietly=T, warn.conflicts=F)
library(patchwork,     quietly=T, warn.conflicts=F)
library(ggplot2,       quietly=T, warn.conflicts=F)
library(vegan,         quietly=T, warn.conflicts=F)
library(RVAideMemoire, quietly=T, warn.conflicts=F)
source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")
theme_RTN <- theme_RTN + theme(panel.spacing=unit(1, "pt"))

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design    <-            read.table(paste(original_data,  "design.txt", sep=""),            sep="\t", header=T,              check.names=F, stringsAsFactors=F)
norm      <-  as.matrix(read.table(paste(processed_data, "zero_centered_DEG.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))

# PCC computation based on zero-centered
cor <- cor(norm)
d <- as.dist(1-cor)

# permanova
set.seed(0)
design$treatment_1 <- factor(design$treatment_1, levels=treatment_1$names)
adonis_strain <- adonis2(d ~ treatment_1, design, sqrt.dist=T, by="term")
sink(paste(stat, "permanova-treatment_1.txt", sep=""))
	print(adonis_strain)
sink()


# pairwise permanova btw lineages
idx <- design$treatment_1 == treatment_1$names[1]
design$lineage[idx] <- "axenic"

idx <- design$treatment_1 %in% treatment_1$names[2:5]
design$lineage[idx] <- "Rhizobiales"

idx <- design$treatment_1 %in% treatment_1$names[6:7]
design$lineage[idx] <- "sister_lineage"

set.seed(0)
pair_permanova_strain <- pairwise.perm.manova(d, design$lineage, p.method="fdr", R2=T)
sink(paste(stat, "pairwise_PERMANOVA-lineage.txt", sep=""))
	print(pair_permanova_strain)
sink()

results <- cbind(melt(pair_permanova_strain$R2.value), melt(pair_permanova_strain$p.value))
results <- na.omit(results)[, c(2, 1, 3, 6)]
colnames(results) <- c("lineage_1", "lineage_2", "R2", "FDR")
write.table(results, file=paste(stat, "pairwise_PERMANOVA-lineage-table.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")










