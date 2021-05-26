
# R script for RNAseq projects 3362 and 3598
# Running GENIE3 Network construction; Fig 5
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,   quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
fastgreedy <- read.table(paste(stat, "GENIE3_WRKY_NAC-fast_greedy_membership.txt", sep=""), sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
cluster    <- read.table(paste(stat, "k_means-zero_centered-9.sorted_hclust.txt", sep=""),  sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
RhiDEG     <- read.table(paste(stat, "Rhi_specific_DEGs.txt", sep=""),                      sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)

# merge
idx <- match(cluster$ID, fastgreedy$ID)
cluster$fastgreedy <- paste("module_", fastgreedy$module[idx], sep="")

idx <- cluster$ID %in% RhiDEG$Gene[RhiDEG$Cluster == "induced"]
cluster$RhiDEG[idx] <- "RhiDEG_induced"

idx <- cluster$ID %in% RhiDEG$Gene[RhiDEG$Cluster == "suppressed"]
cluster$RhiDEG[idx] <- "RhiDEG_suppresed"

idx <- is.na(cluster$RhiDEG)
cluster$RhiDEG[idx] <- ""

# plot
clusters <- data.frame(
	names=paste("Cl_0", 1:9, sep=""),
	fill=c(c_dark_red, c_red, c_cudo_magenta, c_green, c_blue, c_dark_green, c_dark_brown, c_very_dark_green, c_cudo_skyblue),
	stringsAsFactors=F)
cluster$Cluster <- factor(cluster$Cluster, levels=clusters$names)

Rhi <- data.frame(
	names=c("RhiDEG_induced", "RhiDEG_suppresed", ""),
	fill=c(c_dark_green, c_cudo_magenta, c_grey),
	stringsAsFactors=F)
cluster$RhiDEG <- factor(cluster$RhiDEG, levels=Rhi$names)

p <- ggplot(cluster, aes(x=fastgreedy, fill=Cluster)) +
	geom_bar(position="stack") +
	scale_fill_manual(values=clusters$fill) +
	theme_RTN +
	theme(axis.text.x=element_text(angle=75, hjust=1, vjust=1))
ggsave(p, file=paste(fig, "GENIE3_vs_kmeans.pdf", sep=""), width=5, height=6, bg="transparent")

p <- ggplot(cluster, aes(x=fastgreedy, group=Cluster, colour=Cluster)) +
	geom_point(stat="count") +
	geom_line(stat="count") +
	scale_colour_manual(values=clusters$fill) +
	theme_RTN +
	theme(axis.text.x=element_text(angle=75, hjust=1, vjust=1))
ggsave(p, file=paste(fig, "GENIE3_vs_kmeans2.pdf", sep=""), width=5, height=6, bg="transparent")

p <- ggplot(cluster, aes(x=fastgreedy, fill=RhiDEG)) +
	geom_bar(position="stack") +
	scale_fill_manual(values=Rhi$fill) +
	theme_RTN +
	theme(axis.text.x=element_text(angle=75, hjust=1, vjust=1))
ggsave(p, file=paste(fig, "GENIE3_vs_RhiDEG.pdf", sep=""), width=5, height=6, bg="transparent")





