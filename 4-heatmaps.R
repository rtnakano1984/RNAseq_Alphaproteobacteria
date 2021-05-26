
# R script for RNAseq projects 3362 and 3598
# creating various heatmaps; Fig S1
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(reshape2,  quietly=T, warn.conflicts=F)
library(patchwork, quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))
theme_RTN <- theme_RTN + theme(panel.spacing=unit(3, "pt"))


# data import
design  <-           read.table(paste(original_data,  "design.txt",                      sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)
norm    <- as.matrix(read.table(paste(processed_data, "zero_centered_DEG.txt", sep=""),           sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
log2cpm <- as.matrix(read.table(paste(processed_data, "log2cpm_DEG.txt",                 sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
logFC   <- as.matrix(read.table(paste(processed_data, "logFC_DEG.txt",                   sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
sig     <- as.matrix(read.table(paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
cluster <- read.table(paste(stat, "k_means-zero_centered-9.sorted_hclust.txt", sep=""),           sep="\t", header=T,              check.names=F, stringsAsFactors=F)


# saturate function

# saturate parameter
saturate_cut <- .05





# log2cpm #########################################
melt_log2cpm <- melt(log2cpm)

idx <- match(melt_log2cpm$Var2, design$Library)
melt_log2cpm <- cbind(melt_log2cpm, design[idx,])

melt_log2cpm$Var1        <- factor(melt_log2cpm$Var1,        levels=cluster$ID)
melt_log2cpm$treatment_1 <- factor(melt_log2cpm$treatment_1, levels=treatment_1$names)
levels(melt_log2cpm$treatment_1) <- treatment_1$ID

melt_log2cpm$value <- saturate(melt_log2cpm$value, saturate_cut)

idx <- match(melt_log2cpm$Var1, cluster$ID)
melt_log2cpm$Cluster <- factor(cluster$Cluster[idx], levels=unique(cluster$Cluster))
levels(melt_log2cpm$Cluster) <- 1:length(unique(melt_log2cpm$Cluster))

p_log2cpm <- ggplot(melt_log2cpm, aes(x=rep, y=Var1, fill=value)) +
	geom_tile(alpha=1) +
	scale_fill_gradient2(low=c_cudo_skyblue, mid=c_black, high=c_yellow, midpoint=0, guide=guide_colorbar(barwidth=5, barheight=.75, frame.colour="black", frame.linewidth=1)) +
	labs(x="", y="", fill="log2cpm") +
	theme_RTN +
	facet_grid(Cluster ~ treatment_1, drop=TRUE, scales="free", space="free", switch="both") +
	theme(legend.position="top",
		legend.text=element_text(size=6),
		legend.title=element_text(size=8),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_text(angle=180, size=6),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())






#Z scores #########################################
melt_norm <- melt(norm)

idx <- match(melt_norm$Var2, design$Library)
melt_norm <- data.frame(melt_norm, design[idx,], stringsAsFactors=F)

melt_norm$Var1        <- factor(melt_norm$Var1,        levels=cluster$ID)
melt_norm$treatment_1 <- factor(melt_norm$treatment_1, levels=treatment_1$names) 
levels(melt_norm$treatment_1) <- treatment_1$ID

melt_norm$value <- saturate(melt_norm$value, saturate_cut)

idx <- match(melt_norm$Var1, cluster$ID)
melt_norm$Cluster <- factor(cluster$Cluster[idx], levels=unique(cluster$Cluster))

p_norm <- ggplot(melt_norm, aes(x=rep, y=Var1, fill=value)) +
	geom_tile(alpha=1) +
	scale_fill_gradient2(low="#0068B7", mid="white", high="#F39800", midpoint=0, guide=guide_colorbar(barwidth=5, barheight=.75, frame.colour="black", frame.linewidth=1)) +
	labs(x="", y="", fill="Zero-centered\nexpression") +
	theme_RTN +
	facet_grid(Cluster ~ treatment_1, drop=TRUE, scales="free", space="free", switch="both") +
	theme(legend.position="top",
		legend.text=element_text(size=6),
		legend.title=element_text(size=8, hjust=1),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_blank(),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())







# logFC #########################################
melt_logFC <- melt(logFC)

melt_logFC$Var1 <- factor(melt_logFC$Var1, levels=cluster$ID)
melt_logFC$Var2 <- factor(melt_logFC$Var2, levels=treatment_1$names[-1])
levels(melt_logFC$Var2) <- treatment_1$ID[-1]

melt_logFC$value <- saturate(melt_logFC$value, saturate_cut)

idx <- match(melt_logFC$Var1, cluster$ID)
melt_logFC$Cluster <- factor(cluster$Cluster[idx], levels=unique(cluster$Cluster))

p_logFC <- ggplot(melt_logFC, aes(x=Var2, y=Var1, fill=value)) +
	geom_tile(alpha=1) +
	scale_fill_gradient2(low=c_cudo_magenta, mid="white", high=c_dark_green, midpoint=0, guide=guide_colorbar(barwidth=5, barheight=.75, frame.colour="black", frame.linewidth=1)) +
	labs(x="", y="", fill="logFC") +
	facet_grid(Cluster ~ Var2, drop=TRUE, scales="free", space="free", switch="both") +
	theme_RTN +
	theme(legend.position="top",
		legend.text=element_text(size=6),
		legend.title=element_text(size=8),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_blank(),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())





# significance #########################################
idx <- rowSums(sig != 0) > 0
melt_sig <- melt(sig[idx,])

melt_sig$Var1 <- factor(melt_sig$Var1, levels=cluster$ID)
melt_sig$Var2 <- factor(melt_sig$Var2, levels=treatment_1$names[-1])
levels(melt_sig$Var2) <- treatment_1$ID[-1]

idx <- melt_sig$value != 0
melt_sig$value <- as.numeric(idx)

idx <- match(melt_sig$Var1, cluster$ID)
melt_sig$Cluster <- factor(cluster$Cluster[idx], levels=unique(cluster$Cluster))

p_sig <- ggplot(melt_sig, aes(x=Var2, y=Var1, fill=factor(value))) +
	geom_tile(alpha=1) +
	scale_fill_manual(values=c("0"="white", "1"="red"), labels=c("0"="ns", "1"="significant"), guide=guide_legend(override.aes=list(colour=c_black, size=.25))) +
	labs(x="", y="", fill="") +
	facet_grid(Cluster ~ Var2, drop=TRUE, scales="free", space="free", switch="both") +
	theme_RTN +
	theme(legend.position="top",
		legend.key.size=unit(.7, "lines"),
		legend.text=element_text(size=6),
		legend.title=element_text(size=8),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_blank(),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())





# composit and export #########################################
p <- p_log2cpm + p_norm + p_logFC + p_sig + plot_layout(width=c(2, 2, 1, 1)) & theme(plot.background = element_blank())
ggsave(p, file=paste(fig, "FigS1-heatmaps.png", sep=""), width=8, height=4.5, bg="transparent")



