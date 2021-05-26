
# R script for RNAseq projects 3362 and 3598
# GO enrichment analysis of k means clusters (based on logFC); Fig 2
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(ggplot2,         quietly=T, warn.conflicts=F)
library(dplyr,           quietly=T, warn.conflicts=F)
library(stringr,         quietly=T, warn.conflicts=F)
library(reshape2,        quietly=T, warn.conflicts=F)
library(patchwork,        quietly=T, warn.conflicts=F)

# diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
temp           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/temp/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# load data
logFC   <- read.table(paste(processed_data, "logFC_DEG.txt", sep=""),                                  sep="\t", header=T, check.names=F, stringsAsFactors=F, comment.char="", quote="", row.names=1) %>% as.matrix
design  <- read.table(paste(original_data,  "design.txt", sep=""),                                     sep="\t", header=T, check.names=F, stringsAsFactors=F, comment.char="", quote="")
cluster <- read.table(paste(stat, "k_means-zero_centered-9.sorted_hclust.txt", sep=""),                sep="\t", header=T, check.names=F, stringsAsFactors=F)


# =============== for plotting logFC =============== #
# melt
logFC_melt <- melt(logFC)

# merge cluster info
idx <- match(logFC_melt$Var1, cluster$ID)
logFC_melt$Cluster <- cluster$Cluster[idx]

# cluster-wise mean logFC
logFC_mean <- logFC_melt %>% group_by(Cluster, Var2) %>% summarise(mean=mean(value))

# merge meta data
logFC_mean$fill    <- saturate(logFC_mean$mean, 0.1)
logFC_mean$Cluster <- factor(logFC_mean$Cluster, levels=rev(unique(cluster$Cluster)))
logFC_mean$Var2    <- factor(logFC_mean$Var2,    levels=treatment_1$names)
levels(logFC_mean$Var2) <- treatment_1$ID



# =============== for plotting logFC =============== #
# manually create a data frame based on stats/AME-k_means-zero_centered-9/ame_combined_transcription-summary.txt

cl <- unique(cluster$Cluster)

AME <- rbind(
	data.frame(Cluster=cl, family="WRKYs",           idx=as.numeric(cl %in% c("Cl_05", "Cl_06", "Cl_08", "Cl_09"))),
	data.frame(Cluster=cl, family="ANACs",           idx=as.numeric(cl %in% c("Cl_02", "Cl_03"))),
	data.frame(Cluster=cl, family="AHL12/20/25",     idx=as.numeric(cl %in% c("Cl_06", "Cl_09"))),
	data.frame(Cluster=cl, family="AT3G42860(CCHC)", idx=as.numeric(cl %in% c("Cl_06", "Cl_09"))),
	data.frame(Cluster=cl, family="ATHB",            idx=as.numeric(cl %in% c("Cl_06"))),
	data.frame(Cluster=cl, family="FRS9",            idx=as.numeric(cl %in% c("Cl_06")))
)

AME$Cluster <- factor(AME$Cluster, levels=rev(cl))
AME$family  <- factor(AME$family,  levels=unique(AME$family))

write.table(AME, file=paste(stat, "AME-k_means-zero_centered-9/curated_summary.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")



# =============== plotting =============== #
# logFC
p2 <- ggplot(logFC_mean, aes(y=Cluster, x=Var2, fill=mean)) +
	geom_tile(colour="black", size=0.6, alpha=1) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_dark_green, guide=guide_colourbar(frame.colour="black", frame.linewidth=1)) +
	theme_RTN +
	theme(
		axis.text.x=element_text(angle=90, vjust=.5, hjust=1, colour="black", size=9),
		axis.text.y=element_text(size=9, colour="black"),
		legend.text=element_text(size=8),
		legend.title=element_text(size=10),
		axis.line.y=element_blank(),
		axis.line.x=element_blank(),
		axis.ticks=element_blank()) +
	labs(x="", y="Cluster", fill="logFC")

# MEME_data
p3 <- ggplot(AME, aes(x=family, y=Cluster, fill=factor(idx))) +
	geom_tile(colour="black", size=0.6, alpha=1) +
	scale_fill_manual(values=c("1"="black", "0"="white"), labels=c("1"="+", "0"="-"), guide="legend") +
	theme_RTN +
	theme(
		axis.text.x=element_text(angle=90, vjust=.5, hjust=1, colour="black", size=9),
		axis.text.y=element_blank(),
		legend.text=element_text(size=8),
		legend.title=element_text(size=10),
		axis.line.y=element_blank(),
		axis.line.x=element_blank(),
		axis.ticks=element_blank()) +
	labs(x="Enriched motifs", y="", fill="")

p <- p2 + p3 + plot_layout(width=c(1, 1)) & theme(plot.background = element_blank())
ggsave(p, file=paste(fig, "Fig2A-AME_cluster-mean_logFC.pdf", sep=""), height=4.5, width=4.5, bg="transparent")





