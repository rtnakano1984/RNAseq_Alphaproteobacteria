
# R script for RNAseq projects 3362 and 3598
# creating various heatmaps of Rhi-specific DEGs; Fig 4
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(lsa,       quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(patchwork, quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)
library(igraph,  　quietly=T, warn.conflicts=F)
library(dplyr,   　quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))
theme_RTN <- theme_RTN + theme(panel.spacing=unit(3, "pt"))


# data import
design  <- read.table(paste(original_data,  "design.txt",                      sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)
norm    <- read.table(paste(processed_data, "zero_centered_DEG.txt",           sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F) %>% as.matrix
log2cpm <- read.table(paste(processed_data, "log2cpm_DEG.txt",                 sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F) %>% as.matrix
logFC   <- read.table(paste(processed_data, "logFC_DEG.txt",                   sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F) %>% as.matrix
sig     <- read.table(paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F) %>% as.matrix
cluster <- read.table(paste(stat, "k_means-zero_centered-9.sorted_hclust.txt", sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)
go      <- read.table(paste(stat, "topGO-Rhi_curated.txt",                     sep=""), sep="\t", header=F,              check.names=F, stringsAsFactors=F)

# saturate function

# saturate parameter
saturate_cut <- .05


# filtering to Rhi-specific DEGs ##################################
Rhi_sig    <- rowSums(sig[, 1:4]) %in% c(-4, 4)
Sis_nonsig <- rowSums(sig[, 5:6] == 0) == 2
Rhi_specific <- rownames(sig)[Rhi_sig & Sis_nonsig]

# The following gave the same result ##
# Rhi_sig_ind <- rowSums(sig[, 1:4]) ==  4 & rowSums(sig[, 5:6] !=  1) == 2
# Rhi_sig_sup <- rowSums(sig[, 1:4]) == -4 & rowSums(sig[, 5:6] != -1) == 2
# Rhi_specific <- Rhi_sig_ind | Rhi_sig_sup

log2cpm <- log2cpm[rownames(log2cpm) %in% Rhi_specific, ]
norm    <-    norm[rownames(norm)    %in% Rhi_specific, ]
logFC   <-   logFC[rownames(logFC)   %in% Rhi_specific, ]
sig     <-     sig[rownames(sig)     %in% Rhi_specific, ]

sorted <- rownames(logFC)[hclust(as.dist(1 - cosine(t(logFC))), "average")$order]

idx <- rowSums(sig)[match(sorted, rownames(sig))] > 0
Rhi_specific_df <- data.frame(Gene=sorted, Cluster=c("suppressed", "induced")[idx + 1])
write.table(Rhi_specific_df, paste(stat, "Rhi_specific_DEGs.txt", sep=""), sep="\t", quote=F, col.names=T, row.names=F)







# =============== Heatmap plot =============== #

# log2cpm #########################################
melt_log2cpm <- melt(log2cpm)

idx <- match(melt_log2cpm$Var2, design$Library)
melt_log2cpm <- cbind(melt_log2cpm, design[idx,])

melt_log2cpm$Var1        <- factor(melt_log2cpm$Var1,        levels=sorted)
melt_log2cpm$treatment_1 <- factor(melt_log2cpm$treatment_1, levels=treatment_1$names)
levels(melt_log2cpm$treatment_1) <- treatment_1$ID

melt_log2cpm$value <- saturate(melt_log2cpm$value, saturate_cut)

p_log2cpm <- ggplot(melt_log2cpm, aes(x=rep, y=Var1, fill=value)) +
	geom_raster(alpha=1) +
	scale_fill_gradient2(low=c_cudo_skyblue, mid=c_black, high=c_yellow, midpoint=0, guide=guide_colorbar(barwidth=4.5, barheight=.75, frame.colour="black", frame.linewidth=1)) +
	labs(x="", y="", fill="log2cpm") +
	theme_RTN +
	facet_grid(. ~ treatment_1, drop=TRUE, scales="free", space="free", switch="both") +
	theme(legend.position="top",
		legend.text=element_text(size=6),
		legend.title=element_text(size=8),
		axis.text.x=element_blank(),
		axis.text.y=element_text(size=5),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_text(angle=180, size=6),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())






#Z scores #########################################
melt_norm <- melt(norm)

idx <- match(melt_norm$Var2, design$Library)
melt_norm <- data.frame(melt_norm, design[idx,], stringsAsFactors=F)

melt_norm$Var1        <- factor(melt_norm$Var1,        levels=sorted)
melt_norm$treatment_1 <- factor(melt_norm$treatment_1, levels=treatment_1$names) 
levels(melt_norm$treatment_1) <- treatment_1$ID

melt_norm$value <- saturate(melt_norm$value, saturate_cut)

p_norm <- ggplot(melt_norm, aes(x=rep, y=Var1, fill=value)) +
	geom_raster(alpha=1) +
	scale_fill_gradient2(low="#0068B7", mid="white", high="#F39800", midpoint=0, guide=guide_colorbar(barwidth=4.5, barheight=.75, frame.colour="black", frame.linewidth=1)) +
	labs(x="", y="", fill="Zero-centered\nexpression") +
	theme_RTN +
	facet_grid(. ~ treatment_1, drop=TRUE, scales="free", space="free", switch="both") +
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

melt_logFC$Var1 <- factor(melt_logFC$Var1, levels=sorted)
melt_logFC$Var2 <- factor(melt_logFC$Var2, levels=treatment_1$names[-1])
levels(melt_logFC$Var2) <- treatment_1$ID[-1]

melt_logFC$value <- saturate(melt_logFC$value, saturate_cut)

p_logFC <- ggplot(melt_logFC, aes(x=Var2, y=Var1, fill=value)) +
	geom_raster(alpha=1) +
	scale_fill_gradient2(low=c_cudo_magenta, mid="white", high=c_dark_green, midpoint=0, guide=guide_colorbar(barwidth=4.5, barheight=.75, frame.colour="black", frame.linewidth=1)) +
	labs(x="", y="", fill="logFC") +
	facet_grid(. ~ Var2, drop=TRUE, scales="free", space="free", switch="both") +
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

melt_sig$Var1 <- factor(melt_sig$Var1, levels=sorted)
melt_sig$Var2 <- factor(melt_sig$Var2, levels=treatment_1$names[-1])
levels(melt_sig$Var2) <- treatment_1$ID[-1]

idx <- melt_sig$value != 0
melt_sig$value <- as.numeric(idx)

p_sig <- ggplot(melt_sig, aes(x=Var2, y=Var1, fill=factor(value))) +
	geom_raster(alpha=1) +
	scale_fill_manual(values=c("0"="white", "1"="red"), guide=F) +
	labs(x="", y="", fill="") +
	facet_grid(. ~ Var2, drop=TRUE, scales="free", space="free", switch="both") +
	theme_RTN +
	theme(legend.position="top",
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_blank(),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())


# GO assignment #########################################
sorted_go <- unique(go$V2)
go_melt <- lapply(sorted_go, function(x){
	value <- as.numeric(sorted %in% go$V1[go$V2==x])
	out <- data.frame(ID=sorted, GO=x, value=factor(value))
}) %>% do.call(rbind, .)

go_melt$ID <- factor(go_melt$ID, levels=sorted)
go_melt$GO <- factor(go_melt$GO, levels=sorted_go)

p_go <- ggplot(go_melt, aes(x=GO, y=ID, fill=value)) +
	geom_tile(colour="black", size=0.25, alpha=1) +
	scale_fill_manual(values=c("1"="black", "0"="white"), labels=c("1"="+", "0"="-"), guide=F) +
	labs(x="", y="", fill="") +
	facet_grid(. ~ GO, drop=TRUE, scales="free", space="free", switch="both") +
	theme_RTN +
	theme(legend.position="top",
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line=element_blank(),
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.text.y=element_blank(),
		panel.background=element_rect(fill=NA, colour=c_black),
		axis.ticks=element_blank())








# composit and export #########################################
p <- p_log2cpm + p_norm + p_logFC + p_sig + p_go + plot_layout(width=c(2, 2, 1, 1, 2.5)) & theme(plot.background = element_blank())
ggsave(p, file=paste(fig, "Fig3A-heatmaps-Rhi_DEGs.png", sep=""), width=10, height=6, bg="transparent")





# =============== GENIE3 network plot =============== #
g <- readRDS(file=paste(processed_data, "GENIE3_network.RDS", sep=""))
g_layout <- readRDS(file=paste(processed_data, "GENIE3_network_layout.RDS", sep=""))


# mapping Rhi-specific DEGs
genes <- V(g)$name
V(g)$Rhi <- "0"
idx <- genes %in% Rhi_specific_df$Gene[Rhi_specific_df$Cluster == "induced"]
V(g)$Rhi[idx] <- "1"
idx <- genes %in% Rhi_specific_df$Gene[Rhi_specific_df$Cluster == "suppressed"]
V(g)$Rhi[idx] <- "-1"

colours <- c("0"=NA, "-1"=c_cudo_magenta, "1"=c_dark_green)
V(g)$colours <- colours[V(g)$Rhi]

pdf(paste(fig, "Fig3B-GENIE3_RhiDEG.pdf", sep=""), width=12, height=12)
	plot.igraph(g,
		edge.width=E(g)$weight*10,
		vertex.size=V(g)$size,
		vertex.label=NA,
		vertex.color=V(g)$colours,
		layout=g_layout,
		main="Rhizobiales DEGs")
	legend("topright", title="Rhizobiales DEGs",legend=c("Induced", "Suppressed"), col=c(c_dark_green, c_cudo_magenta), pch=20, pt.cex=2.5, bty="n")
dev.off()






# =============== GENIE3 network plot =============== #
g <- readRDS(file=paste(processed_data, "GENIE3_network.DETFs.RDS", sep=""))
g_layout <- readRDS(file=paste(processed_data, "GENIE3_network_layout.DETFs.RDS", sep=""))


# mapping Rhi-specific DEGs
genes <- V(g)$name
V(g)$Rhi <- "0"
idx <- genes %in% Rhi_specific_df$Gene[Rhi_specific_df$Cluster == "induced"]
V(g)$Rhi[idx] <- "1"
idx <- genes %in% Rhi_specific_df$Gene[Rhi_specific_df$Cluster == "suppressed"]
V(g)$Rhi[idx] <- "-1"

colours <- c("0"=NA, "-1"=c_cudo_magenta, "1"=c_dark_green)
V(g)$colours <- colours[V(g)$Rhi]

pdf(paste(fig, "Fig3B-GENIE3_RhiDEG.DETFs.pdf", sep=""), width=12, height=12)
	plot.igraph(g,
		edge.width=E(g)$weight*10,
		vertex.size=V(g)$size,
		vertex.label=NA,
		vertex.color=V(g)$colours,
		layout=g_layout,
		main="Rhizobiales DEGs")
	legend("topright", title="Rhizobiales DEGs",legend=c("Induced", "Suppressed"), col=c(c_dark_green, c_cudo_magenta), pch=20, pt.cex=2.5, bty="n")
dev.off()
















