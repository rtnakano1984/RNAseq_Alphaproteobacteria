
# R script for RNAseq projects 3362 and 3598
# Comparing responses between treatments; Fig 1
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(reshape2,  quietly=T, warn.conflicts=F)
library(patchwork, quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)
library(lsa,       quietly=T, warn.conflicts=F)
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
log2cpm   <-  as.matrix(read.table(paste(processed_data, "log2cpm_DEG.txt", sep=""),       sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
logFC     <-  as.matrix(read.table(paste(processed_data, "logFC_DEG.txt", sep=""),         sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
norm      <-  as.matrix(read.table(paste(processed_data, "zero_centered_DEG.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))

# PCoA; Fig 1A
# PCC computation based on zero-centered
cor <- cor(norm)

# PCoA
pcoa <- cmdscale(as.dist(1-cor), k=3, eig=T)

# variation explained (eigen vectors)
pro.var <-pcoa$eig[1:3] / sum(pcoa$eig)

# coordinates for data points
points <- pcoa$points
colnames(points) <- c("x", "y", "z")

# merge with meta data
idx <- match(rownames(points), design$Library)
points <- data.frame(design[idx,], points)
points$treatment_1 <- factor(points$treatment_1, levels=treatment_1$names)

# plot
p <- ggplot(points, aes(x=x, y=y, colour=treatment_1, shape=treatment_1)) +
	geom_point(alpha=.75, size=3) +
	scale_colour_manual(values=treatment_1$colours, labels=treatment_1$ID) +
	scale_shape_manual(values=treatment_1$shapes, labels=treatment_1$ID) +
	theme_RTN +
	theme(legend.position="right") +
	labs(x=paste("PC1 (", format(pro.var[1]*100, digits=4), "%)", sep=""),
	     y=paste("PC2 (", format(pro.var[2]*100, digits=4), "%)", sep=""),
	     colour="",
	     shape="") 
ggsave(p, file=paste(fig, "Fig1A-PCoA.zero_centered.DEG.pdf", sep=""), bg="transparent", width=4.5, height=3)


# plot
p <- ggplot(points, aes(x=z, y=y, colour=treatment_1, shape=treatment_1)) +
	geom_point(alpha=.75, size=3) +
	scale_colour_manual(values=treatment_1$colours, labels=treatment_1$ID) +
	scale_shape_manual(values=treatment_1$shapes, labels=treatment_1$ID) +
	theme_RTN +
	theme(legend.position="right") +
	labs(x=paste("PC3 (", format(pro.var[3]*100, digits=4), "%)", sep=""),
	     y=paste("PC2 (", format(pro.var[2]*100, digits=4), "%)", sep=""),
	     colour="",
	     shape="") 
ggsave(p, file=paste(fig, "Fig1A-PCoA2.zero_centered.DEG.pdf", sep=""), bg="transparent", width=4.5, height=3)




# Heatmap; Fig 1C
# PCC computation based on logFC
logFC_cor <- cosine(logFC)

# reshape data
logFC_cor[upper.tri(logFC_cor, diag=F)] <- NA
melt <- melt(logFC_cor)
melt <- melt[!is.na(melt$value),]

melt$Var1 <- factor(melt$Var1, levels=treatment_1$names[-1])
melt$Var2 <- factor(melt$Var2, levels=treatment_1$names[-1])

levels(melt$Var1) <- treatment_1$ID[-1]
levels(melt$Var2) <- treatment_1$ID[-1]

melt$label <- format(round(melt$value, digits=2), nsmall=2)

# plot
p <- ggplot(melt, aes(x=Var1, y=Var2, fill=value, label=label)) +
	geom_tile() +
	geom_text(colour=c_black, size=2.5) +
	scale_fill_gradient2(mid="black", high="yellow", 
		midpoint=0) +
	scale_x_discrete(limits=rev(treatment_1$ID[-1])) +
	scale_y_discrete(limits=treatment_1$ID[-1]) +
	theme_RTN +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		axis.line.x=element_blank(),
		axis.line.y=element_blank()) +
	labs(x="", y="", fill="Cosine similarity")
ggsave(p, file=paste(fig, "Fig1B-logFC_cosine-heatmap.pdf", sep=""), bg="transparent", height=4.5, width=4)



# # Heatmap of logFC; Fig 1C
# # reshape data
# melt <- melt(as.matrix(logFC))
# melt$Var2 <- factor(melt$Var2, levels=treatment_1$names)
# levels(melt$Var2) <- treatment_1$ID

# # sort by logFC R129_E
# idx <- order(logFC[, colnames(logFC) == "Rhi129E"])
# sort_R129E <- rownames(logFC)[idx]
# melt$Var1 <- factor(melt$Var1, levels=sort_R129E)

# # saturate for heatmap
# melt$fill <- melt$value

# minq <- quantile(melt$value, .01)
# idx <- melt$value < minq
# melt$fill[idx] <- minq

# maxq <- quantile(melt$value, .99)
# idx <- melt$value > maxq
# melt$fill[idx] <- maxq

# # for axis scaling
# min <- min(melt$value)
# max <- max(melt$value)

# # for filtering strains
# idx <- treatment_1$ID %in% c("Root1497", "Root700")

# # plot
# p1 <- ggplot(melt, aes(x=Var2, y=Var1, fill=fill)) +
# 	geom_raster(alpha=1) +
# 	scale_fill_gradient2(low=c_cudo_magenta, mid="white", high=c_dark_green, midpoint=0) +
# 	labs(x="", y="", fill="logFC") +
# 	facet_grid(. ~ Var2, drop=TRUE, scales="free", space="free", switch="x") +
# 	theme_RTN +
# 	theme(legend.position="top",
# 		axis.text.x=element_blank(),
# 		axis.text.y=element_blank(),
# 		axis.line.x=element_blank(),
# 		axis.line.y=element_blank(),
# 		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
# 		strip.text.y=element_text(angle=180, size=6),
# 		panel.background=element_rect(fill=NA, colour=c_black),
# 		axis.ticks=element_blank())

# p2 <- ggplot(melt[idx,], aes(x=value, y=Var1)) +
# 	geom_vline(xintercept=0, size=.2) +
# 	geom_point(colour=c_black, size=.25, alpha=.3) +
# 	labs(x="Rhizobiales", y="", shape="") +
# 	xlim(min, max) +
# 	theme_RTN +
# 	theme(axis.text.y=element_blank(),
# 		axis.title.x=element_text(size=7),
# 		axis.line.y=element_blank(),
# 		axis.ticks.y=element_blank(),
# 		legend.position="none")

# p3 <- ggplot(melt[!idx,], aes(x=value, y=Var1)) +
# 	geom_vline(xintercept=0, size=.2) +
# 	geom_point(colour=c_black, size=.25, alpha=.3) +
# 	labs(x="Sister lineage", y="", shape="") +
# 	xlim(min, max) +
# 	theme_RTN +
# 	theme(axis.text.y=element_blank(),
# 		axis.title.x=element_text(size=7),
# 		axis.line.y=element_blank(),
# 		axis.ticks.y=element_blank(),
# 		legend.position="none")

# p <- (p1 | p2 | p3) + plot_layout(width=c(4, 1, 1)) & theme(plot.background = element_blank())
# ggsave(p, file=paste(fig, "Fig1C-logFC_compare-R129_E.png", sep=""), width=4.5, height=6.5, bg="transparent")













