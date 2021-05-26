
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
library(GENIE3, 　 quietly=T, warn.conflicts=F)
library(igraph, 　 quietly=T, warn.conflicts=F)
library(patchwork, quietly=T, warn.conflicts=F)
library(ggrepel,   quietly=T, warn.conflicts=F)


#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
fig_temp       <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/igraph_temp/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
design   <- read.table(paste(original_data,  "design.txt",  sep=""),                     sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
norm     <- read.table(paste(processed_data, "zero_centered_DEG.txt", sep=""),           sep="\t", header=T, row.names=1,    check.names=F, stringsAsFactors=F) %>% as.matrix
logFC    <- read.table(paste(processed_data, "logFC_DEG.txt", sep=""),                   sep="\t", header=T, row.names=1,    check.names=F, stringsAsFactors=F) %>% as.matrix
AHL      <- read.table(paste(stat, "AME-k_means-zero_centered-9/AHL_list.txt", sep=""),  sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
TF       <- read.table("/netscratch/dep_psl/grp_psl/ThomasN/resources/Ath_TF_list",      sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
cluster  <- read.table(paste(stat, "k_means-zero_centered-9.sorted_hclust.txt", sep=""), sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)

# filter ANACs and WRKYs
families <- data.frame(
	names=c("WRKY", "NAC", "HSF", "bHLH", "MYB", "MYB_related", "HD-ZIP", "YABBY", "AHL", "Others"),
	colours=c(c_dark_green, c_cudo_magenta, c_dark_brown, c_cudo_skyblue, c_blue, c_blue, c_dark_red, c_red, c_dark_yellow, c_grey),
	stringsAsFactors=F)

# families <- data.frame(
# 	names=c("WRKY", "NAC", "HD-ZIP", "AHL"),
# 	colours=c(c_dark_green, c_cudo_magenta, c_blue, c_dark_red),
# 	stringsAsFactors=F)



# extract TF names
idx1 <- TF$Gene_ID %in% rownames(logFC)
idx2 <- AHL$Gene_ID %in% rownames(logFC)
regulators <- rbind(TF[idx1, c("Gene_ID", "Family")], AHL[idx2, c("Gene_ID", "Family")]) %>% unique

idx <- regulators$Family %in% families$names
regulators$Family[!idx] <- "Others"

# filter to DEGs
idx <- regulators$Gene_ID %in% rownames(norm)
regulators <- regulators[idx,]

# stats
tab <- table(regulators$Family)
write.table(as.data.frame(tab), file=paste(stat, "AME-k_means-zero_centered-9/DEG_TFs.txt", sep=""), col.names=F, row.names=T, quote=F, sep="\t")

# refine family df
idx <- families$names %in% regulators$Family
families <- families[idx,]

idx <- is.na(families$Family)
families$Family <- "Others"


# creating a weighted matrix
set.seed(0)
weight_mat <- GENIE3(norm, regulators$Gene_ID, treeMethod="RF", nCores=40, verbose=T)

# extracting edges
linkList   <- getLinkList(weight_mat, reportMax=3000)
g <- graph.data.frame(linkList, directed=F)
E(g)$weight <- linkList$weight

# plot layout
g_layout <- layout.fruchterman.reingold(g, niter=10000, area=100*vcount(g)^2)

# assign TF families
for(x in families$names){
	idx <- V(g)$name %in% regulators$Gene_ID[regulators$Family==x]
	V(g)$group[idx] <- x
}

# assign colours
colours <- families$colours
names(colours) <- families$names
V(g)$colours <- colours[V(g)$group]

# assign labels only for regulators
V(g)$labels <- ""
idx <- match(V(g)$name, regulators$Gene_ID)
V(g)$labels <- regulators$Gene_ID[idx]

# larger size for regulators
idx <- is.na(V(g)$group)
V(g)$size[idx]  <- 2
V(g)$size[!idx] <- 4


saveRDS(g,        file=paste(processed_data, "GENIE3_network.DETFs.RDS", sep=""))
saveRDS(g_layout, file=paste(processed_data, "GENIE3_network_layout.DETFs.RDS", sep=""))

# g        <- readRDS(file=paste(processed_data, "GENIE3_network.RDS", sep=""))
# g_layout <- readRDS(file=paste(processed_data, "GENIE3_network_layout.RDS", sep=""))

# plot
pdf(paste(fig, "Fig4B-GENIE3_all.DETFs.pdf", sep=""), width=12, height=12)
	# plot.igraph(g1,
	# 	edge.width=E(g1)$weight*10,
	# 	vertex.size=3,
	# 	vertex.label=V(g1)$labels,
	# 	vertex.label.size=2,
	# 	vertex.color=V(g1)$colours,
	# 	layout=g1_layout,
	# 	main="GENIE3_ALL")
	plot.igraph(g,
		edge.width=E(g)$weight*10,
		vertex.size=V(g)$size,
		vertex.label=V(g)$labels,
		vertex.label.size=2,
		vertex.color=V(g)$colours,
		layout=g_layout,
		main="GENIE3_WRKY_NAC")
	legend("topright", title="Regulators", legend=names(colours), col=colours, pch=20, pt.cex=3, bty="n")

	plot.igraph(g,
	edge.width=E(g)$weight*10,
	vertex.size=V(g)$size,
	vertex.label=NA,
	vertex.color=V(g)$colours,
	layout=g_layout,
	main="GENIE3_WRKY_NAC")
	legend("topright", title="Regulators", legend=names(colours), col=colours, pch=20,pt.cex=3, bty="n")

dev.off() 



# mapping log2FC on g
# creating a colour palette
pal <- colorRampPalette(c(c_cudo_magenta, c_black, c_green))
breaks <- 7

# melt
melt_FC  <- melt(logFC)

# renname treatment_1
idx <- match(melt_FC$Var2, treatment_1$names)
melt_FC$Var2 <- treatment_1$ID[idx]

# saturate logFC for plotting
melt_FC$value <- saturate(melt_FC$value, 0.01)

# make 0 at the midpoint
lim <- max(abs(max(melt_FC$value)), abs(min(melt_FC$value)))
cut <- cut(c(-lim, melt_FC$value, lim), breaks=breaks)
cut <- cut[2:(length(cut)-1)]

# make legend labels
cut_lab <- levels(cut)
lim_round <- format(lim, digits=3)
cut_lab[1]      <- str_replace(cut_lab[1],      lim_round, paste(lim_round, "<", sep=""))
cut_lab[breaks] <- str_replace(cut_lab[breaks], lim_round, paste(">", lim_round, sep=""))

# assign color according to logFC
melt_FC$color <- pal(breaks)[as.numeric(cut)]

# plotting
pdf(paste(fig, "Fig2B-GENIE3_logFC.DETFs.pdf", sep=""), width=12, height=12)
	genes <- V(g)$name
	i <- 0
	for(x in treatment_1$ID[-1]){
		i <- i + 1
		idx <- melt_FC$Var2 == x
		temp <- melt_FC[idx,]

		idx <- match(genes, temp$Var1)
		V(g)$colours <- temp$color[idx]
			plot.igraph(g,
				edge.width=E(g)$weight*10,
				vertex.size=4,
				vertex.frame.color=NA,
				vertex.label=NA,
				vertex.color=V(g)$colours,
				layout=g_layout)
			title(x)
			if(i == 1) legend("topleft", title="logFC",legend=rev(cut_lab), col=rev(pal(breaks)), pch=20, pt.cex=2.5)
	}
dev.off()






# mapping cluster on g
# assign cluster
for(x in cluster$Cluster){
	idx <- V(g)$name %in% cluster$ID[cluster$Cluster == x]
	V(g)$group[idx] <- x
}

# assign colours
colours <- c(c_dark_red, c_red, c_cudo_magenta, c_green, c_blue, c_dark_green, c_dark_brown, c_very_dark_green, c_cudo_skyblue)
names(colours) <- unique(cluster$Cluster)
V(g)$colours <- colours[V(g)$group]


# plot
pdf(paste(fig, "FigS2-GENIE3-k_means-zero_centered-9.DETFs.pdf", sep=""), width=12, height=12)

	plot.igraph(g,
	edge.width=E(g)$weight*10,
	vertex.size=V(g)$size,
	vertex.label=NA,
	vertex.color=V(g)$colours,
	layout=g_layout,
	main="GENIE3_WRKY_NAC")
	legend("topright", title="k-means clusters", legend=names(colours), col=colours, pch=20,pt.cex=3, bty="n")

	plot.igraph(g,
		edge.width=E(g)$weight*10,
		vertex.size=V(g)$size,
		vertex.label=V(g)$labels,
		vertex.label.size=2,
		vertex.color=V(g)$colours,
		layout=g_layout,
		main="GENIE3_WRKY_NAC")
	legend("topright", title="k-means clusters", legend=names(colours), col=colours, pch=20, pt.cex=3, bty="n")

dev.off() 


# assign colours
colours <- c(c_dark_red, c_red, c_cudo_magenta, c_green, c_blue, c_dark_green, c_dark_brown, c_very_dark_green, c_cudo_skyblue)
names(colours) <- unique(cluster$Cluster)
V(g)$colours <- colours[V(g)$group]








# network properties
degree  <- degree(g, mode="in")
close   <- closeness(g, mode="all", weights=NA)
between <- betweenness(g, directed=F, weights=NA)

centrality <- cbind(as.data.frame(degree), close) %>% cbind(between)
centrality$ID <- rownames(centrality)

# regulators
idx <- centrality$ID %in% regulators$Gene_ID
centrality <- centrality[idx,]

idx <- match(centrality$ID, regulators$Gene_ID)
centrality$Family <- str_replace(regulators$Family[idx], "MYB_related", "MYB")
centrality$Family <- factor(centrality$Family, levels=families$names)

melt <- melt(centrality, id.vars=c("ID", "Family"))
melt$variable <- factor(melt$variable, levels=c("degree", "close", "between"))
levels(melt$variable) <- c("Degree", "Closeness", "Betweenness")

# for text labels
idx <- centrality$degree > 75
centrality$label[idx]  <- centrality$ID[idx]

# scatter plot
p1 <- ggplot(centrality, aes(x=degree, y=close, colour=Family, label=label)) +
	geom_point() +
	scale_colour_manual(values=families$colours, guide=F) +
	labs(x="Degree", y="Closeness", colour="") +
	theme_RTN
p2 <- ggplot(centrality, aes(x=degree, y=between, colour=Family, label=label)) +
	geom_point() +
	scale_colour_manual(values=families$colours, guide=F) +
	labs(x="Degree", y="Betweenness", colour="") +
	theme_RTN
p <- p1 + p2 + plot_layout(ncol=2) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "Fig2C-GENIE3_centrality_scatter.DETFs.pdf", sep=""), width=6, height=2.5, bg="transparent")

p1 <- p1 + geom_text_repel(segment.size=.5, size=3, show.legend=F)
p2 <- p2 + geom_text_repel(segment.size=.5, size=3, show.legend=F)
p <- p1 + p2 + plot_layout(ncol=2) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "Fig2C-GENIE3_centrality_scatter.DETFs.label.pdf", sep=""), width=6, height=2.5, bg="transparent")


# boxplot
p3 <- ggplot(melt, aes(x=Family, y=value, colour=Family)) +
	geom_boxplot(outlier.shape=NA, fill=NA) +
	geom_point(position=position_jitterdodge(jitter.width=.75)) +
	scale_colour_manual(values=families$colours, guide=F) +
	facet_wrap( ~ variable, scales="free", nrow=1) +
	theme_RTN +
	labs(x="", y="Centrality", colour="") +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(p3, file=paste(fig, "Fig2C-GENIE3_centrality_boxplot.DETFs.pdf", sep=""), width=6, height=3, bg="transparent")







