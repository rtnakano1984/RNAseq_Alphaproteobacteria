
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
library(stringr,   quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))
# theme_RTN <- theme_RTN + theme(panel.spacing=unit(3, "pt"))


# data import
design  <-           read.table(paste(original_data,  "design.txt",                      sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)
go_dat  <-           read.table(paste(fig, "all.ti_pkevnz/Enrichment_GO/GO_AllLists-sub.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")

# gene_list <- lapply(1:nrow(go_dat), function(x) c(unlist(str_split(go_dat$GeneID[x], "\\|"))) )
# names(gene_list) <- str_replace_all(paste(go_dat$GeneList, go_dat$Description, sep="_"), " ", "_")

dat <- lapply(1:nrow(go_dat), function(x){
	genes <- c(unlist(str_split(go_dat$GeneID[x], "\\|")))
	out <- data.frame(go_dat[rep(x, length(genes)), 1:2], GeneID=genes, stringsAsFactors=F)
}) %>% do.call(rbind, .) %>% unique

p_list <- lapply(unique(dat$GeneList), function(x){
	idx <- dat$GeneList == x
	temp <- dat[idx, ]

	dcast <- dcast(temp, Description ~ GeneID, value.var="GeneID", length)

	# sort genes
	cor <- cor(as.matrix(dcast[, 2:ncol(dcast)]))
	cor[is.na(cor)] <- 1
	hclust <- hclust(as.dist(1-cor), "ward.D2")
	sorted_geneID <- hclust$label[hclust$order]
	temp$GeneID <- factor(temp$GeneID, levels=sorted_geneID)

	# sort GOs
	cor <- cor(t(as.matrix(dcast[, 2:ncol(dcast)])))
	cor[is.na(cor)] <- 1
	hclust <- hclust(as.dist(1-cor), "ward.D2")
	sorted_GO <- as.character(dcast$Description)[hclust$order]
	temp$Description <- factor(temp$Description, levels=sorted_GO)

	p <- ggplot(temp, aes(x=GeneID, y=Description)) +
		geom_tile(fill=c_black) +
		theme_RTN +
		theme(axis.text.y=element_text(size=4),
			axis.text.x=element_blank()) +
		labs(x="Genes", y="GO terms", title=x)

	return(p)
})

p <- wrap_plots(p_list, ncol=3, nrow=3)
ggsave(p, file=paste(fig, "all.ti_pkevnz/Enrichment_GO/GO_AllLists-gene_profile.pdf", sep=""), width=12, height=12, bg="transparent")



