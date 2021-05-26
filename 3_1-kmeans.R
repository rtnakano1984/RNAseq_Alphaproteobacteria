
# R script for RNAseq projects 3362 and 3598
# k mean clustering based on log2cpm; Fig S1
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,      quietly=T, warn.conflicts=F)
library(dplyr,        quietly=T, warn.conflicts=F)
library(reshape2,     quietly=T, warn.conflicts=F)
library(ggplot2,      quietly=T, warn.conflicts=F)
library(multcompView, quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# functions
# Compute AIC and BIC
kmeansAIC <- 

# data import
design    <-            read.table(paste(original_data,  "design.txt", sep=""),            sep="\t", header=T,              check.names=F, stringsAsFactors=F)
# log2cpm   <-  as.matrix(read.table(paste(processed_data, "log2cpm_DEG.txt", sep=""),       sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
logFC     <-  as.matrix(read.table(paste(processed_data, "logFC_DEG.txt", sep=""),         sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
norm      <-  as.matrix(read.table(paste(processed_data, "zero_centered_DEG.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))


# =============== appropriate k search based on BIC =============== #
min <- 1
max <- 50
range <- c(min:max)

AIC <- t(sapply(range, function(x){

	set.seed(0)
	fit <- kmeans(norm, centers=x)

	m <- ncol(fit$centers)
	n <- length(fit$cluster)
	k <- nrow(fit$centers)
	D <- fit$tot.withinss
	return(c(
		AIC=(D + 2*m*k),
		BIC=(D + log(n)*m*k)))
}))
AIC <- data.frame(n=range, AIC)

# plot AIC/BIC
melt <- melt(AIC, id.vars="n")
p <- ggplot(melt, aes(x=n, y=value, colour=variable)) +
	geom_line() +
	geom_point() +
	labs(colour="", x="Number of clusters", y="") +
	theme_RTN +
	theme(legend.position="top")
ggsave(p, file=paste(stat, "AIC_BIC-zero_centered.pdf", sep=""), width=4.5, height=3.5)


# output best k
k <- range[which.min(AIC$BIC)]
out <- data.frame("#AIC/BIC Likelihodd test",
	paste("#Tested range of clusters:  from ", min, " to ", max, sep=""),
	"#Best k (with lowest BIC) is:      ",
	k)
write.table(out, file=paste(stat, "AIC_BIC_best_k.zero_centered.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\n")



# =============== kmeans clustering and plotting heatmaps =============== #
# k means clustering with the best k
set.seed(0)
fit <- kmeans(norm, centers=k)

sink(paste(stat, "k_means_summary-zero_centered-", k, ".txt", sep=""))
	print(fit)
sink()

# mean log2cpm score for each cluster
melt_norm <- melt(norm)

idx <- match(melt_norm$Var1, names(fit$cluster))
melt_norm$cluster <- fit$cluster[idx]

melt_norm_mean <- melt_norm %>% group_by(Var2, cluster) %>% summarise(mean=mean(value)) %>% data.frame
norm_mean <- dcast(melt_norm_mean, Var2 ~ cluster, value.var="mean")

rownames(norm_mean) <- norm_mean$Var2
norm_mean <- as.matrix(norm_mean[, -1])

# hclust based on the corerlation of mean Z between clusters
cor <- cor(norm_mean)

hclust <- hclust(as.dist(1-cor), "ward.D2")
sorted_clusters <- rownames(cor)[hclust$order]

pdf(file=paste(stat, "dendrogram_hclust_of_k_mean-zero_centered-clusters_", k, ".pdf", sep=""))
	plot(hclust)
dev.off()

# sort according to sorted clusters
cluster <- data.frame(ID=names(fit$cluster), Cluster=factor(fit$cluster, sorted_clusters))
idx <- order(cluster$Cluster)
cluster <- cluster[idx,]

# rename clusters
levels(cluster$Cluster) <- paste("Cl", str_pad(1:k, 2, "0", side="left"), sep="_")
write.table(cluster, file=paste(stat, "k_means-zero_centered-", k, ".sorted_hclust.txt", sep=""), sep="\t", quote=F, col.names=T, row.names=F)

# export csv for metascape
max_N <- max(table(cluster$Cluster))
metascape <- lapply(levels(cluster$Cluster), function(x){
	idx <- cluster$Cluster == x
	target <- as.character(cluster$ID[idx])
	out <- rep(NA, max_N)
	out[1:length(target)] <- target
	return(out)
}) %>% do.call(cbind, .)
colnames(metascape) <- levels(cluster$Cluster)
write.table(metascape, file=paste(stat, "k_means-zero_centered-", k, ".sorted_hclust.for_metascape.csv", sep=""), sep=",", quote=F, col.names=T, row.names=F)



# =============== compare logFCs among clusters =============== #
melt_logFC <- melt(logFC)

idx <- match(melt_logFC$Var1, cluster$ID)
melt_logFC$cluster <- cluster$Cluster[idx]

melt_logFC$Var2 <- factor(melt_logFC$Var2, levels=treatment_1$names[-1])
levels(melt_logFC$Var2) <- treatment_1$ID[-1]

melt_logFC$value_colour <- saturate(melt_logFC$value, .25)

# boxplot
p <- ggplot(melt_logFC, aes(x=Var2, y=value, colour=value_colour)) +
	geom_hline(yintercept=0, size=1, colour=c_grey, alpha=.5) +
	geom_point(position=position_jitter(width=.25), size=.5) +
	geom_boxplot(fill=c_white, outlier.shape=NA, colour="black", alpha=1) +
	scale_colour_gradient2(low=c_cudo_magenta, mid=c_black, high=c_dark_green, midpoint=0, guide=F) +
	facet_grid( ~ cluster, scales="free") +
	labs(x="", y="logFC") +
	theme_RTN +
	theme(legend.position="none",
		axis.text.x=element_text(colour="black", angle=90, hjust=1, vjust=.5))
ggsave(p, file=paste(fig, "Fig1D-kmeans_boxplot.pdf", sep=""), width=12, height=3.5, bg="transparent")

# stats
kruskal_list <- lapply(levels(cluster$Cluster), function(x){

	idx <- melt_logFC$cluster == x
	temp <- melt_logFC[idx,]

	kruskal <- kruskal.test(value ~ Var2, data=temp)
	if(kruskal$p.value < alpha){

		mean <- temp %>% group_by(Var2) %>% summarise(mean=mean(value))
		sorted_var2 <- as.character(mean$Var2[order(mean$mean)])

		temp$Var2 <- factor(temp$Var2, levels=sorted_var2)

		wilcox <- pairwise.wilcox.test(temp$value, temp$Var2, p.adjust="bonferroni")$p.value
		wilcox <- rbind(matrix(NA, nrow=1, ncol=length(sorted_var2)), cbind(wilcox, matrix(NA, ncol=1, nrow=length(sorted_var2)-1)))

		rownames(wilcox)[1]            <- head(sorted_var2, 1)
		colnames(wilcox)[ncol(wilcox)] <- tail(sorted_var2, 1)

		comb <- expand.grid(sorted_var2, sorted_var2)[,2:1]

		idx <- as.numeric(comb$Var1) > as.numeric(comb$Var2)
		comb <- as.matrix(comb[idx,])
		p.values <- apply(comb, 1, function(x){
			pvalue_1 <- wilcox[x[1], x[2]]
			pvalue_2 <- wilcox[x[2], x[1]]
			pvalue <- sum(c(pvalue_1, pvalue_2), na.rm=T)
			return(pvalue)
		})
		names(p.values) <- apply(comb, 1, paste, collapse="-")

		p_letters <- multcompLetters(p.values, reverse=T)$Letters
		p_letters_df <- as.data.frame(p_letters[treatment_1$ID[-1]])
		write.table(p_letters_df, file=paste(stat, "Fig1D-kmeans_boxplot-kruskal/Fig1-kmeans_boxplot-wilcoxon", x, ".txt", sep=""), col.names=F, row.names=T, quote=F, sep="\t")
	}
	sink(file=paste(stat, "Fig1D-kmeans_boxplot-kruskal/Fig1-kmeans_boxplot-kruskal.", x, ".txt", sep=""))
	print(kruskal)
	sink()
	return(kruskal)
})

# one sample t test
melt_logFC$group <- paste(melt_logFC$cluster, melt_logFC$Var2, sep=":")
t_pvalue <- as.data.frame(tapply(melt_logFC$value, melt_logFC$group, function(x) t.test(x, mu=0)$p.value))

t_pvalue <- cbind(t_pvalue, matrix(unlist(str_split(rownames(t_pvalue), ":")), ncol=2, byrow=T))
colnames(t_pvalue) <- c("p.value", "cluster", "treatment_1")

t_pvalue$p.adjust <- p.adjust(t_pvalue$p.value, "bonferroni")
t_pvalue$sig <- as.numeric(t_pvalue$p.adjust < 0.001)

t_pvalue$treatment_1 <- factor(t_pvalue$treatment_1, levels=treatment_1$ID[-1])
t_pvalue <- t_pvalue[order(t_pvalue$cluster, t_pvalue$treatment_1),]

write.table(t_pvalue, file=paste(stat, "Fig1D-cluster_wise_logFC-one_sample_t_test.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

