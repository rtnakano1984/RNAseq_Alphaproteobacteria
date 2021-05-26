# R script for RNAseq projects 3362 and 3598
# GLM analysis; Fig S1
# 06 Nov 2020
# by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(methods,  quietly=T, warn.conflicts=F)
library(edgeR,    quietly=T, warn.conflicts=F)
library(stringr,  quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/original_data/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# data import
dat_tab <- read.table(paste(original_data, "featureCounts_table-SE.txt", sep=""), sep="\t", stringsAsFactors=F, row.names=1, header=T, check.names=F)
design  <- read.table(paste(original_data, "design.txt", sep=""),                 sep="\t", stringsAsFactors=F,              header=T, check.names=F)

# refine colnames
colnames(dat_tab) <- str_replace(colnames(dat_tab), "_hisat2_sorted.bam", "") %>% str_replace(., ".*hisat2_map\\/", "")

# filter other data
idx <- colnames(dat_tab) %in% design$Library
dat_tab <- dat_tab[, idx]

# sort data table accorindg to AGI code
idx <- order(rownames(dat_tab))
dat_tab <- dat_tab[idx,]

# sort data table according to the design table
idx <- match(design$Library, colnames(dat_tab))
dat_tab <- dat_tab[,idx]


# count list for edgeR
x <- as.matrix(dat_tab)
y <- DGEList(counts=x)

# extract genes with a least 10 counts over all samples AND expressed more than 5 counts at least in two samples 
cpm_count_coef <- mean((y$samples$lib.size) / 1000000)
idx <- rowSums(cpm(y))>(10/cpm_count_coef)  & rowSums(cpm(y)>(5/cpm_count_coef))>=2
y <- y[idx, , keep.lib.sizes=F]


# calculate normalization factor (TMM-normalization)
y <- calcNormFactors(y)

# Create model matrix
Rep       <- factor(design$rep,         levels=c("1", "2", "3"))
group     <- factor(design$treatment_1, levels=treatment_1$names)

model <- model.matrix( ~ 0 + group + Rep)
colnames(model) <- str_replace(colnames(model), "group", "")
colnames(model) <- str_replace(colnames(model), "-", "")

# estimate dispersion
y <- estimateGLMCommonDisp(y, model)
y <- estimateGLMTrendedDisp(y, model)
y <- estimateGLMTagwiseDisp(y, model)

# Generalized linear model fitting
fit <- glmFit(y, model)


# LRT contrasts
contrasts <- makeContrasts(
	Rhi129E     = (     Rhi129E - mock ),
	Rhi13A      = (      Rhi13A - mock ),
	Agro491     = (     Agro491 - mock ),
	Sino142     = (     Sino142 - mock ),
	Sphingo1497 = ( Sphingo1497 - mock ),
	Caulo700    = (    Caulo700 - mock ),
	levels=model)
contrast_names <- colnames(contrasts)
n <- length(contrast_names)

# LRT for each contrasts
LRT.list <- lapply(1:n, function(x) glmLRT(fit, contrast=contrasts[,x]))
names(LRT.list) <- contrast_names

# logFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
	table <- LRT.list[[x]]$table[,c(1,4)]
	table$PValue <- p.adjust(table$PValue, method=p.adj.method)
	colnames(table) <- paste(contrast_names[x], colnames(table), sep="_")
	return(table)
	})
logFC_P <- do.call(cbind, logFC_P.list)
write.table(logFC_P, file=paste(processed_data, "logFC.P.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)

# log2cpm
log2cpm <- cpm(y, prior.count=2, log=TRUE)
write.table(log2cpm, file=paste(processed_data, "log2cpm.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)




# Z-score-like transformation
Z <- t(apply(log2cpm, 1, function(x) {
	m   <- median(x)
	mad <- mad(x, constant=1.4826)	# median absolute deviation
	meanAD <- mean(abs(x-mean(x)))		# mean absolute deviation
	if(mad == 0) {
		y <- (x-m)/(1.253314*meanAD)
	} else {
		y <- (x - m)/mad
	}
	return(y)
}))
write.table(Z, file=paste(processed_data, "mod_Z.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)



# zero-centered expression
norm <- t(apply(log2cpm, 1, function(x) x - mean(x)))
write.table(norm, file=paste(processed_data, "zero_centered.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)




##################### DEG extraction ############################

# Significance picking for each tested model
DE.list <- lapply(1:n, function(x) decideTestsDGE(LRT.list[[x]], adjust.method=p.adj.method, p.value=alpha, lfc=log2(FC_threshold)))
names(DE.list) <- contrast_names

# Number of significant differentially abundant OTUs
total    <- sapply(DE.list, function(x) sum(abs(x)))
induced  <- sapply(DE.list, function(x) sum(x ==  1))
reduced  <- sapply(DE.list, function(x) sum(x == -1))
count <- data.frame(total, induced, reduced)
write.table(count, file=paste(stat, "number_of_DEGs.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

# significance table
DE <- sapply(1:n, function(x) DE.list[[x]][,1])
colnames(DE) <- contrast_names
write.table(DE, file=paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", quote=T, row.names=T, col.names=NA)

# export log2cpm, logFC, pvalue of DEGs
idx <- rowSums(DE != 0) > 0
DEG <- rownames(DE)[idx]
write.table(DEG, file=paste(stat, "DEG_list.txt", sep=""), sep="\n", quote=F, col.names=F, row.names=F)

colnames(logFC_P) <- str_replace(colnames(logFC_P), "_.*", "")
logFC <- logFC_P[, (1:n)*2-1]
P     <- logFC_P[, (1:n)*2]

idx <- rownames(logFC_P) %in% DEG
write.table(logFC[idx,],   file=paste(processed_data, "logFC_DEG.txt", sep=""),   sep="\t", quote=F, row.names=T, col.names=NA)
write.table(P[idx,],       file=paste(processed_data, "P_DEG.txt", sep=""),       sep="\t", quote=F, row.names=T, col.names=NA)

idx <- rownames(log2cpm) %in% DEG
write.table(log2cpm[idx,], file=paste(processed_data, "log2cpm_DEG.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)

idx <- rownames(Z) %in% DEG
write.table(Z[idx,], file=paste(processed_data, "mod_Z_DEG.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)

idx <- rownames(norm) %in% DEG
write.table(norm[idx,], file=paste(processed_data, "zero_centered_DEG.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)







