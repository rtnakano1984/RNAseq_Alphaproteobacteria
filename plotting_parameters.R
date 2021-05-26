
require(ggplot2)
source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")

# parameter
p.adj.method <- "fdr"
FC_threshold <- 2
alpha <- 0.05

# colour
treatment_1 <- data.frame(names=c("mock", "Rhi129E", "Rhi13A", "Agro491", "Sino142", "Sphingo1497", "Caulo700"),
	ID=c("axenic", "R129_E", "R13_A", "Root491", "Root142", "Root1497", "Root700"),
  group=c("Control", "Rhizobium spp.", "Rhizobium spp.", "Agrobacterium", "Sinorhizobium", "Sphingomonadales", "Caulobacterales"),
  shapes=c(8, 15, 16, 17, 18, 1, 2),
  stringsAsFactors=F)
idx <- match(treatment_1$group, rhizobia_colour$Taxonomy)
treatment_1$colours <- rhizobia_colour$Colour[idx]

num.rep = 3

saturate <- function(value, cut){
	max <- quantile(value, (1-cut))
	idx <- value > max
	value[idx] <- max

	min <- quantile(value, cut)
	idx <- value < min
	value[idx] <- min

	return(value)
}
