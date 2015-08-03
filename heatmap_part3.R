# ----------------------------------------------------------
# set up system
# ----------------------------------------------------------

# clear memory
rm(list = ls())

# load packages (install if necessary)
if (!require(ALL)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("ALL")
}
if (!require("genefilter")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("genefilter")
}
if (!require("gplots")) {
  install.packages("gplots")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}

# ----------------------------------------------------------
# Step 1: Prepare the data
# ----------------------------------------------------------

# load accute lymphoblastic leukemia data
data(ALL)
# get information on ALL data set
?ALL
# get basic information on data
class(ALL)
dim(ALL)
# check for NA values
any(is.na(exprs(ALL)))
# data formating following 
# http://www.bioconductor.org/help/publications/books/bioconductor-case-studies/
# subset data set to B-cell samples with genotype "NEG" or "BCR/ABL" translocation
# "BT" is the type of the cancer (B-cell or T-cell)
# "mol.biol" is the genotypic classification of the cancer
# get indices of cancers with either no cytogenetic abnormalities (NEG) or 
# the BCR-ABL translocation (BCR/ABL)
bcrabl <- ALL$mol.biol %in% c("NEG", "BCR/ABL")
# get indices cancers originating from B-cells
bcell <- grepl("^B", ALL$BT)
# subset the ALL data set
all <- ALL[, bcell & bcrabl]
# re-adjust the factor levels to reflect the subset
all$mol.biol <- droplevels(all$mol.biol)
all$mol.biol <- relevel(all$mol.biol, ref = "NEG")
# get dimensions again
dim(all)

# # determine standard deviation of all genes
# all_sd <- apply(exprs(all), 1, sd)
# there is an optimized function in genefilter package "rowSds"
all_sd <- rowSds(exprs(all))
# get names of 200 most variable genes
top200 <- names(sort(all_sd, decreasing = TRUE))[1:200]
all_var <- all[top200, ]

# ----------------------------------------------------------
# Step 2: Decide on a distance metric
# ----------------------------------------------------------

# distance function
# remember: cor computes across columns
dist_cor <- function(x) {
  as.dist(1 - cor(t(x), method = "pearson"))
}

# ----------------------------------------------------------
# Step 3: Decide on a clustering method
# ----------------------------------------------------------

# clustering function
clus_wd2 <- function(x) {
  hclust(x, method = "ward.D2")
}

# ----------------------------------------------------------
# Step 4: Plot a microarray heatmap
# ----------------------------------------------------------

# it is customary to use a red (up-regulated), black (neutral), 
# green (down-regulated) color scheme for expression microarrays
# create a red-black-green color palette
redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)
# generate genotype class labels
# "NEG" is light grey, "BCR/ABL" is dark grey
class_labels <- ifelse(all_var$mol.biol == "NEG", "grey80", "grey20")

# plot microarray data using "heatmap.2" from "gplots" library
heatmap.2(exprs(all_var), 
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = "", 
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")

# ----------------------------------------------------------
# Step 5: A "better" way of selecting genes
# ----------------------------------------------------------

# non-specific filtering
# the shortest interval containing half of the data
# reasonable estimate of the "peak" of the distribution
sh <- shorth(all_sd)
all_sh <- all[all_sd >= sh, ]
dim(all_sh)
# row-wise t-tests using "rowttests" from "genefilter" package
tt <- rowttests(all_sh, all_sh$mol.biol)
# adjust p-values for multiple testing 
# using "Benjamini-Hochberg" method
tt$p.adj <- p.adjust(tt$p.value, method = "BH")
all_sig <- all_sh[tt$p.adj < 0.05, ]
# how many genes are we left with
dim(all_sig)

# plot heatmap for differentially expressed genes
heatmap.2(exprs(all_sig), 
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = "", 
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")


# ----------------------------------------------------------
# Step 6: Have mercy with the color-challenged
# ----------------------------------------------------------

# blue to black to yellow color map
yellowblackblue <- colorRampPalette(c("dodgerblue", "black", "gold"))(n = 100)

heatmap.2(exprs(all_sig), 
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = yellowblackblue, 
          # labels
          labRow = "", 
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")
