############################################################################
# Bioramble
# Heatmap: Part 2 - How to create a simple heatmap with R
# by Jesse Lipp
# Jul 31, 2015
############################################################################

# ----------------------------------------------------------
# set up system
# ----------------------------------------------------------

# clear memory
rm(list = ls())

# load packages (install if necessary)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}

# ----------------------------------------------------------
# Step 1: Get to know your data
# ----------------------------------------------------------

# get information on mtcars data set
?mtcars
# load mtcars data set
data(mtcars)
# investigate mtcars data set
dim(mtcars)
class(mtcars)
summary(mtcars)

# ----------------------------------------------------------
# Step 2: Get to know the "heatmap" function
# ----------------------------------------------------------

# convert mtcars data.frame into matrix object
mat <- as.matrix(mtcars)
# plot heatmap using default parameters
heatmap(mat)
# plot heatmap without clustering and scaling
heatmap(mat, Rowv = NA, Colv = NA, scale = "none")
# plot heatmap without clustering but scaling of column features
heatmap(mat, Rowv = NA, Colv = NA, scale = "column")

# clustering of car attributes (rows) using default parameters
row_dist <- dist(mat, method = "euclidean")
row_clus <- hclust(row_dist, method = "complete")
# clustering of cars (columns) using default parameters
col_dist <- dist(t(mat), method = "euclidean")
col_clus <- hclust(col_dist, method = "complete")
# plot heatmap with clustering (default parameters) and scaling of column features
heatmap(mat, 
        # clustering
        Rowv = as.dendrogram(row_clus), 
        Colv = as.dendrogram(col_clus), 
        # scaling
        scale = "column")

# ----------------------------------------------------------
# Step 3: Final tweaks
# ----------------------------------------------------------

# set up custom color palette
yellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
# plot final heatmap
heatmap(mat, 
        # clustering
        Rowv = as.dendrogram(row_clus), 
        Colv = as.dendrogram(col_clus), 
        # scaling
        scale = "column", 
        # color
        col = yellowred)
