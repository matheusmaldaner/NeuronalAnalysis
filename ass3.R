# CHANGE PATH VARIABLE
setwd("C:/Users/twinb/OneDrive/Desktop/school/uf_fall_23/cgs4144/BioinformaticsProject/")

library(readr)
library(dplyr)
library(matrixStats)

## Grab 5000 most variable genes -----------------------------
# load diff expr results file
differential_expression_df <- readr::read_tsv("results/SRP094496_diff_expr_results.tsv")
head(differential_expression_df)

# store gene IDs in vars and drop them from main df
# this is so rowVars works (doesn't work with columsn of val <char>)
symbol <- differential_expression_df$symbol
Ensembl <- differential_expression_df$Ensembl
differential_expression_df$symbol <- NULL
differential_expression_df$Ensembl <- NULL

# calculate variances for each gene
gene_variability <- rowVars(as.matrix(differential_expression_df))

# bind the variances and gene IDs
unsorted_gene_vars <- cbind(1:length(symbol), symbol, gene_variability)
sorted_gene_vars <- unsorted_gene_vars[order(as.numeric(unsorted_gene_vars[, 3]), decreasing=TRUE), ] 
sorted_gene_vars

# drop NA vals
sorted_gene_vars_no_na <- sorted_gene_vars[complete.cases(sorted_gene_vars), ]
sorted_gene_vars_no_na

# grab first 5000
most_var_5000 <- sorted_gene_vars_no_na[1:5000, ]
most_var_5000

data_to_cluster <- differential_expression_df[most_var_5000[,1], ]
# ----------------------------------------------------------------------------

# Clustering

# K MEANS

k <- 3
kmeans_result <- kmeans(data_to_cluster, centers = k)
kmeans_result


# cluster assignments for each sample
kmeans_cluster_assignments <- kmeans_result$cluster

table(kmeans_cluster_assignments)

cluster_results <- cbind(1:5000, kmeans_cluster_assignments)


# ----------------------------------------------------------------------------

# HIERARCHICAL CLUSTERING

# Calculate the distance matrix
dist_matrix <- dist(data_to_cluster, method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "complete")
plot(hclust_result, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_cluster_assignments <- cutree(hclust_result, k=10)

cluster_results <- cbind(cluster_results, hclust_cluster_assignments)
table(hclust_cluster_assignments)

# test with 10 genes
hc <- hclust(dist(data_to_cluster[1:10, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
table(cutree(hc,k=10))

# test with 100 genes
hc <- hclust(dist(data_to_cluster[1:100, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
table(cutree(hc,k=10))

# test w 1000 genes
hc <- hclust(dist(data_to_cluster[1:1000, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
table(cutree(hc,k=10))


## sankey/alluvial plot





# ----------------------------------------------------------------------------

# CONSESUS CLUSTER PLUS








# ----------------------------------------------------------------------------

# GAUSSIAN MIXTURE MODELS

# installs and loads mclust package
install.packages("mclust")
library(mclust)


# performs GMM clustering
gmm_result <- Mclust(data_to_cluster)


# prints the clustering results
summary(gmm_result)


# number of clusters identified
num_clusters <- gmm_result$G
cat("Number of clusters identified by GMM:", num_clusters, "\n")


# cluster assignments for each sample
gmm_cluster_assignments <- gmm_result$classification
table(gmm_cluster_assignments)

# ----------------------------------------------------------------------------

# RANDOM STUFF - JUST PLOTTING A PCA TO VISUALIZE CLUSTERS

install.packages(c("ggplot2", "FactoMineR"))
library(ggplot2)
library(FactoMineR)

# performs PCA
data_matrix <- as.matrix()
pca_result <- PCA(data_matrix, graph = FALSE)


# extracts first two principal components
pca_data <- data.frame(
  PC1 = pca_result$ind$coord[, 1],
  PC2 = pca_result$ind$coord[, 2],
  Cluster = factor(kmeans_cluster_assignments) # WE CAN CHANGE THIS TO DIFF TYPES
)

# plot using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(title = "PCA Plot of Clusters", x = "PC1", y = "PC2") +
  scale_color_discrete(name = "Cluster")


# ----------------------------------------------------------------------------

# HEATMAP

# Install and load necessary packages
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("dendextend", quietly = TRUE)) {
  install.packages("dendextend")
}

library(pheatmap)
library(dendextend)

# Load your data
# Replace 'your_data_file.csv' with the actual file path to your gene expression data
data <- data_to_cluster

data

# Create clustering data
# Replace 'cluster_data.csv' with the actual file path to your clustering results
cluster_data <- as.data.frame(cluster_results)

# Create a heatmap
heatmap_data <- data  # Select the first 5000 genes
sample_groups <- cluster_data$ # Sample groups from Assignment 1

# Create row and column dendrograms
row_dend <- as.dendrogram(hclust(dist(heatmap_data)))
col_dend <- as.dendrogram(hclust(dist(t(heatmap_data))))

# Create annotation data
annotation_data <- data.frame(
  Method1 = cluster_data$kmeans_cluster_assignments,
  Method2 = cluster_data$hclust_cluster_assignments
)

# Create the heatmap
pheatmap(
  heatmap_data,
  row_dendrogram = row_dend,
  col_dendrogram = col_dend,
  annotation_col = annotation_data,
  show_colnames = FALSE,  # You can set this to TRUE if you want to display column names
  legend = TRUE,
  main = "Heatmap of 5,000 Genes",
  filename = "heatmap.png"  # You can specify the file name and format
)


# different types, just label yours...


