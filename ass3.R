# CHANGE PATH VARIABLE
#setwd("C:/Users/Matheus/OneDrive/University of Florida/JUNIOR FALL/CGS4144/project")

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
most_var_10000 <- sorted_gene_vars_no_na[1:10000, ]
most_var_5000

# ----------------------------------------------------------------------------

# Clustering
install.packages("ggalluvial")
library(cluster)
library(ggplot2)
library(ggalluvial)

# K MEANS
data_to_cluster <- differential_expression_df[most_var_5000[,1], ]
data_to_cluster_10000 <- differential_expression_df[most_var_10000[,1], ]
#Set to 3 clusters
k <- 3
kmeans_result <- kmeans(data_to_cluster, centers = k)
kmeans_result

# cluster assignments for each sample
kmeans_cluster_assignments <- kmeans_result$cluster

table(kmeans_cluster_assignments)

#Set different k values
k_values <- c(2, 3, 4, 5, 6)

kmeans_results <- list()

# Run K-means for each k and store the results
for (k in k_values) {
  set.seed(123) # Set seed for reproducibility
  kmeans_results[[paste("k", k, sep = "_")]] <- kmeans(data_to_cluster, centers = k)
}

# Calculate silhouette widths for each k
sil_widths <- sapply(kmeans_results, function(km) {
  silhouette_scores <- silhouette(km$cluster, dist(data_to_cluster))
  mean(silhouette_scores[, 3]) # average silhouette width
})

# Plot average silhouette width for each k
plot(k_values, sil_widths, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters 'k'", ylab = "Average silhouette width",
     main = "Silhouette Analysis of k-means clustering")

#Test different number of genes
k <- 3
gene_values <- c(10, 100, 1000)
kmeans_genes_results <- list()
for (g in gene_values) {
  set.seed(123) # Set seed for reproducibility
  kmeans_genes_results[[paste(g, "genes", sep = "_")]] <- kmeans(data_to_cluster[1:g, ], centers = k)
}
set.seed(123) # Set seed for reproducibility
kmeans_genes_results[[paste("10000", "genes", sep = "_")]] <- kmeans(data_to_cluster_10000, centers = k)

sample_ids <- rownames(data_to_cluster)

alluvial_data <- do.call(rbind, lapply(names(kmeans_genes_results), function(genes) {
  data.frame(
    sample = sample_ids,
    cluster = kmeans_genes_results[[genes]]$cluster,
    genes = genes
  )
}))
# Convert 'genes' to factor and specify levels/ordering if needed
alluvial_data$genes <- factor(alluvial_data$genes, levels = c("10_genes", "100_genes", "1000_genes", "10000_genes"))

# Plot
ggplot(data = alluvial_data,
       aes(axis1 = genes, axis2 = cluster)) +
  geom_alluvium(aes(fill = cluster), width = 1/12) +  # you might adjust width based on your preference
  geom_stratum(width = 1/12) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), min.segment.length = 0) +
  theme_minimal() +
  labs(title = "Changes in cluster membership across different gene counts",
       x = "Number of genes used in clustering",
       y = "Sample count",
       fill = "Cluster")  # to add legend title
# ----------------------------------------------------------------------------

# HIERARCHICAL CLUSTERING

# Calculate the distance matrix
data_to_cluster <- differential_expression_df[most_var_5000[,1], ]
dist_matrix <- dist(data_to_cluster, method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "complete")
plot(hclust_result, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_cluster_assignments <- cutree(hclust_result, k=10)

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








# ----------------------------------------------------------------------------

# CONSESUS CLUSTER PLUS








# ----------------------------------------------------------------------------

# GAUSSIAN MIXTURE MODELS

# installs and loads mclust package
install.packages("mclust")
library(mclust)


# performs GMM clustering
data_to_cluster <- differential_expression_df[most_var_5000[,1], ]
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
data_matrix <- as.matrix(differential_expression_df[most_var_5000[,1], ])
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

# Install and load the pheatmap package
install.packages("pheatmap")
library(pheatmap)


# different types, just label yours...


