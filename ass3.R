# CHANGE PATH VARIABLE
setwd("C:/Users/twinb/OneDrive/Desktop/school/uf_fall_23/cgs4144/BioinformaticsProject/")

library(readr)
library(dplyr)
library(matrixStats)

## Grab 5000 most variable genes -----------------------------
# load diff expr results file
expression_df <- readr::read_tsv("results/mapped_df.tsv")
metadata_df <- readr::read_tsv("data/SRP094496/metadata_SRP094496.tsv")
expression_df

gene_expression <- expression_df

# store gene IDs in vars and drop them from main df
# this is so rowVars works (doesn't work with columsn of val <char>)
symbol <- gene_expression$first_mapped_hugo_id
Ensembl <- expression_df$Ensembl
gene_expression$first_mapped_hugo_id <- NULL
gene_expression$all_hugo_ids <- NULL
gene_expression$Ensembl <- NULL

# calculate variances for each gene
gene_variability <- rowVars(as.matrix(gene_expression))

# bind the variances and gene IDs
unsorted_gene_vars <- cbind(1:length(symbol), symbol, gene_variability)
unsorted_gene_vars
sorted_gene_vars <- unsorted_gene_vars[order(as.numeric(unsorted_gene_vars[, 3]), decreasing=TRUE), ] 
sorted_gene_vars

# drop NA vals
sorted_gene_vars_no_na <- sorted_gene_vars[complete.cases(sorted_gene_vars), ]
sorted_gene_vars_no_na

# grab first 5000
most_var_5000 <- sorted_gene_vars_no_na[1:5000, ]
most_var_10000 <- sorted_gene_vars_no_na[1:10000, ]
most_var_5000

data_to_cluster <- gene_expression[most_var_5000[,1], ]
# ----------------------------------------------------------------------------

# Clustering
install.packages("ggalluvial")
library(cluster)
library(ggplot2)
library(ggalluvial)

# K MEANS
data_to_cluster <- gene_expression[most_var_5000[,1], ]
data_to_cluster_10000 <- gene_expression[most_var_10000[,1], ]
#Set to 3 clusters
k <- 5
kmeans_result <- kmeans(data_to_cluster, centers = k)
kmeans_result

# cluster assignments for each sample
kmeans_cluster_assignments <- kmeans_result$cluster

# Convert to data frame for easier CSV writing
table_df <- as.data.frame(table(kmeans_cluster_assignments))

# Write to CSV
write.csv(table_df, file = "~/BioinformaticsProject/plots/tables/kmeans_cluster_table.csv", row.names = FALSE)

cluster_results <- cbind(1:5000, kmeans_cluster_assignments)

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
png("~/BioinformaticsProject/plots/k_means_silhouette.png")

plot(k_values, sil_widths, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters 'k'", ylab = "Average silhouette width",
     main = "Silhouette Analysis of k-means clustering")

dev.off()

#Test different number of genes
k <- 5
gene_values <- c(10, 100, 1000)
kmeans_genes_results <- list()
for (g in gene_values) {
  set.seed(123) # Set seed for reproducibility
  kmeans_genes_results[[paste(g, "genes", sep = "_")]] <- kmeans(data_to_cluster[1:g, ], centers = k)
}
set.seed(123) # Set seed for reproducibility
kmeans_genes_results[[paste("10000", "genes", sep = "_")]] <- kmeans(data_to_cluster_10000, centers = k)

cluster_data <- lapply(kmeans_genes_results, function(x) x$cluster)
names(cluster_data) <- c("kmeans_10", "kmeans_100", "kmeans_1000", "kmeans_10000")  # renaming for clarity
cluster_df <- as.data.frame(cluster_data)

library(tidyr)

# Reshaping data for plotting
cluster_df_long <- cluster_df %>% 
  pivot_longer(cols = starts_with("kmeans"),
               names_to = "gene_count",
               values_to = "cluster")

# Converting cluster to factor
cluster_df$kmeans_10 <- as.factor(cluster_df$kmeans_10)
cluster_df$kmeans_100 <- as.factor(cluster_df$kmeans_100)
cluster_df$kmeans_1000 <- as.factor(cluster_df$kmeans_1000)
cluster_df$kmeans_10000 <- as.factor(cluster_df$kmeans_10000)

kmeans_alluvial <- ggplot(data = cluster_df,
                          aes(axis1 = kmeans_10, axis2 = kmeans_100, 
                              axis3 = kmeans_1000, axis4 = kmeans_10000)) +
  geom_alluvium(aes(fill = kmeans_10000), width = 1/12) +  # Fill based on final clustering, adjust width as needed
  geom_stratum(width = 1/12) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), min.segment.length = 0) +
  theme_minimal() +
  labs(title = "Changes in K-means Cluster Membership Across Different Gene Counts",
       x = "Number of Genes Used in Clustering",
       y = "Gene Count",
       fill = "Cluster")  # To add legend title

print(kmeans_alluvial)

ggsave(filename = "~/BioinformaticsProject/plots/k_means_alluvial.png", plot = kmeans_alluvial, width = 10, height = 8, dpi = 400)

# ----------------------------------------------------------------------------

# HIERARCHICAL CLUSTERING

# Calculate the distance matrix
data_to_cluster <- gene_expression[most_var_5000[,1], ]
dist_matrix <- dist(data_to_cluster, method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "complete")
plot(hclust_result, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_cluster_assignments <- cutree(hclust_result, k=50)

cluster_results <- cbind(1:5000, hclust_cluster_assignments)
table(hclust_cluster_assignments)

# test with 10 genes
hc <- hclust(dist(data_to_cluster[1:10, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_10 <- cutree(hc, h=3000)
table(hclust_10)
cluster_results <- cbind(cluster_results, hclust_10)

# test with 100 genes
hc <- hclust(dist(data_to_cluster[1:100, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_100 <- cutree(hc, h=3000)
table(hclust_100)
cluster_results <- cbind(cluster_results, hclust_100)

# test w 1000 genes
hc <- hclust(dist(data_to_cluster[1:1000, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_1000 <- cutree(hc, h=3000)
table(hclust_1000)
cluster_results <- cbind(cluster_results, hclust_1000)

# test w 10000 genes
hc <- hclust(dist(data_to_cluster[1:10000, ],method="euclidean"),method="complete")
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "Samples")
hclust_10000 <- cutree(hc, k=50)
table(hclust_10000)
cluster_results <- cbind(cluster_results, hclust_10000)

library(ggplot2)
library(ggalluvial)

## sankey/alluvial plot
cluster_results <- as.data.frame(cluster_results)
cluster_results
hclust_alluvial <- ggplot(data = cluster_results,
                                aes(axis1 = hclust_10, axis2 = hclust_100, axis3=hclust_1000, axis4=hclust_10000)) +
  geom_alluvium(aes(fill = hclust_cluster_assignments), width = 1/12) +  # you might adjust width based on your preference
  geom_stratum(width = 1/12) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Changes in cluster membership across different gene counts",
       x = "Number of genes used in clustering",
       y = "Sample count",
       fill = "Cluster")  # to add legend title

hclust_alluvial

# ----------------------------------------------------------------------------

# CONSESUS CLUSTER PLUS, RAMA
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

data_to_cluster <- gene_expression[most_var_5000[,1], ]
data_matrix = as.dist(1 - cor(t(data_to_cluster), method="pearson"))
#results for top 5000
library(ConsensusClusterPlus)
results_5000 <- ConsensusClusterPlus(data_matrix,
                                maxK=3, 
                                reps=50, 
                                pItem=0.8, 
                                pFeature=1, 
                                clusterAlg="hc", 
                                distance="pearson",
                                seed=12345)


#results for top 10
data_top_10 <- gene_expression[most_var_5000[1:10, 1], ]
matrix_top_10 = as.dist(1 - cor(t(data_top_10), method="pearson"))
results_10 <- ConsensusClusterPlus(matrix_top_10,
                                     maxK=3,
                                     reps=50,
                                     pItem=0.8,
                                     pFeature=1,
                                     clusterAlg="hc",
                                     distance="pearson",
                                     seed=12345)

#results for top 100
data_top_100 <- gene_expression[most_var_5000[1:100, 1], ]
matrix_top_100 = as.dist(1 - cor(t(data_top_100), method="pearson"))
results_100 <- ConsensusClusterPlus(matrix_top_100,
                                   maxK=3, 
                                   reps=50, 
                                   pItem=0.8, 
                                   pFeature=1, 
                                   clusterAlg="hc", 
                                   distance="pearson",
                                   seed=12345)

#results for top 1000
data_top_1000 <- gene_expression[most_var_5000[1:1000, 1], ]
matrix_top_1000 = as.dist(1 - cor(t(data_top_1000), method="pearson"))
results_1000 <- ConsensusClusterPlus(matrix_top_1000,
                                    maxK=3, 
                                    reps=50, 
                                    pItem=0.8, 
                                    pFeature=1, 
                                    clusterAlg="hc", 
                                    distance="pearson",
                                    seed=12345)


#results for top 10000
most_var_10000 <- sorted_gene_vars_no_na[1:10000, ]
data_top_10000 <- gene_expression[most_var_10000[1:10000, 1], ]
matrix_top_10000 = as.dist(1 - cor(t(data_top_10000), method="pearson"))
results_10000 <- ConsensusClusterPlus(matrix_top_10000,
                                     maxK=3, 
                                     reps=1, 
                                     pItem=0.5, 
                                     pFeature=1, 
                                     clusterAlg="hc", 
                                     distance="pearson",
                                     seed=12345)


#alluvial plot
if (!require("ggalluvial", quietly = TRUE))
  install.packages("ggalluvial")

if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2")


library(ggalluvial)
library(reshape2)


cluster_assignments_10 <- results_10[[3]]$consensusClass ##SHOWS SAMPLES IN CLUSTERING
cluster_assignments_100 <- results_100[[3]]$consensusClass
cluster_assignments_1000 <- results_1000[[3]]$consensusClass
cluster_assignments_5000 <- results_5000[[3]]$consensusClass
cluster_assignments_10000 <- results_10000[[3]]$consensusClass


data <- data.frame(
  id = 1:length(cluster_assignments_10),
  assignment_10 = cluster_assignments_10,
  assignment_100 = cluster_assignments_100,
  assignment_1000 = cluster_assignments_1000,
  assignment_5000 = cluster_assignments_5000,
  assignment_10000 = cluster_assignments_10000
)

# 2. Plot with ggalluvial
library(ggplot2)
library(ggalluvial)

ccp_alluvial <- ggplot(data = data, 
       aes(axis1 = assignment_10, axis2 = assignment_100, axis3 = assignment_1000, axis4 = assignment_5000, axis5 = assignment_10000)) +
  geom_alluvium(aes(fill = assignment_5000), width = 1/12) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  labs(title = "Changes in cluster membership from 10 to 10,000",
       x = "Number of genes used in clustering",
       y = "gene count",
       fill = "Cluster")

print(ccp_alluvial)

ggsave(filename = "~/BioinformaticsProject/plots/ccp_alluvial.png", plot = ccp_alluvial, width = 10, height = 8, dpi = 400)


# ----------------------------------------------------------------------------

# GAUSSIAN MIXTURE MODELS

# installs and loads mclust package
install.packages("mclust")
library(mclust)


# performs GMM clustering
data_to_cluster <- gene_expression[most_var_5000[,1], ]
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

# install.packages(c("ggplot2", "FactoMineR"))
library(ggplot2)
library(FactoMineR)

# performs PCA
data_matrix <- as.matrix(gene_expression[most_var_5000[,1], ])
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

# Create a heatmap
data_to_cluster <- gene_expression[most_var_5000[,1], ]
heatmap_data <- data_to_cluster

# Create row and column dendrograms
row_dend <- as.dendrogram(hclust(dist(heatmap_data)))
col_dend <- as.dendrogram(hclust(dist(t(heatmap_data))))

# Create annotation data
annotation_data <- data.frame(
  hclust = hclust_cluster_assignments,
  kmeans = kmeans_cluster_assignments
)
colnames(heatmap_data)

# Create the heatmap
pheatmap(
  as.matrix(heatmap_data,  rownames.force = TRUE),
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  show_colnames = FALSE,  # You can set this to TRUE if you want to display column names
  legend = TRUE,
  annotation_row = annotation_data, 
  main = "Heatmap of 5,000 Genes",
  filename = "heatmap6.png"  # You can specify the file name and format
)

# ----------------------------------------------------------------------------

# Statistics 

# Chi-squared

# CCP
