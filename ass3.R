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
most_var_5000

# ----------------------------------------------------------------------------

# Clustering

# K MEANS
data_to_cluster <- differential_expression_df[most_var_5000[,1], ]

k <- 3
kmeans_result <- kmeans(data_to_cluster, centers = k)
kmeans_result


# cluster assignments for each sample
kmeans_cluster_assignments <- kmeans_result$cluster


table(kmeans_cluster_assignments)


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

# CONSESUS CLUSTER PLUS, RAMA
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")

data_to_cluster <- differential_expression_df[most_var_5000[,1], ]
data_matrix <- as.matrix(data_to_cluster)


#results for top 5000
library(ConsensusClusterPlus)
results_5000 <- ConsensusClusterPlus(data_matrix,
                                maxK=4, 
                                reps=1500, 
                                pItem=0.8, 
                                pFeature=1, 
                                clusterAlg="hc", 
                                distance="pearson",
                                seed=1262118388.71279)

#results for top 10
data_top_10 <- differential_expression_df[most_var_5000[1:10, 1], ]
matrix_top_10 <- as.matrix(data_top_10)
results_10 <- ConsensusClusterPlus(matrix_top_10,
                                     maxK=4, 
                                     reps=1500, 
                                     pItem=0.8, 
                                     pFeature=1, 
                                     clusterAlg="hc", 
                                     distance="pearson",
                                     seed=1262118388.71279)

#results for top 100
data_top_100 <- differential_expression_df[most_var_5000[1:100, 1], ]
matrix_top_100 <- as.matrix(data_top_100)
results_100 <- ConsensusClusterPlus(matrix_top_100,
                                   maxK=4, 
                                   reps=1500, 
                                   pItem=0.8, 
                                   pFeature=1, 
                                   clusterAlg="hc", 
                                   distance="pearson",
                                   seed=1262118388.71279)

#results for top 1000
data_top_1000 <- differential_expression_df[most_var_5000[1:1000, 1], ]
matrix_top_1000 <- as.matrix(data_top_1000)
results_1000 <- ConsensusClusterPlus(matrix_top_1000,
                                    maxK=4, 
                                    reps=1500, 
                                    pItem=0.8, 
                                    pFeature=1, 
                                    clusterAlg="hc", 
                                    distance="pearson",
                                    seed=1262118388.71279)


#results for top 10000
most_var_10000 <- sorted_gene_vars_no_na[1:10000, ]
data_top_10000 <- differential_expression_df[most_var_10000[1:10000, 1], ]
matrix_top_10000 <- as.matrix(data_top_10000)
results_10000 <- ConsensusClusterPlus(matrix_top_10000,
                                     maxK=4, 
                                     reps=1500, 
                                     pItem=0.8, 
                                     pFeature=1, 
                                     clusterAlg="hc", 
                                     distance="pearson",
                                     seed=1262118388.71279)


#alluvial plot
if (!require("ggalluvial", quietly = TRUE))
  install.packages("ggalluvial")

if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2")


library(ggalluvial)
library(reshape2)

cluster_assignments_10 <- results_10[[4]]$consensusClass
cluster_assignments_100 <- results_100[[4]]$consensusClass
cluster_assignments_1000 <- results_1000[[4]]$consensusClass
cluster_assignments_5000 <- results_5000[[4]]$consensusClass
cluster_assignments_10000 <- results_10000[[4]]$consensusClass
head(cluster_assignments_10)
all_clusters <- data.frame(sample_id = names(cluster_assignments_10000), 
                           most_var_10 = cluster_assignments_10,
                           most_var_100 = cluster_assignments_100,
                           most_var_1000 = cluster_assignments_1000,
                           most_var_5000 = cluster_assignments_5000,
                           most_var_10000 = cluster_assignments_10000
)
colnames(all_clusters)
head(all_clusters)
melted_data <- melt(all_clusters, id.vars = "sample_id", variable.name = "variation", value.name = "cluster")
colnames(melted_data)
library(ggalluvial)
ccp_alluvial <- ggplot(data = melted_data,
       aes(axis1 = sample_id, axis2 = variation, y = cluster)) +
  geom_alluvium(aes(fill = cluster)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()


ggsave(filename = "~/BioinformaticsProject/plots/ccp_alluvial.png", plot = ccp_alluvial, width = 10, height = 8, dpi = 400)


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


