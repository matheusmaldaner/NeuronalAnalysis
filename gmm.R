# CHANGE PATH VARIABLE
setwd("C:/Users/Matheus/OneDrive/University of Florida/JUNIOR FALL/CGS4144/project")

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



# -----------------------------------------------------------------------------------------------





# GAUSSIAN MIXTURE MODELS

most_var_1000 <- sorted_gene_vars_no_na[1:1000, ]


library(mclust)


# defines gene values
gene_values <- c(10, 100, 1000, 10000)


# performs GMM clustering
data_to_cluster <- gene_expression[most_var_1000[,1], ]
gmm_result <- Mclust(data_to_cluster)


# prints the clustering results
summary(gmm_result)


# number of clusters identified
num_clusters <- gmm_result$G
cat("Number of clusters identified by GMM:", num_clusters, "\n")


# cluster assignments for each sample
gmm_cluster_assignments <- gmm_result$classification
table(gmm_cluster_assignments)


# Rerun the clustering method using different numbers of genes
gmm_genes_results <- list()
for (g in gene_values) {
  gmm_genes_results[[paste(g, "genes", sep = "_")]] <- Mclust(data_to_cluster[1:g, ])
}


# Print the number of clusters for each gene count
lapply(gmm_genes_results, function(result) {
  num_clusters <- result$G
  return(num_clusters)
})


# Data Transformation for Visualization
cluster_data_gmm <- lapply(gmm_genes_results, function(x) x$classification)
names(cluster_data_gmm) <- c("gmm_10", "gmm_100", "gmm_1000")  # renaming for clarity
cluster_df_gmm <- as.data.frame(cluster_data_gmm)


library(tidyr)
library(ggplot2)
library(ggalluvial)


# Reshaping data for plotting
cluster_df_gmm_long <- cluster_df_gmm %>% 
  pivot_longer(cols = starts_with("gmm"),
               names_to = "gene_count",
               values_to = "cluster")


# Converting cluster to factor for each gene count
cluster_df_gmm$gmm_10 <- as.factor(cluster_df_gmm$gmm_10)
cluster_df_gmm$gmm_100 <- as.factor(cluster_df_gmm$gmm_100)
cluster_df_gmm$gmm_1000 <- as.factor(cluster_df_gmm$gmm_1000)


# Alluvial Plot Visualization
gmm_alluvial <- ggplot(data = cluster_df_gmm,
                       aes(axis1 = gmm_10, axis2 = gmm_100, axis3 = gmm_1000)) +
  geom_alluvium(aes(fill = gmm_1000), width = 1/12) +
  geom_stratum(width = 1/12) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), min.segment.length = 0) +
  theme_minimal() +
  labs(title = "Changes in GMM Cluster Membership Across Different Gene Counts",
       x = "Number of Genes Used in Clustering",
       y = "Gene Count",
       fill = "Cluster")


print(gmm_alluvial)


ggsave(filename = "GMM_alluvial.png", plot = gmm_alluvial, width = 10, height = 8, dpi = 400)


