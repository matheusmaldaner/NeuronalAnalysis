library(readr)
library(dplyr)
library(matrixStats)

diff_expr_df <- readr::read_tsv("results/SRP094496_diff_expr_results.tsv")
head(diff_expr_df)

genes <- diff_expr_df$symbol
diff_expr_df$symbol <- NULL
diff_expr_df$Ensembl <- NULL


# Assuming your data frame is named 'differential_expression_df'
# Calculate the variance for each gene
gene_variability <- rowVars(as.matrix(diff_expr_df))
gene_variability
# Sort genes based on variability in descending order
varsID <- cbind(seq(1, length(gene_variability), 1), gene_variability)
genesID <- cbind(seq(1, length(gene_variability), 1), genes)
varsID
genesID

?names

result
genes[1]
seq(1, length(gene_variability), 1)
dim(result)
result[,2]

sorted_genes <- sort(varsID, decreasing = TRUE)
sorted_genes

# Select the top 5000 most variable genes
top_5000_variable_genes <- sorted_genes[1:5000]

head(top_5000_variable_genes)
top
