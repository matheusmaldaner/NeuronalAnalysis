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

## Clustering -----------------------------

# testing clustering w K-Means
k <- 3
kmeans_result <- kmeans(differential_expression_df[most_var_5000[,1], ], centers = k)
kmeans_result

table(kmeans_result$cluster)
