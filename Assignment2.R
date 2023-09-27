# 1 - install and load dependencies to read tsv file
install.packages("readr")
install.packages("AnnotationDbi")
install.packages("BiocManager")
install.packages("sessioninfo")
library(sessioninfo)
library(readr)



# 1.1 - installing House Mouse Annotation Package
if (!("org.Mm.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}
library(org.Mm.eg.db) 
library(magrittr) # needed to use pipe: %>%

# defines the file path to directories
data_dir <- file.path("data", "SRP094496")
data_file <- file.path(data_dir, "SRP094496.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP094496.tsv")
results_dir <- "results"
mapped_file <- file.path(results_dir, "mapped_df.tsv")
plots_dir <- "plots"



# 1.2 Import and set up data

# reads metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# reads data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # tucks away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# reorders the data according to Metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# checks if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# prepares for Gene ID Mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")




# 1.3 maps ENSEMBL to SYMBOL (Hugo)
mapped_list <- mapIds(
  org.Mm.eg.db, # annotation package for mus_musculus
  keys = expression_df$Gene,
  keytype = "ENSEMBL",
  column = "SYMBOL", 
  multiVals = "list"
)

# turns list into dataframe for better usage
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Hugo") %>%
  tidyr::unnest(cols = Hugo)

# shows database preview and summary (14109 NA entries)
head(mapped_df)
summary(as.factor(mapped_df$Hugo), maxsum = 10)



# 1.4 counts # of times each Ensembl ID appears and arranges in descending order
# basically groups same Ensembl IDs - reduces dimensionality of matrix 
multi_mapped <- mapped_df %>%
  dplyr::count(Ensembl, name = "hugo_id_count") %>%
  dplyr::arrange(desc(hugo_id_count))

# previews the newly created object
head(multi_mapped)


multi_mapped



# 1.5 groups by Ensembl IDs and collapses them into one 'all_hugo_ids' column
collapsed_mapped_df <- mapped_df %>%
  dplyr::group_by(Ensembl) %>%
  dplyr::summarize(all_hugo_ids = paste(Hugo, collapse = ";"))

collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_entrez_ids` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_hugo_ids, ";")) %>%
  # We only need a preview here
  head()

collapsed_mapped_df



# 1.6 writes mapped results to file
readr::write_tsv(collapsed_mapped_df, file.path(
  results_dir,
  "SRP040561_Hugo_IDs.tsv" # can change this file name if needed
))

# prints session info
sessioninfo::session_info()


# end of step 1 


# reads in data for steps 2, 3, 4...
results_file <- file.path(results_dir, "SRP040561_Hugo_IDs.tsv")
df <- readr::read_tsv(results_file)


# ANYTHING BELOW THIS IS A WORK IN PROGRESS 
# FEEL FREE TO DELETE, CHANGE, ETC
head(df)

# calculates dimensions of matrix
matrix_size <- dim(df)
matrix_size

# num of genes in the dataframe
num_genes <- nrow(df)
num_genes

# log scaled data
log_scaled_data <- log(df + 1)

per_gene_median_range <- apply(log_scaled_data, 1, function(x) max(x) - min(x))

library(ggplot2)

ggplot(data.frame(PerGeneMedianRange = per_gene_median_range), aes(x = PerGeneMedianRange)) +
  geom_density(fill = "blue", alpha = 0.5) +
  xlab("Per-Gene Median Expression Range (log-scaled)") +
  ylab("Density") +
  ggtitle("Density Plot of Per-Gene Median Expression Ranges")






# < IN PROGRESS >
# R Code for Expression Matrix and PCA Plot

# Log-scale the data
#>log_expression_data <- log(expression_data, 2)

# Calculate per-gene median expression ranges and make a density plot
# ...

# Generate PCA plot using DESeq2
# ...

# Color the plot by groups (e.g., cancer vs normal)
# ...

# Generate t-SNE or UMAP plot if applicable
# ...




# R Code for Differential Analysis
# Perform differential analysis

# Install these packages if they aren't installed yet

if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

#attaching libraries
library(DESeq2)
library(ggplot2)
library(magrittr)

#set random seed
set.seed(50341)


metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(mapped_file) %>%
  tibble::column_to_rownames("Ensembl")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

expression_df <- expression_df %>%
  dplyr::select(first_mapped_hugo_id, all_hugo_ids, metadata$refinebio_accession_code)
rownames(expression_df)


# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

#since for out two groups were using titles pm (85 unique, 211 instances) and all others (~1420)
head(metadata$refinebio_title)

# Set up metadata
metadata <- metadata %>%
  dplyr::mutate(title_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "PM-\\d+") ~ "PM",
    TRUE ~ "other"
  ))
table(metadata$title_status)
any(is.na(metadata$title_status))


dplyr::select(metadata, refinebio_title, title_status)
str(metadata$title_status)


metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    title_status = factor(title_status, levels = c("PM", "other"))
  )
any(is.na(metadata$title_status))

table(metadata$title_status)

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)


# Create a DESeq2Dataset
# round all expression counts
gene_matrix <- round(filtered_expression_df)

any(is.na(metadata$title_status))
table(metadata$title_status)



ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~title_status
)

deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

head(deseq_results)

deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Ensembl") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

deseq_df$symbol <- mapIds(org.Mm.eg.db, 
                          keys = deseq_df$Ensembl,
                          keytype = "ENSEMBL",
                          column = "SYMBOL",
                          multiVals = "first")

head(deseq_df)

#plotCounts(ddset, gene = "ENSMUSG00000021743", intgroup = "title_status")

readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "SRP094496_diff_expr_results.tsv"
  )
)

# Create a volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$symbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
volcano_plot
ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "SRP123625_volcano_plot.png")
) # Replace with a plot name relevant to your data
 

# Create a table of differentially expressed genes
install.packages("devtools")
library(devtools)
install_github("jokergoo/ComplexHeatmap", force = TRUE)
library(ComplexHeatmap)

if (!("org.Mm.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}
library(org.Mm.eg.db)

diff_expr_df <- readr::read_tsv("results/SRP094496_diff_expr_results.tsv")

significant_df <- subset(diff_expr_df, (threshold == T))
diff_expr_df <- subset(diff_expr_df, ((diff_expr_df$baseMean >= 7) & abs(diff_expr_df$log2FoldChange) >= 1) & (threshold == T))
diff_expr_df

diff_expr_df$symbol <- mapIds(org.Mm.eg.db, 
                            keys = diff_expr_df$Ensembl,
                            keytype = "ENSEMBL",
                            column = "SYMBOL",
                            multiVals = "first")


# video says to use normalized = T in the counts() func
mat <- counts(ddset)[diff_expr_df$Ensembl,]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(metadata)

hm_plot <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z),
                   row_labels = diff_expr_df$symbol, name = "Z-score")
hm_plot


# R Code for Summary Tables and Plots
# # Create tables of enriched processes


# TopGo
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")
library(topGO)
library(org.Mm.eg.db)

deseq_df <- file.path(results_dir, "SRP094496_diff_expr_results.tsv")
deseq_data <- readr::read_tsv(deseq_df)

#named vector of p-values
all_genes <- setNames(deseq_data$pvalue, deseq_data$Ensembl)
#statistically sig genes
geneSelectionFunc <- function(pvalues) {
  return(pvalues < 0.01)
}



selected_genes <- sum(geneSelectionFunc(all_genes))
print(selected_genes)

library(org.Mm.eg.db)

# Pick a random gene from your list
keytypes(org.Mm.eg.db)
head(names(all_genes))

valid_keys <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
head(intersect(valid_keys, names(all_genes)))
length(valid_keys)

invalid_keys <- setdiff(names(all_genes), valid_keys)
head(invalid_keys)
length(invalid_keys)

all_genes_valid <- all_genes[names(all_genes) %in% valid_keys]


sample_gene <- sample(names(all_genes), 1)
select(org.Mm.eg.db, keys = sample_gene, columns = "GO", keytype = "ENSEMBL")


# Create topGOdata Object:
GOdata <- new("topGOdata",
              description = "My GO Analysis",
              ontology = "BP", 
              allGenes = all_genes, 
              geneSel = geneSelectionFunc, 
              nodeSize = 10, 
              annot = annFUN.org,
              mapping = "org.Mm.eg.db",
              ID = "ENSEMBL")

GOdata

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

tableFisher <- GenTable(GOdata, 
                        classicFisher = resultFisher, 
                        topNodes = 10, 
                        orderBy = "classicFisher", 
                        ranksOf = "classicFisher")
print(tableFisher)


allRes <- GenTable(GOdata, classicFisher = resultFisher,topNodes = 10)

# Write summaries for each plot/table
# ...
