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
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

#since for out two groups were using titles pm (85) and all others (~1520)
head(metadata$refinebio_title)

# Create a volcano plot
# ...

# Create a table of differentially expressed genes
# ...





# R Code for Summary Tables and Plots
# # Create tables of enriched processes
# ...

# Write summaries for each plot/table
# ...
