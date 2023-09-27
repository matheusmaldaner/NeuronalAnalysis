library(dplyr)
library(tibble)
library(sessioninfo)
library(readr)
library(magrittr)
library(ggplot2)
library(org.Mm.eg.db)

# ====================================================
# step 1: Matheus Kunzler Maldaner
# ====================================================

# defines file paths
data_dir <- file.path("data", "SRP094496")
data_file <- file.path(data_dir, "SRP094496.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP094496.tsv")
results_dir <- "results"
plots_dir <- "plots"

# defines working directory
working_directory <- "C:/Users/Matheus/OneDrive/University of Florida/JUNIOR FALL/CGS4144/project"
setwd(working_directory)

metadata <- readr::read_tsv(metadata_file)


# loads the results file after running 'ensembl-to-hugo.R'
mapped_file <- file.path("results", "mapped_df.tsv")
mapped <- readr::read_tsv(mapped_file)


# drops the two non numeric columns
drop <- c("Ensembl","all_hugo_ids")
df_numerical <- mapped[,!(names(mapped) %in% drop)]


# drops all 14109 NAs
df_non_na <- na.omit(df_numerical) 

#df_averaged <- df_non_na %>%
#group_by(first_mapped_hugo_id) %>%
#  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = 'drop')


# keeps only unique hugo ids, averaging duplicate ones (below function is deprecated)
df_averaged <- df_non_na %>%
  group_by(first_mapped_hugo_id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop')


# turns first_mapped_hugo_id from a column to the row names
expression_df <- df_averaged %>%
  tibble::column_to_rownames("first_mapped_hugo_id")


# how many genes does it include (27080 genes and 1692 samples)
cat('It includes', dim(expression_df)[1], 'genes and', dim(expression_df)[2], 'samples.')


# how much variation do you see in the data?


# log scales the data without categorical columns
log_scaled_mapped <- log(expression_df[-1] + 1)


median_ranges <- apply(log_scaled_mapped, 1, function(row) {range(row)[2] - range(row)[1]})


density_plot <- ggplot(data.frame(median_ranges), aes(x=median_ranges)) + 
  geom_density(fill='blue') + 
  labs(title='Density Plot of Per-Gene Median Expression Ranges', 
       x='Median Expression Range', 
       y='Density')


plot(density_plot)

# -----------------------------------------------------------------------------

# ====================================================
# step 2: Ethan
# ====================================================

install.packages('umap')

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
library(DESeq2)
library(ggplot2)
library(umap)
library(cowplot)

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





expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

all.equal(colnames(expression_df), metadata$refinebio_accession_code)

rounded_expression_df <- round(expression_df)

dds <- DESeqDataSetFromMatrix(
  countData = rounded_expression_df,
  colData = metadata,
  design = ~title_status
)
# gc()
# rm(mapped_df)
# rm(expression_df)
# rm(df)
dds <- DESeq(dds)
vst_data <- vst(dds)
#Create pca plot
pcaData <- plotPCA(vst_data, intgroup = "title_status")
ggsave(
  plot = pcaData,
  file.path(plots_dir, "SRP123625_pca_plot.png")
)
# Extract the assay data
data_matrix <- assay(vst_data)
# Run UMAP
umap_results <- umap(t(data_matrix))

umap_df <- data.frame(UMAP1 = umap_results$layout[,1], 
                      UMAP2 = umap_results$layout[,2],
                      title_status = metadata$title_status)

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = title_status)) + 
  geom_point(size = 1, alpha = 0.8) + 
  theme_minimal() +
  ggtitle("UMAP plot") +
  theme(legend.position = "bottom")

plot(umap_plot)

ggsave(filename = file.path(plots_dir, "UMAP_plot.png"), plot = umap_plot, width = 6, height = 5)


# -----------------------------------------------------------------------------


# ====================================================
# step 3: Rama
# ====================================================


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
  lab = deseq_df$Ensembl,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
volcano_plot
ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "SRP123625_volcano_plot.png")
) # Replace with a plot name relevant to your data



# -----------------------------------------------------------------------------

# ====================================================
# step 4: Nathan
# ====================================================



# ====================================================
# step 5: Everyone
# ====================================================

# topGo




# ----------------------------------------------------
# clustProfiler




# ----------------------------------------------------
# gProfiler2




# ----------------------------------------------------
# GenomicSuperSignature






