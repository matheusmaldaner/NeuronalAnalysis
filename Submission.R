library(dplyr)
library(tibble)
library(sessioninfo)
library(readr)
library(magrittr)
library(ggplot2)
library(org.Mm.eg.db)
library(DESeq2)
library(topGO)
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


# keeps only unique hugo ids, averaging duplicate ones (below function is deprecated)
df_averaged <- df_non_na %>%
group_by(first_mapped_hugo_id) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = 'drop')




# turns first_mapped_hugo_id from a column to the row names
expression_df <- df_averaged %>%
  tibble::column_to_rownames("first_mapped_hugo_id")


# how many genes does it include (27080 genes and 1692 samples)
cat('It includes', dim(expression_df)[1], 'genes and', dim(expression_df)[2], 'samples.')


# how much variation do you see in the data?


# log scales the data without categorical columns
log_scaled_mapped <- log(expression_df[-1] + 1)


median_ranges <- apply(log_scaled_mapped, 1, function(row) {range(row)[2] - range(row)[1]})

# creates density plot 
density_plot <- ggplot(data.frame(median_ranges), aes(x = median_ranges)) + 
  geom_density(fill = 'blue') + 
  labs(title = 'Density Plot of Per-Gene Median Expression Ranges', 
       x = 'Median Expression Range', y = 'Density') +
  theme_minimal(base_size = 15)

# saves density plot png
ggsave(filename = file.path(plots_dir, "density_plot.png"), 
       plot = density_plot, 
       width = 7, 
       height = 7, 
       bg = "white")

# ====================================================
# step 2: Ethan Fan
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


dds <- DESeq(dds)
vst_data <- vst(dds)


# Extracting PCA data
pca_data <- plotPCA(vst_data, intgroup = "title_status", returnData = TRUE)

# Extracting principal components
pca_object <- DESeq2::plotPCA(vst_data, intgroup = "title_status", returnData = TRUE)

# Extracting the percent variance explained by each principal component.
percentVar <- round(100 * attr(pca_object, "percentVar"))


# Creating ggplot object
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = title_status)) +
  geom_point(size = 0.5, alpha = 0.8) +
  theme_minimal(base_size = 15) +
  ggtitle("PCA Plot") +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance"), 
       color = "Title Status") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("red", "blue"))  # adjust colors according to your factor levels


# Saving the plot
ggsave(filename = file.path(plots_dir, "SRP123625_pca_plot.png"), plot = pca_plot, width = 7, height = 7, bg = "white")


# extracts the assay data
data_matrix <- assay(vst_data)


# runs UMAP
umap_results <- umap(t(data_matrix))

# creates UMAP dataframe
umap_df <- data.frame(UMAP1 = umap_results$layout[,1], 
                      UMAP2 = umap_results$layout[,2],
                      title_status = metadata$title_status)


# creates UMAP plot
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = title_status)) + 
  geom_point(size = 0.5, alpha = 0.8) + 
  theme_minimal(base_size = 15) +
  ggtitle("UMAP Plot") +
  labs(color = "Labels") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("red", "blue"))


# saves UMAP plot png
ggsave(filename = file.path(plots_dir, "UMAP_plot.png"),
       plot = umap_plot,
       width = 7,
       height = 7,
       bg = "white")


# ====================================================
# step 3: Rama Janco
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

#set random seed
set.seed(50341)


# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)


# Create a DESeq2Dataset
# round all expression counts
gene_matrix <- round(filtered_expression_df)

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
  tibble::rownames_to_column("Hugo") %>%
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
  lab = deseq_df$Hugo,
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

install.packages("devtools")
library(devtools)
install_github("jokergoo/ComplexHeatmap", force = TRUE)
library(ComplexHeatmap)

if (!("org.Mm.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}
library(org.Mm.eg.db)

diff_expr_df <- readr::read_tsv("results/SRP094496_diff_expr_results.tsv")

# filter differential expression data to make a manageable subset for the heatmap
significant_df <- subset(diff_expr_df, (threshold == T))
significant_df <- subset(significant_df, ((significant_df$baseMean >= 7) & abs(significant_df$log2FoldChange) >= 1))
significant_df

# calculate z-scores for the genes provided
mat <- counts(ddset)[significant_df$symbol,]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(metadata)

# create heatmap given z-score data
hm_plot <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z),
                   row_labels = diff_expr_df$symbol, name = "Z-score")
hm_plot


# ====================================================
# step 5: Everyone
# ====================================================

# topGo, BP
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")


deseq_df <- file.path(results_dir, "SRP094496_diff_expr_results.tsv")
deseq_data <- readr::read_tsv(deseq_df)

# named vector of p-values
all_genes <- setNames(deseq_data$pvalue, deseq_data$Hugo)

# statistically sig genes
geneSelectionFunc <- function(pvalues) {
  return(pvalues < 0.01)
}

selected_genes <- sum(geneSelectionFunc(all_genes))
print(selected_genes)


# Pick a random gene from your list
keytypes(org.Mm.eg.db)
head(names(all_genes))

valid_keys <- keys(org.Mm.eg.db, keytype = "SYMBOL")
head(intersect(valid_keys, names(all_genes)))
length(valid_keys)

invalid_keys <- setdiff(names(all_genes), valid_keys)
head(invalid_keys)
length(invalid_keys)

all_genes_valid <- all_genes[names(all_genes) %in% valid_keys]

# Create topGOdata Object:
GOdata <- new("topGOdata",
              description = "My GO Analysis",
              ontology = "BP", 
              allGenes = all_genes, 
              geneSel = geneSelectionFunc, 
              nodeSize = 10, 
              annot = annFUN.org,
              mapping = "org.Mm.eg.db",
              ID = "SYMBOL")

GOdata

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
## other enrichment tests
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

tableKS <- GenTable(GOdata, 
                    classicKS = resultKS, 
                    topNodes = 10, 
                    orderBy = "classicKS", 
                    ranksOf = "classicKS")
print(tableKS)

tableKSelim <- GenTable(GOdata, 
                        elimKS = resultKS.elim, 
                        topNodes = 10, 
                        orderBy = "elimKS", 
                        ranksOf = "elimKS")
print(tableKSelim)
## other^
tableFisher <- GenTable(GOdata, 
                        classicFisher = resultFisher, 
                        topNodes = 15, 
                        orderBy = "classicFisher", 
                        ranksOf = "classicFisher")
print(tableFisher)


readr::write_tsv(
  tableFisher,
  file.path(plots_dir, "topGO_table.tsv")
)


allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

readr::write_tsv(
  allRes,
  file.path(tables_dir, "allRes.tsv")
)

print(allRes)


# ----------------------------------------------------
# clustProfiler




# ----------------------------------------------------
# gProfiler2




# ----------------------------------------------------
# topGO, ontology = MF
# Create topGOdata Object:
mfGOdata <- new("topGOdata",
              description = "My MF GO Analysis",
              ontology = "MF", 
              allGenes = all_genes, 
              geneSel = geneSelectionFunc, 
              nodeSize = 10, 
              annot = annFUN.org,
              mapping = "org.Mm.eg.db",
              ID = "SYMBOL")

mfGOdata

resultFisher <- runTest(mfGOdata, algorithm = "classic", statistic = "fisher")
## other enrichment tests
resultKS <- runTest(mfGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(mfGOdata, algorithm = "elim", statistic = "ks")

tableKS <- GenTable(mfGOdata, 
                    classicKS = resultKS, 
                    topNodes = 10, 
                    orderBy = "classicKS", 
                    ranksOf = "classicKS")
print(tableKS)

tableKSelim <- GenTable(mfGOdata, 
                        elimKS = resultKS.elim, 
                        topNodes = 10, 
                        orderBy = "elimKS", 
                        ranksOf = "elimKS")
print(tableKSelim)
## other^
tableFisher <- GenTable(mfGOdata, 
                        classicFisher = resultFisher, 
                        topNodes = 15, 
                        orderBy = "classicFisher", 
                        ranksOf = "classicFisher")
print(tableFisher)


readr::write_tsv(
  tableFisher,
  file.path(plots_dir, "MFtopGO_table.tsv")
)


allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

readr::write_tsv(
  allRes,
  file.path(tables_dir, "MFallRes.tsv")
)

print(allRes)





