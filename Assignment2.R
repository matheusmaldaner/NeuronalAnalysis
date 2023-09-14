# install and load dependencies to read tsv file
install.packages("readr")
install.packages("AnnotationDbi")
install.packages("BiocManager")
install.packages("sessioninfo")
library(sessioninfo)
library(readr)


# installing House Mouse Annotation Package
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
plots_dir <- "plots"


# (4.2) Import and set up data

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



# (4.3) Map Ensembl IDs to Entrez IDs
mapped_list <- mapIds(
  org.Mm.eg.db, # annotation package for mus_musculus
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "ENTREZID", # The type of gene identifiers you would like to map to
  multiVals = "list"
)


# 4.4 Explore gene ID conversion
# Let's use the `head()` function for a preview of our mapped list
head(mapped_list)


# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Entrez") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Entrez)


head(mapped_df)


# Use the `summary()` function to show the distribution of Entrez values
# We need to use `as.factor()` here to get the count of unique values
# `maxsum = 10` limits the summary to 10 distinct values
summary(as.factor(mapped_df$Entrez), maxsum = 10)



multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "entrez_id_count") %>%
  # Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(entrez_id_count))

# Let's look at the first 6 rows of our `multi_mapped` object
head(multi_mapped)




# 4.4.1 Store all mapped values
collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
  # Collapse the Entrez IDs `mapped_df` into one column named `all_entrez_ids`
  dplyr::summarize(all_entrez_ids = paste(Entrez, collapse = ";"))



collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_entrez_ids` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_entrez_ids, ";")) %>%
  # We only need a preview here
  head()



# 4.4.2 Map Ensembl IDs to Entrez – keeping only the first mapped value
final_mapped_df <- data.frame(
  "first_mapped_entrez_id" = mapIds(
    org.Mm.eg.db, # mus_musculus
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the gene identifiers used in your data
    column = "ENTREZID", # The type of gene identifiers you would like to map to
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))

final_mapped_df %>%
  # Filter `final_mapped_df` to rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_entrez_ids, ";")) %>%
  head()


# 4.5 Write mapped results to file
# Write mapped and annotated data frame to output file
readr::write_tsv(final_mapped_df, file.path(
  results_dir,
  "SRP040561_Entrez_IDs.tsv" # Replace with a relevant output file name
))



# Print session info
sessioninfo::session_info()






# R Code for Expression Matrix and PCA Plot

# Load the data into R
# Assuming the data is in a CSV file called 'expression_data.csv'

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
# ...

# Create a volcano plot
# ...

# Create a table of differentially expressed genes
# ...





# R Code for Summary Tables and Plots
# # Create tables of enriched processes
# ...

# Write summaries for each plot/table
# ...
