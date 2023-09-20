# -----------------------------------------------------------------------------
# Script Name:        ensembl-to-hugo.R
# Authors:            Matheus, Rama, Nathan, Ethan 
# Created Date:       2023-09-20
# Last Modified Date: 2023-09-20
# Version:            1.0
# Purpose:            Provides conversion from Ensembl to Hugo IDs 
# -----------------------------------------------------------------------------


# =====================
# Instructions for Use:
# =====================

# 1. Set the working directory to where your data files are located.
working_directory <- "C:/Users/Matheus/OneDrive/University of Florida/JUNIOR FALL/CGS4144/project"

# 2. If your data files have different names, update the `data_file` and 
#    `metadata_file` variables and comment out lines [48] and [49].
data_file <- "your_data_file.tsv"
metadata_file <- "your_metadata_file.tsv"

# 3. If you want the results to be saved in a different directory, update the 
#    `results_dir` variable and comment out line [50].
results_dir <- "path/to/save/results"


# -----------------------------------------------------------------------------

setwd(working_directory)

# installs and load required packages
required_packages <- c("readr", "AnnotationDbi", "BiocManager", "sessioninfo", "magrittr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
library(sessioninfo)
library(readr)
library(magrittr)
# installs House Mouse Annotation Package if not already installed
if (!("org.Mm.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}
library(org.Mm.eg.db)


# defines file paths
data_dir <- file.path("data", "SRP094496")
data_file <- file.path(data_dir, "SRP094496.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP094496.tsv")
results_dir <- "results"


# imports and sets up data
metadata <- readr::read_tsv(metadata_file)
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")
expression_df <- expression_df %>% dplyr::select(metadata$refinebio_accession_code)
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")


# maps ENSEMBL to SYMBOL (Hugo)
mapped_list <- mapIds(
  org.Mm.eg.db,
  keys = expression_df$Gene,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "list")


# turns list into dataframe
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Hugo") %>%
  tidyr::unnest(cols = Hugo)


# counts occurrences of each Ensembl ID
multi_mapped <- mapped_df %>%
  dplyr::count(Ensembl, name = "hugo_id_count") %>%
  dplyr::arrange(desc(hugo_id_count))


# collapses multiple Hugo IDs into a single column
collapsed_mapped_df <- mapped_df %>%
  dplyr::group_by(Ensembl) %>%
  dplyr::summarize(all_hugo_ids = paste(Hugo, collapse = ";"))


# maps only the first instance of each Ensembl ID
final_mapped_df <- data.frame(
  "first_mapped_hugo_id" = mapIds(
    org.Mm.eg.db, 
    keys = expression_df$Gene,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
) %>%
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))


# writes results to file
readr::write_tsv(final_mapped_df, file.path(results_dir, "mapped_df"))

# prints session info
sessioninfo::session_info()

# -----------------------------------------------------------------------------

# Thank you!
# Project created for course CGS4144 at the University of Florida in Fall 2023
