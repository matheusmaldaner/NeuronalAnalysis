ibrary(e1071)
library(tidymodels)  # tidymodels includes dplyr and broom, which are used for data manipulation and modeling
library(caret)       # caret is used for creating consistent interfaces for model training

# Read the metadata
metadata <- readr::read_tsv(metadata_file)

# Data processing, assuming this has already been done as per the previous steps you have
metadata <- metadata %>%
  dplyr::mutate(title_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "PM-\\d+") ~ 1,
    TRUE ~ 0
  ))
labels <- metadata$title_status

# Data setup for clustering and model training
df <- data.frame(t(data_to_cluster), title_status = labels)
df <- df[sample(nrow(df)), ]  # Shuffle the rows
split_indices <- createDataPartition(df$title_status, p = 0.2, list = FALSE)
train_data <- df[-split_indices, ]
test_data <- df[split_indices, ]

# trains a Naive Bayes model
naive_bayes_model <- naiveBayes(title_status ~ ., data = train_data)
predictions <- predict(naive_bayes_model, newdata = test_data, type = "raw")

df_predictions <- data.frame(predictions = ifelse(predictions[,2] > 0.5, 1, 0))
table(df_predictions$predictions)

# AUC calculation for different gene sets
auc_list <- list()
data_to_cluster <- gene_expression[most_var_10000[,1], ]
for (x in c(10,100,1000,10000)) {
  df <- data.frame(t(data_to_cluster[1:x,]), title_status = labels)
  df <- df[sample(nrow(df)), ]
  
  split_indices <- createDataPartition(df$title_status, p = 0.2, list = FALSE)
  train_data <- df[-split_indices, ]
  test_data <- df[split_indices, ]
  
  naive_bayes_model <- naiveBayes(title_status ~ ., data = train_data)
  predictions <- predict(naive_bayes_model, newdata = test_data, type = "raw")
  
  df_predictions <- data.frame(predictions = ifelse(predictions[,2] > 0.5, 1, 0))
  print(table(df_predictions$predictions))
  
  # Use pROC library for AUC
  roc_obj <- roc(test_data$title_status, predictions[,2])
  auc_value <- auc(roc_obj)
  auc_list <- append(auc_list, auc_value)
  
}

# print this in a fancy way like Nathan did
auc_list


# Heatmap plotting, assuming that 'test_data' and 'predictions' are defined as before
library(pheatmap)
library(dendextend)

# only uses positive predictions
positive_predictions <- predictions[,2] > 0.5
scaled <- log(test_data[positive_predictions,] + 1)
annotation_data_heatmap <- data.frame(status = round(scaled$title_status))

# removes nulls
scaled[is.na(scaled)] <- 0

# plots heatmap
pheatmap(
  as.matrix(scaled,  rownames.force = NA),
  cluster_rows = TRUE,
  show_colnames = FALSE,  # You can set this to TRUE if you want to display column names
  legend = TRUE, 
  main = "Heatmap of Genes",
  filename = "plots/ass4_naivebayes.png"  # Specify the filename and format
)