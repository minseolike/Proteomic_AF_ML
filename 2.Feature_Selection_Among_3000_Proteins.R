# SPDX-License-Identifier: GPL-3.0-or-later

# Feature Selection (Among 3000 Proteins)

# Load necessary libraries
library(doParallel) 
library(caret)
library(Boruta) 
library(tidyverse)

# Set up parallel processing for efficiency
cl <- makeCluster(20)  # Adjust based on system capabilities
registerDoParallel(cl)
set.seed(123)  # For reproducibility

# Load the training dataset (Replace with actual file path)
trainingSet <- readRDS("path/to/UKB_training_validation_set.rds")

# Ensure 'incidentAF' is a factor
trainingSet <- trainingSet %>%
  filter(!is.na(incidentAF)) %>%
  mutate(
    incidentAF = as.factor(incidentAF),
    AF = ifelse(incidentAF == "2", "Case", "Control"),
    AF = factor(CAD, levels = c("Control", "Case"))
  )

# Define feature columns (assuming proteins start at column 103 and end at 3013)
col_names <- colnames(trainingSet)[103:3013]  # Adjust if needed

# Split into continuous and categorical datasets
continuous_features <- col_names
categorical_features <- c()  # Define categorical features if needed

# Select continuous and categorical features
Continuous <- trainingSet %>% select(sample_id, all_of(continuous_features), AF)
Categorical <- trainingSet %>% select(sample_id, all_of(categorical_features), AF)

# Store original data for repeated feature selection
Cont <- Continuous
Cat <- Categorical

# Feature Selection Loop (100 Iterations)
for (h in 1:100) {
  
  Continuous <- Cont
  Categorical <- Cat
  
  print(paste("...starting iteration", h))
  set.seed(h)
  
  # Sample 80% of cases and 2x controls for training
  I <- Continuous$sample_id[sample(which(Continuous$AF == "Case"), size = round(length(which(Continuous$AF == "Case")) * 0.8))]
  I <- c(I, sample(Continuous$sample_id[which(Continuous$AF == "Control")], size = length(I) * 2))
  Continuous_train <- Continuous %>% filter(sample_id %in% I)
  Categorical_train <- Categorical %>% filter(sample_id %in% I)
  
  # Use remaining 20% for testing
  J <- setdiff(Continuous$sample_id, I)
  I <- J[J %in% Continuous$sample_id[Continuous$AF == "Case"]]
  I <- c(I, sample(Continuous$sample_id[Continuous$AF == "Control"], size = length(I) * 2))
  Continuous_test <- Continuous %>% filter(sample_id %in% I)
  Categorical_test <- Categorical %>% filter(sample_id %in% I)
  
  # Scale continuous data based on training set
  cat_names <- c("sample_id", "AF")
  fit_scale <- preProcess(Continuous_train %>% select(-all_of(cat_names)), method = c("center", "scale"))
  Continuous_train[,-which(names(Continuous_train) %in% cat_names)] <- predict(fit_scale, Continuous_train[,-which(names(Continuous_train) %in% cat_names)])
  Continuous_test[,-which(names(Continuous_test) %in% cat_names)] <- predict(fit_scale, Continuous_test[,-which(names(Continuous_test) %in% cat_names)])
  
  # Feature selection using Boruta
  BB <- Boruta(x = Continuous_train %>% select(-all_of(cat_names)), y = as.factor(Continuous_train$AF))
  selected_features <- names(BB$finalDecision[BB$finalDecision != "Rejected"])
  
  feature_selection <- data.frame(feature = selected_features)
  
  # Save feature selection results for this iteration
  saveRDS(feature_selection, file = paste0("path/to/feature_selection/feature_selection_proteins_", h, ".rds"))
}

# Aggregate Feature Selections Across Iterations
feature_files <- list.files("path/to/feature_selection/", pattern = "feature_selection_proteins_.*\\.rds", full.names = TRUE)
feature_selections <- lapply(feature_files, readRDS)

# Combine feature selection results
featimp_out <- bind_rows(feature_selections)

# Process feature selection results
featimp_out <- featimp_out %>%
  mutate(feature = rownames(featimp_out)) %>%
  rename(confirmed = feature) %>%
  filter(confirmed == "Confirmed")

# Count occurrences of each feature across iterations
feature_counts <- featimp_out %>% count(feature)

# Retain features appearing in at least 20% of iterations
frequent_features <- feature_counts %>% filter(n >= 20) %>% select(feature)

# Save final feature list
feature_list <- as.character(frequent_features$feature)
saveRDS(feature_list, "path/to/feature_selection/final_protein_list.rds")

