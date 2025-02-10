# Machine Learning for Hold-Out Set (Final Model)

# Load necessary libraries
library(doParallel) 
library(caret)
library(caretEnsemble)
library(MLmetrics)
library(data.table)
library(tidyverse)
library(cvAUC)
library(pROC)
library(skimr)
library(randomForest)

# Set up parallel processing
cl <- makeCluster(20) 
registerDoParallel(cl)
set.seed(123)  # For reproducibility

# Load training and validation datasets
trainingSet <- readRDS("path/to/UKB_training_validation_set.rds")
validationSet <- readRDS("path/to/UKB_holdout_set.rds")

# Rename identifier column
trainingSet <- trainingSet %>% rename(sample_id = id)
validationSet <- validationSet %>% rename(sample_id = id)

# Filter and process datasets
trainingSet <- trainingSet %>%
  filter(!is.na(incidentAF)) %>%
  mutate(AF = factor(ifelse(incidentAF == "2", "Case", "Control"), levels = c("Control", "Case")))

validationSet <- validationSet %>%
  filter(!is.na(incidentAF)) %>%
  mutate(AF = factor(ifelse(incidentAF == "2", "Case", "Control"), levels = c("Control", "Case")))

# Load selected feature list
feature_list <- readRDS("path/to/feature_selection/final_protein_list.rds")

# Define continuous and categorical feature sets
continuous_features <- c(feature_list) 
categorical_features <- c()

# Ensure proper column naming
trainingSet$sample_id <- as.character(trainingSet$sample_id)
validationSet$sample_id <- as.character(validationSet$sample_id)

# Split dataset into continuous and categorical subsets
Cont_t <- trainingSet %>% select(sample_id, all_of(continuous_features), AF)
Cont_v <- validationSet %>% select(sample_id, all_of(continuous_features), AF)

# Machine Learning Loop (100 Iterations)
for (h in 1:100) {
  
  Continuous_t <- Cont_t
  Continuous_v <- Cont_v
  
  print(paste("...starting iteration", h))
  set.seed(h)
  
  # Use entire holdout set for testing
  Continuous_train <- Continuous_t
  Continuous_test <- Continuous_v
  
  # Scale continuous features using training set parameters
  cat_names <- c("sample_id", "AF")
  fit_scale <- readRDS(paste0("path/to/training_validation/fit_scale_", h, ".RDS"))
  Continuous_train[,-which(names(Continuous_train) %in% cat_names)] <- predict(fit_scale, Continuous_train[,-which(names(Continuous_train) %in% cat_names)])
  Continuous_test[,-which(names(Continuous_test) %in% cat_names)] <- predict(fit_scale, Continuous_test[,-which(names(Continuous_test) %in% cat_names)])
  
  # Train Random Forest model using pre-trained model
  print(paste("...loading trained RF model for", h))
  fit_model <- readRDS(paste0("path/to/training_validation/fit_model_", h, ".RDS"))
  
  # Evaluate testing performance
  YtestPredProb <- predict(fit_model, Continuous_test %>% select(-sample_id), type = "prob")[,2]
  YtestPredRaw <- predict(fit_model, Continuous_test %>% select(-sample_id), type = "raw")
  YtestTrue <- factor(Continuous_test$AF, levels = c("Control", "Case"))
  ConfMtest <- confusionMatrix(data = YtestPredRaw, reference = YtestTrue, positive = "Case")
  
  # Calculate and save results
  metrics <- data.frame(AUC = AUC(predictions = as.numeric(YtestPredProb), labels = as.numeric(as.character(YtestTrue))))
  write.table(metrics, file = paste0("path/to/holdout/metrics_", h, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Extract and save feature importance
  importance_values <- as.data.frame(importance(fit_model$finalModel))
  importance_values$feature <- rownames(importance_values)
  importance_values$iter <- h
  write.table(importance_values, file = paste0("path/to/holdout/feature_importance_", h, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Save ISAF probabilities
  isaf_out <- data.frame(IID = Continuous_test$sample_id, Prob = YtestPredProb)
  write.table(isaf_out, file = paste0("path/to/holdout/isaf_probabilities_", h, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  print(paste("***** done with", h))
}

# Cleanup
stopCluster(cl)
