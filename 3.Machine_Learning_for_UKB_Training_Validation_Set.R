# Machine Learning for UKB (Training/Validation Set) Based on Selected Features

# Load necessary libraries
library(doParallel) 
library(caret)
library(caretEnsemble)
library(MLmetrics)
library(data.table)
library(Boruta) 
library(tidyverse)
library(cvAUC)
library(ggplot2)
library(pROC)
library(skimr)
library(randomForest)

# Set up parallel processing
cl <- makeCluster(20) 
registerDoParallel(cl)
set.seed(123)  # For reproducibility

# Load training dataset
trainingSet <- readRDS("path/to/UKB_training_validation_set.rds")

# Ensure 'incidentAF' is a factor and define case/control labels
trainingSet <- trainingSet %>%
  filter(!is.na(incidentAF)) %>%
  mutate(
    incidentAF = as.factor(incidentAF),
    AF = factor(ifelse(incidentAF == "2", "Case", "Control"), levels = c("Control", "Case"))
  )

# Load selected feature list
feature_list <- readRDS("path/to/feature_selection/final_protein_list.rds")

# Define continuous and categorical feature sets
continuous_features <- c(feature_list) 
categorical_features <- c()  

# Ensure proper column naming
trainingSet <- trainingSet %>%
  rename(sample_id = id) %>%
  mutate(sample_id = as.character(sample_id))

# Split dataset into continuous and categorical subsets
Cont <- trainingSet %>% select(sample_id, all_of(continuous_features), AF)
Cat <- trainingSet %>% select(sample_id, all_of(categorical_features), AF)

# Machine Learning Loop (100 Iterations)
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
  
  # Scale continuous features
  cat_names <- c("sample_id", "AF")
  fit_scale <- preProcess(Continuous_train %>% select(-all_of(cat_names)), method = c("center", "scale"))
  Continuous_train[,-which(names(Continuous_train) %in% cat_names)] <- predict(fit_scale, Continuous_train[,-which(names(Continuous_train) %in% cat_names)])
  Continuous_test[,-which(names(Continuous_test) %in% cat_names)] <- predict(fit_scale, Continuous_test[,-which(names(Continuous_test) %in% cat_names)])
  
  # Merge continuous and categorical data
  Train <- full_join(Continuous_train, Categorical_train, by = "sample_id") %>%
    select(-AF.x) %>%
    rename(AF = AF.y) %>%
    relocate(AF, .after = last_col())
  
  Test <- full_join(Continuous_test, Categorical_test, by = "sample_id") %>%
    select(-AF.x) %>%
    rename(AF = AF.y) %>%
    relocate(AF, .after = last_col())
  
  # Define 10-fold cross-validation
  fitControl_10CV <- trainControl(method = "cv", number = 10, savePredictions = "final", classProbs = TRUE)
  
  # Train Random Forest model
  print(paste("...training RF model for", h))
  fit_model <- train(AF ~ ., data = Train %>% select(-sample_id), method = "rf", trControl = fitControl_10CV, importance = TRUE)
  print(paste("...finished RF model for", h))
  
  # Extract feature importance
  importance_values <- as.data.frame(importance(fit_model$finalModel))
  importance_values$feature <- rownames(importance_values)
  importance_values$iter <- h
  
  # Save feature importance
  write.table(importance_values, file = paste0("path/to/training_validation/feature_importance_", h, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Save ISAF probabilities
  YtestPredProb <- predict(fit_model, Test %>% select(-sample_id), type = "prob")[,2]
  isaf_out <- data.frame(IID = Test$sample_id, Prob = YtestPredProb)
  write.table(isaf_out, file = paste0("path/to/training_validation/isaf_probabilities_", h, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  print(paste("***** done with", h))
}

# Cleanup
stopCluster(cl)