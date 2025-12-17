# SPDX-License-Identifier: GPL-3.0-or-later

# Split Cohort into Training/Validation and Holdout Sets

# Load necessary library
library(caret)

# Load the dataset (Replace 'your_data_path' with the actual file path)
protein_imputed <- readRDS("path/to/your_dataset.rds")

# Ensure 'incidentAF' exists in the dataset before proceeding
if (!"incidentAF" %in% colnames(protein_imputed)) {
  stop("Error: The dataset does not contain 'incidentAF'. Please check the data.")
}

# Set seed for reproducibility
set.seed(123)

# Create a 70/30 split (70% Training/Validation, 30% Holdout)
splitIndex <- createDataPartition(protein_imputed$incidentAF, p = 0.7, list = FALSE)

# Define training/validation and holdout sets
trainingSet <- protein_imputed[splitIndex, ]
validationSet <- protein_imputed[-splitIndex, ]

# Display class distribution to ensure balance
cat("Training Set Distribution:\n")
print(table(trainingSet$incidentAF))

cat("\nValidation Set Distribution:\n")
print(table(validationSet$incidentAF))

# Save the split datasets (Replace 'your_save_path' with actual save directory)
saveRDS(trainingSet, "path/to/training_validation_set.rds")
saveRDS(validationSet, "path/to/holdout_set.rds")

cat("\nDatasets successfully saved.\n")


