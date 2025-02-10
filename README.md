# Proteomic_AF_ML

This repository provides machine learning code used for a project "A Machine Learning-Based Plasma Protein Risk Score Improves Atrial Fibrillation Prediction Over Clinical and Genomic Models" 

![image](https://github.com/user-attachments/assets/516278d4-8500-464a-8b0a-3f7a65724bf3)

Depicted is an overview of training and testing of Pro-AF, a random forest-based machine learning model to predict incident AF. A nested cross-validation procedure was used to train models and select hyperparameters. The best performing models were then evaluated in a multiply resampled internal test set comprising n=30,632 individuals (internal test set), and a hold out test set comprising n=13,998 individuals (hold out test set). The sampling method ensures that no model is tested on individuals used in model training or validation.

# The code consists of the following components:
1. Cohort Splitting – Dividing the dataset into training/validation and holdout sets.
2. Feature Selection – Selecting relevant features from nearly 3,000 proteins.
3. Machine Learning Model Training – Training and validating the model on the UKB dataset using the selected features.
4. Model Application – Applying the trained model to the holdout set for evaluation.
