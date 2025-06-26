# Author: cxf
# Date: 2025-06-24
# Objective: SVR age prediction model building

# Environment Setup -------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")

# Load Packages -----------------------------------------------------------
library(dplyr)
library(regnet)
library(caret)
library(glmnet)
library(e1071)
library(splitTools)
library(doParallel)
library(tidyverse)
library(ggpubr)

# Data Preparation --------------------------------------------------------
load("../Rdata/cpg_info.RData")
set.seed(42)
exclude_samples <- c(91, 93)
retain_idx <- setdiff(1:nrow(combined_cpg), exclude_samples)
combined_cpg <- as.matrix(combined_cpg[retain_idx,])
age_info <- as.numeric(combined_info[retain_idx,]$Age)
# Age Grouping ------------------------------------------------------------
age_groups <- cut(
  age_info,
  breaks = c(0, 18, 30, 40, 50, 60, 70, Inf),
  labels = c("A", "B", "C", "D", "E", "F", "G"),
  right = TRUE
)

age_threshold <- 50

# Data Splitting ----------------------------------------------------------
svr_grid <- expand.grid(
  C = 10^seq(-3, 1, length.out = 10),         
  sigma = 10^seq(-1, -3, length.out = 10)  
)

train_idx <- createDataPartition(
  age_groups,   
  p = 0.6,        
  list = FALSE,   
  times = 1     
)

test_idx <- setdiff(1:nrow(combined_cpg), train_idx)
# Data Assignment ---------------------------------------------------------
X_train <- combined_cpg[train_idx, ]
X_test <- combined_cpg[test_idx, ]
y_train <- age_info[train_idx]
y_test <- age_info[test_idx]

train_young <- y_train < age_threshold
train_old <- y_train >= age_threshold
test_young <- y_test < age_threshold
test_old <- y_test >= age_threshold

# Model Training ----------------------------------------------------------
ctrl <- trainControl(method = "cv", number = 10)

# Global Model
model_global <- train(
  x = X_train,
  y = y_train,
  method = "svmRadial",
  tuneGrid = svr_grid,
  trControl = ctrl
)

train_model <- function(x, y, group_name) {
  cat("\nTraining", group_name, "model with all features...")
  train(
    x = x,
    y = y,
    method = "svmRadial",
    tuneGrid = svr_grid,
    trControl = ctrl
  )
}

# Age-specific Models
model_young <- train_model(
  X_train[train_young, ], 
  y_train[train_young],
  "Young"
)

model_old <- train_model(
  X_train[train_old, ], 
  y_train[train_old],
  "Old"
)

# Model Evaluation --------------------------------------------------------
compare_models_full <- function(global_model, 
                                segmented_model_young = NULL, 
                                segmented_model_old = NULL,
                                X_train, y_train,
                                X_test, y_test,
                                train_young_idx, train_old_idx,
                                test_young_idx, test_old_idx) {
  
  results <- data.frame()
  
  evaluate_performance <- function(model, x, y, model_name, group_name, dataset_type) {
    pred <- predict(model, x)
    data.frame(
      Model = model_name,
      Group = group_name,
      Dataset = dataset_type,
      MAE = round(MAE(pred, y), 2),
      RMSE = round(RMSE(pred, y), 2),
      R2 = round(R2(pred, y), 3),
      stringsAsFactors = FALSE
    )
  }
  
  # Global Model Evaluation
  results <- rbind(results,
                   evaluate_performance(global_model, X_train, y_train, "Overall", "All", "Train"),
                   evaluate_performance(global_model, X_test, y_test, "Overall", "All", "Test"),
                   evaluate_performance(global_model, X_train[train_young_idx, ], y_train[train_young_idx],
                                        "Overall", "Young", "Train"),
                   evaluate_performance(global_model, X_test[test_young_idx, ], y_test[test_young_idx],
                                        "Overall", "Young", "Test"),
                   evaluate_performance(global_model, X_train[train_old_idx, ], y_train[train_old_idx],
                                        "Overall", "Old", "Train"),
                   evaluate_performance(global_model, X_test[test_old_idx, ], y_test[test_old_idx],
                                        "Overall", "Old", "Test")
  )
  
  # Segmented Model Evaluation
  if (!is.null(segmented_model_young) && !is.null(segmented_model_old)) {
    results <- rbind(results,
                     evaluate_performance(segmented_model_young, X_train[train_young_idx, ], y_train[train_young_idx],
                                          "Segmented", "Young", "Train"),
                     evaluate_performance(segmented_model_young, X_test[test_young_idx, ], y_test[test_young_idx],
                                          "Segmented", "Young", "Test"),
                     evaluate_performance(segmented_model_old, X_train[train_old_idx, ], y_train[train_old_idx],
                                          "Segmented", "Old", "Train"),
                     evaluate_performance(segmented_model_old, X_test[test_old_idx, ], y_test[test_old_idx],
                                          "Segmented", "Old", "Test")
    )
    
    # Segmented Model Combined Performance
    segmented_all_train_pred <- numeric(length(y_train))
    segmented_all_train_pred[train_young_idx] <- predict(segmented_model_young, X_train[train_young_idx, ])
    segmented_all_train_pred[train_old_idx] <- predict(segmented_model_old, X_train[train_old_idx, ])
    
    segmented_all_test_pred <- numeric(length(y_test))
    segmented_all_test_pred[test_young_idx] <- predict(segmented_model_young, X_test[test_young_idx, ])
    segmented_all_test_pred[test_old_idx] <- predict(segmented_model_old, X_test[test_old_idx, ])
    
    results <- rbind(results,
                     data.frame(
                       Model = "Segmented",
                       Group = "All",
                       Dataset = "Train",
                       MAE = round(MAE(segmented_all_train_pred, y_train), 2),
                       RMSE = round(RMSE(segmented_all_train_pred, y_train), 2),
                       R2 = round(R2(segmented_all_train_pred, y_train), 3),
                       stringsAsFactors = FALSE
                     ),
                     data.frame(
                       Model = "Segmented",
                       Group = "All",
                       Dataset = "Test",
                       MAE = round(MAE(segmented_all_test_pred, y_test), 2),
                       RMSE = round(RMSE(segmented_all_test_pred, y_test), 2),
                       R2 = round(R2(segmented_all_test_pred, y_test), 3),
                       stringsAsFactors = FALSE
                     )
    )
  }
  
  # Sort results
  results <- results[order(results$Model, results$Group, results$Dataset), ]
  cat("\n===== Model Performance Summary =====\n")
  print(results, row.names = FALSE)
  return(results)
}

# Run Evaluation
performance_results_full <- compare_models_full(
  global_model = model_global,
  segmented_model_young = model_young,
  segmented_model_old = model_old,
  X_train = X_train, y_train = y_train,
  X_test = X_test, y_test = y_test,
  train_young_idx = train_young,
  train_old_idx = train_old,
  test_young_idx = test_young,
  test_old_idx = test_old
)

write.csv(performance_results_full, "model_performance_comparison.csv", row.names = FALSE)

# Generate Prediction Tables ----------------------------------------------
generate_result_tables <- function(global_model, 
                                   segmented_model_young = NULL, 
                                   segmented_model_old = NULL,
                                   X_train, y_train,
                                   X_test, y_test,
                                   train_young, train_old,
                                   test_young, test_old) {
  
  # Overall Model Results
  overall_df <- rbind(
    # Training Set
    data.frame(
      True_Age = y_train,
      Predicted_Age = predict(global_model, X_train),
      Dataset = "Train",
      Age_Group = ifelse(y_train < age_threshold, "Young", "Old"),
      Model = "Overall",
      stringsAsFactors = FALSE
    ),
    # Test Set
    data.frame(
      True_Age = y_test,
      Predicted_Age = predict(global_model, X_test),
      Dataset = "Test",
      Age_Group = ifelse(y_test < age_threshold, "Young", "Old"),
      Model = "Overall",
      stringsAsFactors = FALSE
    )
  )
  
  # Segmented Model Results
  segmented_df <- data.frame()
  
  if (!is.null(segmented_model_young) && !is.null(segmented_model_old)) {
    # Training Predictions
    train_pred <- numeric(length(y_train))
    train_pred[train_young] <- predict(segmented_model_young, X_train[train_young, ])
    train_pred[train_old] <- predict(segmented_model_old, X_train[train_old, ])
    
    # Test Predictions
    test_pred <- numeric(length(y_test))
    test_pred[test_young] <- predict(segmented_model_young, X_test[test_young, ])
    test_pred[test_old] <- predict(segmented_model_old, X_test[test_old, ])
    
    segmented_df <- rbind(
      # Training Set
      data.frame(
        True_Age = y_train,
        Predicted_Age = train_pred,
        Dataset = "Train",
        Age_Group = ifelse(y_train < age_threshold, "Young", "Old"),
        Model = "Segmented",
        stringsAsFactors = FALSE
      ),
      # Test Set
      data.frame(
        True_Age = y_test,
        Predicted_Age = test_pred,
        Dataset = "Test",
        Age_Group = ifelse(y_test < age_threshold, "Young", "Old"),
        Model = "Segmented",
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(list(
    overall_results = overall_df,
    segmented_results = segmented_df
  ))
}

# Generate and Save Results
results <- generate_result_tables(
  global_model = model_global,
  segmented_model_young = model_young,
  segmented_model_old = model_old,
  X_train = X_train, y_train = y_train,
  X_test = X_test, y_test = y_test,
  train_young = train_young,
  train_old = train_old,
  test_young = test_young,
  test_old = test_old
)

# get_overall_predicted_age -------------------------------------------------
y_train_pred<-predict(model_global,X_train)
y_test_pred<-predict(model_global,X_test)
svr_train_info<-data_frame('Age'=y_train,'pAge'=y_train_pred)
svr_test_info<-data_frame('Age'=y_test,'pAge'=y_test_pred)
save(svr_train_info,svr_test_info,model_global, model_young, model_old, results, file = "SVR.RData")
  