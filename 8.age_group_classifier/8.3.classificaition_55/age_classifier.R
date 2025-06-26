# Load required libraries
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(dplyr)
library(doParallel)
rm(list = ls())
source('analysis_tool_kits.R')
# =============================================
# DATA PREPARATION
# =============================================
# Set reproducible seed and load data
set.seed(42)
load("../../Rdata/cpg_info.Rdata")
# Process CpG data
colnames(combined_cpg) <- paste0("CpG_", seq_len(ncol(combined_cpg)))
combined_cpg <- as.matrix(combined_cpg)
# Remove excluded samples
exclude_samples <- c(91, 93)
retained_id <- setdiff(seq_len(nrow(combined_cpg)), exclude_samples)
combined_cpg <- combined_cpg[retained_id, ]
combined_info <- combined_info[retained_id, ]
# Create age groups with proper factor levels
age_groups <- factor(
  ifelse(combined_info$Age < 55, "Young", "Old"),
  levels = c("Young", "Old")  # Young = positive class (<55)
)
# Final dataset structure
data <- data.frame(
  Age = combined_info$Age,
  Age_Group = age_groups,
  combined_cpg,
  check.names = FALSE
)
# =============================================
# PARALLEL PROCESSING SETUP
# =============================================
cl <- makePSOCKcluster(30)
clusterEvalQ(cl, {
  .libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")
})
registerDoParallel(cl)
# =============================================
# MODEL TRAINING CONFIGURATION
# =============================================
ctrl <- trainControl(
  method = "LOOCV",
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  allowParallel = TRUE
)
# =============================================
# MODEL 1: ELASTIC NET REGRESSION
# =============================================
glm_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.1),
  lambda = 10^seq(-3, 0, 0.1)
)
set.seed(42)
model_glm <- train(
  Age_Group ~ . - Age,
  data = data,
  method = "glmnet",
  family = "binomial",
  tuneGrid = glm_grid,
  trControl = ctrl,
  metric = "ROC"
)
# =============================================
# MODEL 2: RANDOM FOREST
# =============================================
rf_grid <- expand.grid(
  mtry = seq(1, floor(sqrt(ncol(data) - 2)) + 3, by = 1),
  splitrule = c("gini", "extratrees"),  
  min.node.size = c(1, 3, 5, 10, 15)
)
set.seed(42)
model_rf <- train(
  Age_Group ~ . - Age,
  data = data,
  method = "ranger",
  importance = "permutation",
  num.trees = 1500,
  replace = FALSE,
  sample.fraction = 0.8,
  tuneGrid = rf_grid,
  trControl = ctrl,
  metric = "ROC"
)
# =============================================
# MODEL 3: SUPPORT VECTOR MACHINE
# =============================================
svm_grid <- expand.grid(
  sigma = exp(seq(-5, 1, length.out = 20)),
  C = 2^seq(-3, 1, length.out = 20)
)
set.seed(42)
model_svm <- train(
  Age_Group ~ . - Age,
  data = data,
  method = "svmRadial",
  tuneGrid = svm_grid,
  trControl = ctrl,
  metric = "ROC",
  prob.model = TRUE
)
# =============================================
# CLEAN UP PARALLEL PROCESSING
# =============================================
stopCluster(cl) 
# =============================================
# MODEL EVALUATION
# =============================================
# Store all models for comparison
models <- list(
  LogisticR = model_glm,
  RandomForest = model_rf,
  SVM = model_svm
)
# Generate ROC comparison plot
roc_result <- plot_roc_comparison(
  model_list = models,
  data = data,
  true_col = "Age_Group",
  prob_col = "Young",
  positive_class = "Young"
)
# View ROC results
roc_result$plot
ggsave("ROC_curves.pdf",roc_result$plot,width = 120, height = 120, units = "mm") 
roc_result$roc_list
# =============================================
# MODEL-SPECIFIC EVALUATIONS
# =============================================
# Evaluate ElasticNet model
glm_eval <- plot_confusion_matrix(
  model_pred = model_glm, 
  data = data,
  model_name = "GLM"
)
# Evaluate SVM model
svm_eval <- plot_confusion_matrix(
  model_pred = model_svm, 
  data = data,
  model_name = "SVM"
)
# Evaluate Random Forest model
rf_eval <- plot_confusion_matrix(
  model_pred = model_rf, 
  data = data,
  model_name = "RF"
)
ggsave("confusion_matrix_lr.pdf",glm_eval$plot,width = 150, height = 120, units = "mm") 
ggsave("confusion_matrix_rf.pdf",rf_eval$plot,width = 150, height = 120, units = "mm") 
ggsave("confusion_matrix_svm.pdf",svm_eval$plot,width = 150, height = 120, units = "mm") 
# =============================================
# check misclassified samples
# =============================================
write.csv(glm_eval$error_samples,"misclassified_samples_glm.csv", row.names = FALSE)
write.csv(rf_eval$error_samples,"misclassified_samples_rf.csv", row.names = FALSE)
write.csv(svm_eval$error_samples,"misclassified_samples_svm.csv", row.names = FALSE)

save.image('age_classifier.RData')
