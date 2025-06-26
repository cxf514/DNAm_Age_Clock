# Author: cxf
# Date: 2025-06-12
# Objective: Robust variable selection for age prediction using DNA methylation data
# Clear workspace and set options
rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")
# Load required packages --------------------------------------------------
library(dplyr)
library(regnet)
library(caret)
library(glmnet)
library(e1071)
library(splitTools)
library(doParallel)
library(tidyverse)
library(ggpubr)
# Load and prepare data ---------------------------------------------------
load("../../0.raw_data/cpg_info.Rdata")
combined_cpg <- as.matrix(combined_cpg)

# Create age groups -------------------------------------------------------
age_info <- as.numeric(combined_info$Age)
age_groups <- cut(
  age_info,
  breaks = c(0, 18, 30, 40, 50, 60, 70, Inf),
  labels = c("A", "B", "C", "D", "E", "F", "G"),
  right = TRUE
)

# Configuration parameters ------------------------------------------------
n_outer_folds <- 5
n_repeats <- 50
alpha_grid <- seq(0.1, 0.9, by = 0.1)
lambda_seq <- 10^seq(-3, 3, length.out = 50)
sigma_vals <- 10^seq(-3, -1, length.out = 10)
c_vals <- 10^seq(1, 3, length.out = 10)
svr_grid <- expand.grid(
  C = c_vals,         
  sigma = sigma_vals  
)
exclude_samples <- c(91, 93)
set.seed(42)
cv_folds <- create_folds(
  y = age_groups,
  k = n_outer_folds,
  m_rep = n_repeats,
  type = "stratified",
  invert = FALSE
)

# Initialize result storage -----------------------------------------------
all_results <- list()

selected_counts <- list(
  lasso = numeric(0),
  rlasso = numeric(0)
)
feature_counts_by_model <- list(
  lasso = list(),
  rlasso = list()
)

# add parallel ------------------------------------------------------------
n_cores <- 50
cl <- makePSOCKcluster(n_cores)
clusterEvalQ(cl, {
  .libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")
  library(caret) 
  library(kernlab)
})
registerDoParallel(cl)

# Main cross-validation loop ----------------------------------------------

for(i in 1:length(cv_folds)) {
  
##1.Data partitioning----
  train_idx <- cv_folds[[i]]
  test_idx <- setdiff(1:nrow(combined_cpg), train_idx)
##2.Outlier removal----
  train_idx <- setdiff(train_idx, exclude_samples)
  test_idx <- setdiff(test_idx, exclude_samples)
##3.Data assignment----
  X_train <- combined_cpg[train_idx, ]
  X_test <- combined_cpg[test_idx, ]
  y_train <- age_info[train_idx]
  y_test <- age_info[test_idx]
  
# Model Training and Evaluation-----------------------------------
  ## 1. LASSO Regression-----
  cv_lasso <- cv.glmnet(X_train, y_train,
                        alpha = 1,
                        lambda = lambda_seq,
                        nfolds = 10,
                        type.measure = "mae")
  best_lambda_lasso <- cv_lasso$lambda.min
  model_lasso <- glmnet(X_train, y_train,
                        alpha = 1,
                        lambda = best_lambda_lasso)
  y_pred_lasso <- predict(model_lasso, newx = X_test)
  mae_lasso <- mean(abs(y_pred_lasso - y_test))
  coefs_lasso <- coef(model_lasso)[-1, ]
  nonzero_lasso <- names(coefs_lasso[coefs_lasso != 0])
  ## 2. Robust LASSO Regression-----
  cv_model <- cv.regnet(
    X_train, y_train,
    penalty = "lasso",
    response = "continuous",
    robust = TRUE,
    folds = 10,
    lamb.1 = lambda_seq,
    verbo = FALSE
  )
  lambda_rlasso <- cv_model$lambda[1]
  final_model <- regnet(
    X_train, y_train,
    lamb.1 = lambda_rlasso,
    penalty = "lasso",
    response = "continuous",
    robust = TRUE
  )
  X_train_mean <- apply(X_train, 2, mean)
  X_train_sd <- apply(X_train, 2, sd)
  X_test_scaled <- scale(X_test,
                         center = X_train_mean,
                         scale = X_train_sd)
  intercept <- final_model$coeff["Intercept"]  
  coefs_rlasso <- final_model$coeff[-1] 
  y_pred <- intercept + X_test_scaled %*% coefs_rlasso
  mae_rlasso <- mean(abs(y_pred - y_test))
  nonzero_rlasso <- names(coefs_rlasso)[coefs_rlasso != 0]
  ## 3. Ridge Regression-----
  cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0, lambda = lambda_seq, nfolds = 10, type.measure = "mae")
  y_pred_ridge <- as.vector(predict(cv_ridge, newx = X_test, s = "lambda.min"))
  best_lambda_ridge <- cv_ridge$lambda.min
  mae_ridge <- mean(abs(y_pred_ridge - y_test))
  ## 4. Elastic Net-----
  results <- data.frame(alpha = numeric(), lambda = numeric(), mae = numeric())
  for (a in alpha_grid) {
    cv_elastic <- cv.glmnet(
      X_train, 
      y_train, 
      alpha = a, 
      lambda = lambda_seq,
      nfolds = 10,
      type.measure = "mae"
    )
    results <- rbind(results, data.frame(
      alpha = a,
      lambda = cv_elastic$lambda.min,
      mae = min(cv_elastic$cvm)
    ))
  }
  best_index <- which.min(results$mae)
  best_alpha <- results$alpha[best_index]
  best_lambda <- results$lambda[best_index]
  model_elastic <- glmnet(X_train, y_train, alpha = best_alpha, lambda = best_lambda)
  y_pred_elastic <- predict(model_elastic, newx = X_test)
  mae_elastic <- mean(abs(y_pred_elastic - y_test))
  
  ## 5. Support Vector Regression-----
  svr_model <- train(
    x = X_train, y = y_train,
    method = "svmRadial",
    tuneGrid = svr_grid,
    metric = "MAE",       
    maximize = FALSE,  
    trControl = trainControl(method = "cv", number = 10)
  )
  y_pred_svr <- predict(svr_model, X_test)
  svr_C = svr_model$bestTune$C
  svr_sigma = svr_model$bestTune$sigma
  mae_svr <- mean(abs(y_pred_svr - y_test))
  
  # result summary ----------------------------------------------------------
  all_results[[i]] <- list(
    iteration = i,
    # 1.Lasso
    mae_lasso = mae_lasso,
    lambda_lasso = best_lambda_lasso,
    features_lasso = nonzero_lasso,
    # 2.robust_Lasso
    mae_rlasso = mae_rlasso,
    lambda_rlasso = lambda_rlasso,
    features_rlasso = nonzero_rlasso,
    # 3.Ridge
    mae_ridge = mae_ridge,
    lambda_ridge = best_lambda_ridge,
    # 4.Elastic
    mae_elastic = mae_elastic,
    lambda_elastic = best_lambda,
    alpha_elastic = best_alpha,
    # 5.svr
    mae_svr = mae_svr,
    svr_C = svr_model$bestTune$C,
    svr_sigma = svr_model$bestTune$sigma
  )
  
  
  # feature_frequency -------------------------------------------------------
  for (model_type in c("lasso", "rlasso")) {
    features <- all_results[[i]][[paste0("features_", model_type)]]
    selected_counts[[model_type]] <- c(selected_counts[[model_type]], length(features))
    for (feature in features) {
      if (feature %in% names(feature_counts_by_model[[model_type]])) {
        feature_counts_by_model[[model_type]][[feature]] <- feature_counts_by_model[[model_type]][[feature]] + 1
      } else {
        feature_counts_by_model[[model_type]][[feature]] <- 1
      }
    }
  }
}

stopCluster(cl)



# summary report ----------------------------------------------------------
# 1. Pairwise t-tests for Model Performance (MAE)----
mae_lasso <- sapply(all_results, function(x) x$mae_lasso)
mae_rlasso <- sapply(all_results, function(x) x$mae_rlasso)
mae_ridge <- sapply(all_results, function(x) x$mae_ridge)
mae_elastic <- sapply(all_results, function(x) x$mae_elastic)
mae_svr <- sapply(all_results, function(x) x$mae_svr)

## Create a data frame of MAE values for all models----
mae_data <- data.frame(
  Lasso = mae_lasso,
  rLasso = mae_rlasso,
  Ridge = mae_ridge,
  Elastic_Net = mae_elastic,
  SVR = mae_svr
)

## Perform paired t-tests for each model pair----
model_pairs <- combn(names(mae_data), 2, simplify = FALSE)
pairwise_tests <- lapply(model_pairs, function(pair) {
  test_result <- t.test(
    mae_data[[pair[1]]], 
    mae_data[[pair[2]]], 
    paired = TRUE
  )
  p_value_formatted <- ifelse(
    test_result$p.value < 1e-16, 
    "<1e-16",
    format(test_result$p.value, scientific = TRUE, digits = 3)
  )
  significance <- ifelse(
    test_result$p.value < 0.001, "***",
    ifelse(test_result$p.value < 0.01, "**",
           ifelse(test_result$p.value < 0.05, "*", "(not sig)")
    )
  )
  data.frame(
    Model1 = pair[1],
    Model2 = pair[2],
    Mean_Diff = round(test_result$estimate, 3),
    CI_Lower = round(test_result$conf.int[1], 3),
    CI_Upper = round(test_result$conf.int[2], 3),
    t_statistic = round(test_result$statistic, 3),
    df = test_result$parameter,
    p_value = p_value_formatted,
    Significance = significance 
  )
})
## Print a nicely formatted table---- 
results_table <- do.call(rbind, pairwise_tests)

write.table(results_table,file="paired-t-tests-rm.txt",sep="\t",row.names = F,col.names = T)

knitr::kable(results_table, align = "c")

# 2. Summarize feature selection frequency for Lasso and rLasso----
lasso_summary <- data.frame(
  feature = names(feature_counts_by_model[["lasso"]]),
  count = unlist(feature_counts_by_model[["lasso"]]),
  frequency = unlist(feature_counts_by_model[["lasso"]]) / 250,
  stringsAsFactors = FALSE
) %>% arrange(desc(count))
rlasso_summary <- data.frame(
  feature = names(feature_counts_by_model[["rlasso"]]),
  count = unlist(feature_counts_by_model[["rlasso"]]),
  frequency = unlist(feature_counts_by_model[["rlasso"]]) / 250,
  stringsAsFactors = FALSE
) %>% arrange(desc(count))

write.table(lasso_summary,file="feature_frequency_lasso_rm.txt",sep="\t",row.names = F,col.names = T)
write.table(rlasso_summary,file="feature_frequency_rlasso_rm.txt",sep="\t",row.names = F,col.names = T)

feature_stats <- list(
  lasso = c(
    mean = mean(selected_counts$lasso), 
    sd = sd(selected_counts$lasso)
  ),
  rlasso = c(
    mean = mean(selected_counts$rlasso), 
    sd = sd(selected_counts$rlasso)
  )
)

# Create boxplot with statistics
mae_data <- tibble(
  MAE = c(mae_rlasso,mae_lasso, mae_ridge, mae_elastic,mae_svr),
  Model = factor(rep(c("rLasso","Lasso","Ridge","Elastic Net","SVR"), 
                     times = c(length(mae_lasso),length(mae_rlasso),
                               length(mae_ridge), length(mae_elastic),
                               length(mae_svr))),
                 levels = c("rLasso","Lasso","Elastic Net", "Ridge","SVR")) 
)

model_colors <- c(
  "rLasso" = "#E41A1C",        # Original red
  "Lasso" = "#377EB8",         # Original blue
  "Elastic Net" = "#4DAF4A",   # Complementary green
  "Ridge" = "#984EA3",         # Purple (analogous to blue)
  "SVR" = "#FF7F00"            # Orange (complementary to blue)
)

ggplot(mae_data, aes(x = Model, y = MAE, color = Model)) +
  geom_boxplot(width = 0.7, alpha = 0.8) +
  stat_compare_means(comparisons = list(c("rLasso", "Lasso"),
                                        c("Elastic Net", "Lasso"),
                                        c("Elastic Net", "Ridge"),
                                        c("Lasso","Ridge"),
                                        c("SVR", "Ridge")),
                     method = "t.test", 
                     paired = TRUE,
                     label = "p.signif") +
  scale_color_manual(values  = model_colors) +
  labs(y = "Mean Absolute Error",x='') +
  theme_classic()+
  theme(legend.position='none',panel.grid=element_blank())+
  theme(
    title = element_text(size=20,color = 'black'),
    legend.text = element_text(size=14,color = 'black'),
    legend.key.size = unit(25, "pt"),
    legend.title = element_text(size = 14),
    axis.text=element_text(size=14,color = 'black'),
    axis.title.x =element_text(size=14), 
    axis.title.y=element_text(size=14))

ggsave("cpg_formal_removal.pdf",width = 190, height = 110, units = "mm") 
save.image(file = 'cpg_formal_removal.Rdata')
