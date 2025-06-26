# =========================================================================
# DNA Methylation Age Prediction: Model Comparison
# 
# Key Features:
# - Strict reproducibility controls
# - Clean code structure while preserving original variable names
# - Maintained parallel processing with proper seeding
# - Preserved statistical and plotting functions
# =========================================================================

# 1. Initialization with Reproducibility Controls -------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
seed_1=42
# Set consistent random seeds
set.seed(seed_1)  # Master seed for overall reproducibility
.libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")
source('../analysis_funcs.R')
# Load required packages
library(caret)       # For SVR and CV folds
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping
library(rstatix)     # Statistical tests
library(ggpubr)      # Publication-quality plots
library(splitTools)  # Data splitting
library(patchwork)   # picture composition
library(broom)       # tidy output
# =========================================================================
# 2. Data Preparation -----------------------------------------------------
# =========================================================================
load("../Rdata/cpg_info.Rdata")

# Prepare methylation data matrix (preserving original variable names)
colnames(combined_cpg) <- paste0("CpG_", 1:17)
combined_cpg <- as.matrix(combined_cpg)
age_info <- as.numeric(combined_info$Age)

# Create age groups
age_groups <- cut(age_info, 
                  breaks = c(0, 18, 30, 40, 50, 60, 70, Inf),
                  labels = c("A", "B", "C", "D", "E", "F", "G"),
                  right = TRUE)
# Experimental parameters
n_outer_folds <- 5      # Number of CV folds
n_repeats <- 50         # Repetitions per fold
age_threshold <- 55     # Threshold for age grouping
exclude_samples <- c(91, 93)  # Samples to exclude
# 3. Generate Reproducible CV Folds ---------------------------------------
set.seed(seed_1)  # Reset seed for consistent fold generation
cv_folds <- create_folds(
  y = age_groups,
  k = n_outer_folds,
  m_rep = n_repeats,
  type = "stratified"
)
# 4. Single-threaded Analysis ---------------------------------------------
all_results <- list()
for (i in seq_along(cv_folds)) {
  set.seed(seed_1 + i)
  # Unpack current fold indices
  fold_indices <- cv_folds[[i]]
  
  # Create train/test splits
  train_idx <- setdiff(fold_indices, exclude_samples)
  test_idx <- setdiff(1:nrow(combined_cpg), c(fold_indices, exclude_samples))
  
  X_train <- combined_cpg[train_idx, ]
  X_test <- combined_cpg[test_idx, ]
  y_train <- age_info[train_idx]
  y_test <- age_info[test_idx]
  
  # Age subgroup indices
  test_young_idx <- which(y_test < age_threshold)
  test_old_idx <- which(y_test >= age_threshold)
  
  # Control object with fixed CV seeds
  ctrl <- trainControl(
    method = "cv",
    number = 10,
    seeds = lapply(1:11, function(x) rep(seed_1, 3)) 
  )
  
  # Model training --------------------------------------------------------
  svr_overall <- train(x = X_train, y = y_train, 
                       method = "svmRadial", trControl = ctrl)
  
  # Age-specific models (if sufficient samples)
  train_young_idx <- which(y_train < age_threshold)
  if (length(train_young_idx) > 1) {
    svr_young <- train(x = X_train[train_young_idx, ], 
                       y = y_train[train_young_idx],
                       method = "svmRadial", trControl = ctrl)
  }
  
  train_old_idx <- which(y_train >= age_threshold)
  if (length(train_old_idx) > 1) {
    svr_old <- train(x = X_train[train_old_idx, ], 
                     y = y_train[train_old_idx],
                     method = "svmRadial", trControl = ctrl)
  }
  
  # MAE calculation helper
  calc_mae <- function(model, X, y) {
    if (length(y) == 0) return(NA_real_)
    pred <- predict(model, X)
    mean(abs(pred - y))
  }
  
  # Store results (maintain original structure)
  all_results[[i]] <- list(
    SVR = list(
      Overall = list(
        Total = calc_mae(svr_overall, X_test, y_test),
        Young = if (length(test_young_idx)) calc_mae(svr_overall, X_test[test_young_idx, ], y_test[test_young_idx]) else NA,
        Old = if (length(test_old_idx)) calc_mae(svr_overall, X_test[test_old_idx, ], y_test[test_old_idx]) else NA
      ),
      Segmented = list(
        Young = if (exists("svr_young") && length(test_young_idx)) calc_mae(svr_young, X_test[test_young_idx, ], y_test[test_young_idx]) else NA,
        Old = if (exists("svr_old") && length(test_old_idx)) calc_mae(svr_old, X_test[test_old_idx, ], y_test[test_old_idx]) else NA
      )
    )
  )
  
  # Progress tracking
  if (i %% 10 == 0) message(sprintf("Completed %d/%d folds", i, length(cv_folds)))
}
# 5. Results Processing ---------------------------------------------------
# Convert to tidy format
result_df <- do.call(rbind, lapply(seq_along(all_results), function(i) {
  res <- all_results[[i]]
  data.frame(
    Fold = i,
    Algorithm = "SVR",
    Model_Type = rep(c("Overall", "Segmented"), c(3, 2)),
    Age_Group = c("Total", "Young", "Old", "Young", "Old"),
    MAE = c(res$SVR$Overall$Total, res$SVR$Overall$Young, res$SVR$Overall$Old,
            res$SVR$Segmented$Young, res$SVR$Segmented$Old),
    stringsAsFactors = FALSE
  )
}))
# 6. Visualization & Statistics -------------------------------------------
mae_data <- result_df %>%
  filter(Age_Group %in% c("Young", "Old")) %>%
  mutate(
    Age_Group = factor(Age_Group, levels = c("Young", "Old"),
                       labels = c("Young (<55)", "Old (≥55)")),
    Model = factor(Algorithm),
    Strategy = factor(Model_Type)
  )
# Generate plots
plots <- list(
  p_young_svr = make_mae_plot(mae_data, "Young (<55)", "SVR"),
  p_old_svr = make_mae_plot(mae_data, "Old (≥55)", "SVR")
)
final_plot <- plots$p_young_svr / plots$p_old_svr +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))
# Perform statistical tests
test_results <- calculate_pairwise_tests(result_df)
# 7. Save Outputs ---------------------------------------------------------
write.table(test_results, "SVR_overall_vs_segmented.txt", sep = "\t", row.names = FALSE)
ggsave("SVR_overall_vs_segmented.pdf", final_plot, width = 7, height = 9)
save.image('SVR_overall_vs_segmented.RData')