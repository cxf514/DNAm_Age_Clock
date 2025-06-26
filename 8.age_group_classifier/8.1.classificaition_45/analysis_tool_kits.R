#' @title ROC Curve Plot (Youden Optimized)
#' @description Plot ROC curves using LOOCV results with Youden optimal threshold
#' @param model_list List of models containing LOOCV results
#' @param data Original data (must contain true labels)
#' @export
plot_roc_comparison <- function(model_list, data, true_col, prob_col, positive_class) {
  require(pROC)
  require(ggplot2)
  
  # Extract predictions from each model (keeping original index handling)
  pred_list <- lapply(model_list, function(m) {
    m$pred[m$pred$rowIndex %in% 1:nrow(data), ]
  })
  
  # Calculate ROC curves (with original levels definition)
  roc_list <- lapply(pred_list, function(pred) {
    roc(
      response = pred$obs,
      predictor = pred[[prob_col]],
      levels = c("Old", "Young") 
    )
  })
  names(roc_list) <- names(model_list)
  
  # Generate the plot (with original visualization parameters)
  p <- ggroc(roc_list, aes = "color", size = 1) +
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color = "grey", linetype = "dashed") +
    labs(
      title = "ROC Curves (LOOCV)",
      subtitle = paste("Positive Class:", "Young")
    ) +
    scale_color_manual(
      name = "Model (AUC)",
      values = c("#E64B35", "#4DBBD5", "#00A087"), 
      labels = sprintf(
        "%s = %.3f", 
        names(roc_list), 
        sapply(roc_list, auc)
      )
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.7, 0.2),
      title = element_text(size = 14, color = 'black'),
      legend.text = element_text(size = 12, color = 'black'),
      legend.key.size = unit(20, "pt"),
      legend.title = element_text(size = 12),
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 12)
    )
  return(list(plot = p, roc_list = roc_list))
}

#' @title Generate Youden-Optimized Confusion Matrix
#' @description Create confusion matrix using model-specific Youden optimal threshold
#' @param model_pred Model prediction object from LOOCV
#' @param data Original dataset
#' @param model_name Name of the model for plot title
#' @export
plot_confusion_matrix <- function(model_pred, data, model_name = "GLM") {
  require(pROC)
  require(caret)
  require(ggplot2)
  require(dplyr)
  
  # 1. Extract predicted probabilities and true labels
  pred_prob <- model_pred$pred[, "Young"]
  true_label <- model_pred$pred$obs
  
  # 2. Calculate ROC curve and optimal threshold (Youden index)
  roc_obj <- roc(response = true_label, predictor = pred_prob, levels = c("Old", "Young"))
  optimal_threshold <- coords(roc_obj, "best", best.method = "youden")$threshold
  
  # 3. Generate predicted classes
  pred_class <- factor(
    ifelse(pred_prob > optimal_threshold, "Young", "Old"),
    levels = c("Young", "Old")
  )
  
  # 4. Compute confusion matrix
  conf_matrix <- confusionMatrix(
    data = pred_class,
    reference = true_label,
    positive = "Young"
  )
  
  # 5. Convert to dataframe for visualization
  conf_matrix_df <- as.data.frame(conf_matrix$table)
  colnames(conf_matrix_df) <- c("Predicted", "Actual", "Count")
  
  # Add readable labels (assuming cutoff at age 55)
  conf_matrix_df <- conf_matrix_df %>%
    mutate(
      Predicted = case_when(
        Predicted == "Young" ~ "Predicted <55",
        TRUE ~ "Predicted ≥55"
      ),
      Actual = case_when(
        Actual == "Young" ~ "Actual <55",
        TRUE ~ "Actual ≥55"
      )
    )
  
  # 6. Create heatmap plot
  p <- ggplot(conf_matrix_df, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "white", linewidth = 1, alpha = 0.8) +
    geom_text(aes(label = Count), color = "black", size = 6, fontface = "bold") +
    scale_fill_gradient(low = "#F5F5F5", high = "#006837") +
    labs(
      title = paste("Confusion Matrix (", model_name, " Model)", sep = ""),
      x = "Actual Age Group",
      y = "Predicted Age Group"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  # 7. Extract misclassified samples
  error_samples <- data.frame(
    ID = model_pred$pred$rowIndex[true_label != pred_class],
    True_Label = true_label[true_label != pred_class],
    Pred_Label = pred_class[true_label != pred_class],
    Probability = pred_prob[true_label != pred_class],
    True_Age = data$Age[model_pred$pred$rowIndex[true_label != pred_class]]
  )
  
  # 8. Return all results
  return(list(
    plot = p,
    confusion_matrix = conf_matrix,
    optimal_threshold = optimal_threshold,
    error_samples = error_samples,
    roc_object = roc_obj
  ))
}
