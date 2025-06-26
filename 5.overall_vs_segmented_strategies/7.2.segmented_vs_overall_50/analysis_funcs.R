#' Create a Paired MAE Comparison Plot
#'
#' Generates a publication-quality plot comparing MAE (Mean Absolute Error) between 
#' two modeling strategies (Overall vs. Segmented) for a specific age group and model type,
#' showing paired samples with statistical significance.
#'
#' @param data Dataframe containing MAE results (must contain columns: Age_Group, Model, Strategy, MAE, Fold)
#' @param age_group Character string specifying age group to plot (e.g. "Young (<55)")
#' @param model_name Character string specifying model type (e.g. "SVR")
#' @param color_palette Vector of two colors for the strategies (default: red/blue)
#' @param jitter_width Numeric controlling horizontal point spread (default: 0.1)
#' @param segment_alpha Numeric (0-1) for line transparency (default: 0.8)
#' @param base_font_size Numeric base font size (default: 14)
#'
#' @return A ggplot object showing:
#'   - Boxplots of MAE distributions
#'   - Paired points connected by lines
#'   - p-value from paired t-test
#'
#' @details
#' Key visual features:
#' 1. Jittered points show individual fold results
#' 2. Gray lines connect paired observations
#' 3. Transparent boxplots show group distributions
#' 4. Color-coded by strategy
#' 5. Automatically calculates and displays significance
#'
#' @examples
#' make_mae_plot(mae_data, "Young (<55)", "SVR")
#' make_mae_plot(mae_data, "Old (≥55)", "SVR", color_palette = c("darkgreen", "orange"))

make_mae_plot <- function(data, age_group, model_name, 
                          color_palette = c("#E41A1C", "#377EB8"),
                          jitter_width = 0.1, 
                          segment_alpha = 0.3,
                          base_font_size = 14) {
  
  # Preprocess and validate input data
  plot_data <- data %>%
    filter(Age_Group == age_group, Model == model_name) %>%
    mutate(Strategy = factor(Strategy, levels = c("Overall", "Segmented")))
  
  if (nrow(plot_data) == 0) {
    stop("No data available for the specified age group and model")
  }
  
  # Create reproducible jittering
  set.seed(123)  # Ensures consistent jitter across runs
  jittered_points <- plot_data %>%
    group_by(Strategy) %>%
    mutate(
      x_jitter = as.numeric(Strategy) + runif(n(), -jitter_width, jitter_width),
      y_jitter = MAE + runif(n(), -0.02, 0.02)  # Small vertical jitter
    ) %>%
    ungroup()
  
  # Prepare paired connection lines
  line_data <- jittered_points %>%
    pivot_wider(
      id_cols = Fold,
      names_from = Strategy,
      values_from = c(x_jitter, y_jitter)
    ) %>%
    drop_na()  # Remove incomplete pairs
  
  # Run statistical test (paired t-test)
  stat_test <- plot_data %>%
    t_test(MAE ~ Strategy, paired = TRUE) %>%
    add_xy_position(x = "Strategy")
  
  # Build the plot
  ggplot(plot_data, aes(x = Strategy, y = MAE)) +
    
    # Boxplot layer (semi-transparent)
    geom_boxplot(
      aes(color = Strategy),
      width = 0.4,
      alpha = 0.3,
      outlier.shape = NA
    ) +
    
    # Paired connection lines
    geom_segment(
      data = line_data,
      aes(x = x_jitter_Overall, xend = x_jitter_Segmented,
          y = y_jitter_Overall, yend = y_jitter_Segmented,
          group = Fold),
      color = "grey60",
      alpha = segment_alpha,
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    
    # Individual data points (jittered)
    geom_point(
      data = jittered_points,
      aes(x = x_jitter, y = y_jitter, color = Strategy),
      size = 1.5,
      na.rm = TRUE
    ) +
    
    # Statistical annotation
    stat_pvalue_manual(
      stat_test,
      label = "p = {p}",
      tip.length = 0.01,
      bracket.nudge.y = 0.5  # Adjust as needed
    ) +
    
    # Visual styling
    scale_color_manual(values = color_palette) +
    labs(
      title = paste(age_group, ":", model_name),
      x = "Model Strategy",
      y = "Mean Absolute Error (Years)"
    ) +
    theme_minimal(base_size = base_font_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
}

#' Perform Paired T-tests Comparing Segmented vs Overall Modeling Strategies
#' 
#' Performs paired t-tests to compare model performance between Segmented and Overall strategies.
#' Positive mean difference indicates Segmented performs better (lower MAE).
#'
#' @param result_df Dataframe containing performance metrics. Required columns:
#'   Algorithm, Age_Group, Fold, Model_Type, and metric_col.
#' @param metric_col Name of the performance metric column (default = "MAE").
#' @return A dataframe with test results containing:
#'   \itemize{
#'     \item Algorithm, Age_Group: Grouping variables
#'     \item Group_A/Group_B: Compared groups (Overall vs Segmented)
#'     \item Mean_Diff: Performance difference (Group_A - Group_B)
#'     \item CI_Lower/CI_Upper: 95% confidence intervals
#'     \item t, df, p_value: Test statistics
#'     \item Sig: Significance markers
#'   }
#' @examples
#' # Example usage:
#' test_results <- calculate_pairwise_tests(mae_data)
calculate_pairwise_tests <- function(result_df, metric_col = "MAE") {
  
  # Data preparation - convert to wide format
  result_wide <- result_df %>%
    filter(Algorithm == "SVR", Age_Group %in% c("Young", "Old")) %>%
    pivot_wider(
      names_from = Model_Type,
      values_from = all_of(metric_col),
      id_cols = c(Algorithm, Age_Group, Fold)
    )
  
  # Calculate performance difference (Overall - Segmented)
  # Positive values indicate Segmented performs better (lower MAE)
  result_diff <- result_wide %>%
    mutate(improvement = Overall - Segmented)
  
  # Perform paired t-tests for each group
  result_diff %>%
    group_by(Algorithm, Age_Group) %>%
    summarise(
      Group_A = "Overall",
      Group_B = "Segmented",
      Mean_Diff = mean(improvement),
      CI_Lower = t.test(improvement)$conf.int[1],
      CI_Upper = t.test(improvement)$conf.int[2],
      t = t.test(improvement)$statistic,
      df = t.test(improvement)$parameter,
      p_value = t.test(improvement)$p.value,
      .groups = "drop"
    ) %>%
    # Add significance markers
    mutate(
      Sig = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "NS"
      )
    )
}

