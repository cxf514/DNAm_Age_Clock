#' @param combined_cpg Dataframe of CpG methylation values (sites in columns)
#' @param combined_info Metadata dataframe containing chronological age (must include "Age" column)
#' @param neg_sites Character vector of negatively correlated sites (for option selection)
#' @param pos_sites Character vector of positively correlated sites (for option selection)
#' @param option Analysis mode: "all" (default), "neg" (negative-correlation only), or "pos" (positive-correlation only)
#' @return Dataframe with aggregated methylation values across specified sites by age
#' 
#' @details
#' This function performs the following analytical pipeline:
#'   1. Z-score normalization of all CpG sites
#'   2. Option-based site selection (all/negative-correlation/positive-correlation)
#'   3. Cubic polynomial fitting for each selected site
#'   4. Aggregation of predicted trajectories across sites
#' 
#' For "all" mode, negatively correlated sites are sign-flipped before aggregation.
#' @examples
#' result_all<- analyze_cpg_inflection(
#'   combined_cpg = combined_cpg,
#'   combined_info = combined_info,
#'   neg_sites = c("cg123","cg456"),
#'   pos_sites = c("cg789","cg101"),
#'   option = "all"
#' )




analyze_cpg_inflection <- function(combined_cpg, combined_info, neg_sites = NULL, pos_sites = NULL, option = "all") {
  
  # --------------------------
  # Input validation
  # --------------------------
  if (!"Age" %in% colnames(combined_info)) {
    stop("Input metadata must contain 'Age' column")
  }
  if (!option %in% c("all", "neg", "pos")) {
    stop("Mode must be either 'all', 'neg', or 'pos'")
  }
  
  # --------------------------
  # Data normalization
  # --------------------------
  # Column-wise Z-score standardization
  df_zscore <- combined_cpg %>%
    mutate(across(everything(), scale))
  
  # --------------------------
  # Option-specific processing
  # --------------------------
  if (option == "neg") {
    # Negative-correlation mode: select only specified negative sites
    if (is.null(neg_sites)) {
      stop("Negative-correlation mode requires neg_sites parameter")
    }
    df_zscore <- df_zscore %>% 
      dplyr::select(all_of(neg_sites))
    
  } else if (option == "pos") {
    # Positive-correlation mode: select only specified positive sites
    if (is.null(pos_sites)) {
      stop("Positive-correlation mode requires pos_sites parameter")
    }
    df_zscore <- df_zscore %>% 
      dplyr::select(all_of(pos_sites))
    
  } else {
    # Default mode: sign-flip negative-correlation sites
    if (!is.null(neg_sites)) {
      df_zscore <- df_zscore %>%
        dplyr::mutate(across(any_of(neg_sites), ~ -.))
    }
  }
  
  # --------------------------
  # Model preparation
  # --------------------------
  # Merge with age information
  df_zscore_all <- df_zscore %>%
    dplyr::mutate(age = combined_info$Age)
  
  # Generate prediction grid
  predictions_df <- data.frame(
    age = seq(
      from = min(combined_info$Age),
      to = max(combined_info$Age),
      length.out = max(combined_info$Age) - min(combined_info$Age) + 1
    )
  )
  
  # --------------------------
  # Modeling pipeline
  # --------------------------
  # Fit cubic polynomial for each CpG site
  for (site in colnames(df_zscore)) {
    # Model formula: cubic polynomial with raw coefficients
    formula <- as.formula(paste0(site, " ~ poly(age, 3, raw=TRUE)"))
    
    # Fit linear model
    model <- lm(formula, data = df_zscore_all)
    
    # Store predictions
    predictions_df[, site] <- predict(model, newdata = predictions_df)
  }
  
  # --------------------------
  # Result aggregation
  # --------------------------
  # Convert to long format for aggregation
  long_df <- predictions_df %>% 
    pivot_longer(
      cols = -age,
      names_to = "site",
      values_to = "value"
    )
  
  # Sum values across all sites at each age point
  sum_df <- long_df %>%
    group_by(age) %>%
    summarise(value = sum(value)) %>%
    dplyr::mutate(site = "aggregated")
  
  # --------------------------
  # Output
  # --------------------------
  return(list(
    aggregated = sum_df,
    by_site = long_df))
}

#' Calculate Correlation Between CpG Sites and Age
#' Simplified version that matches your original code structure.
#' @param cpg_data Dataframe/matrix of methylation data (CpGs in columns)
#' @param sample_info Dataframe containing an "Age" column
#' @return Named list with: site names, correlation coefficients, and p-values
#'
#' @examples
#' results <- cpg_age_correlation(combined_cpg, combined_info)
cpg_age_correlation <- function(cpg_data, sample_info) {
  correlation_pvalues <- list()
  
  for (col in colnames(cpg_data)) {
    test_result <- cor.test(sample_info$Age, cpg_data[, col])
    correlation_pvalues[[col]] <- list(
      sites = col,
      cor = test_result$estimate,
      P = test_result$p.value
    )
  }
  
  return(do.call(rbind.data.frame, correlation_pvalues))
}



#' Visualize change points from multiple detection methods
#'
#' @param data Dataframe with columns `age` (x-axis) and `value` (y-axis)
#' @param kneedle_point Numeric or NULL. Detected knee point from Kneedle method.
#' @param cpt_ages Numeric vector or NULL. Change points from CPM (PELT) method.
#' @param seg_psi Numeric or NULL. Breakpoint from Segmented Regression (ignored if `seg_model` is provided).
#' @param seg_model Fitted `segmented` model object (optional). If provided, overrides `seg_psi` and plots the fitted curve.
#' @param output_file Character or NULL. File path to save plot as PDF (e.g., "output.pdf"). If NULL, plot displays in active device.
#'
#' @return No return value. Generates a plot.
#'
#' @examples
#' # Plot
#' plot_change_points(
#'   data = dat,
#'   kneedle_point = 55,
#'   cpt_ages = cpt_ages,
#'   seg_psi = 50.6,
#'   output_file = "change_points.pdf"
#' )

plot_change_points <- function(data, kneedle_point = NULL, cpt_ages = NULL, seg_psi = NULL, 
                               seg_model = NULL, output_file = NULL) {
  
  # Validate input data
  if (!all(c("age", "value") %in% colnames(data))) {
    stop("Input data must contain columns 'age' and 'value'.")
  }
  
  # Open PDF device if output file specified
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 6)
    on.exit(dev.off())  # Ensure device closure on exit
  }
  
  # Set plot margins and layout
  par(mar = c(5, 5, 4, 2), mgp = c(3, 0.8, 0))
  
  # Base scatter plot
  plot(data$age, data$value, 
       pch = 16, col = "gray50", cex = 0.8,
       main = "Change Point Detection Comparison",
       xlab = "Age (years)", ylab = "Measurement Value")
  
  # Add Segmented Regression fit (if model provided)
  if (!is.null(seg_model)) {
    plot(seg_model, add = TRUE, col = "red", lwd = 2, rug = FALSE)
    seg_psi <- seg_model$psi[, "Est."]  # Extract breakpoints from model
  }
  
  # Plot Segmented breakpoints
  if (!is.null(seg_psi)) {
    abline(v = seg_psi, lty = 2, col = "red", lwd = 1.5)
    text(x = seg_psi, y = max(data$value) * 0.95, 
         labels = paste0(round(seg_psi, 1), "y (Seg)"), 
         pos = 4, col = "red", cex = 0.8)
  }
  
  # Plot CPM breakpoints
  if (!is.null(cpt_ages)) {
    abline(v = cpt_ages, lty = 3, col = "blue", lwd = 1.5)
    text(x = cpt_ages, y = max(data$value) * 0.85, 
         labels = paste0(round(cpt_ages, 1), "y (CPM)"), 
         pos = 4, col = "blue", cex = 0.8)
  }
  
  # Plot Kneedle breakpoint
  if (!is.null(kneedle_point)) {
    abline(v = kneedle_point, lty = 4, col = "darkgreen", lwd = 1.5)
    text(x = kneedle_point, y = max(data$value) * 0.75, 
         labels = paste0(round(kneedle_point, 1), "y (Kneedle)"), 
         pos = 4, col = "darkgreen", cex = 0.8)
  }
  
  # Dynamically generate legend entries
  legend_items <- list(
    "Data" = list(pch = 16, col = "gray50", lty = NA)
  )
  
  if (!is.null(seg_psi)) {
    legend_items$`Segmented Fit` <- list(col = "red", lty = 1, pch = NA)
    legend_items$`Breakpoint (Seg)` <- list(col = "red", lty = 2, pch = NA)
  }
  if (!is.null(cpt_ages)) {
    legend_items$`Breakpoint (CPM)` <- list(col = "blue", lty = 3, pch = NA)
  }
  if (!is.null(kneedle_point)) {
    legend_items$`Breakpoint (Kneedle)` <- list(col = "darkgreen", lty = 4, pch = NA)
  }
  
  # Draw legend
  legend("bottomright", 
         legend = names(legend_items),
         col = sapply(legend_items, `[[`, "col"),
         pch = sapply(legend_items, `[[`, "pch"),
         lty = sapply(legend_items, `[[`, "lty"),
         lwd = 1.5, cex = 0.8, bg = "white")
}

#' @param cpg_result     Dataframe: Long-format DNAm z-scores by CpG site (columns: `age`, `value`, `site`).
#' @param agg_result     Dataframe: Aggregated DNAm z-scores (columns: `age`, `value`).
#' @param kneedle_point  Numeric: Age value where the breakpoint occurs.
#' @param output_file    Character: Output filename/path (default: "break_points.pdf").
#' @param panel_height   Numeric vector: Relative heights of top/bottom panels (default: c(1, 1.5)).
#' @param width          Numeric: Plot width in inches (default: 8).
#' @param height         Numeric: Plot height in inches (default: 6).
#'
#' @return A combined ggplot object (invisibly). Output is saved to `output_file`.
#' @example
#' plot_dnam_results(
#'   cpg_result = your_cpg_data,
#'   agg_result = your_agg_data,
#'   kneedle_point = 42.5,  # Specific age for breakpoint
#'   output_file = "custom_breakpoint.pdf"
plot_dnam_results <- function(cpg_result, agg_result, kneedle_point, 
                              output_file = "break_points.pdf",
                              panel_height = c(1, 1), width = 4, height = 6) {
  
  # --- Panel 1: Individual CpG Site Trajectories ---
  p1 <- ggplot(cpg_result, aes(x = age, y = value, group = site)) +
    geom_line(size = 0.8, alpha = 0.5, color = "black") +
    labs(x = "", y = "DNAm levels (z-score)") +
    theme_bw(base_rect_size = 1) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 5, b = 0, l = 5, r = 5)
    )
  
  # --- Panel 2: Aggregated Profile with Breakpoint Annotation ---
  p2 <- ggplot(agg_result, aes(x = age, y = value)) +
    geom_line(size = 1, color = "black") +
    # Add vertical line at knee point
    geom_vline(
      xintercept = kneedle_point,
      linetype = "dashed",
      color = "#ffbb78",  # Orange color for visibility
      size = 1,
      alpha = 0.8
    ) +
    # Add knee point annotation (using annotate() to avoid warning)
    annotate(
      "text",
      x = kneedle_point* 0.65,
      y = max(agg_result$value) * 0.95,  # Position at 95% of y-axis
      label = paste0('knee point=',round(kneedle_point, 1), " years"),  # Clean label format
      color = "#ffbb78",
      size = 3.5,
      hjust = -0.1  # Right adjustment
    ) +
    labs(x = "Age (years)", y = "Sum of DNAm (z-score)") +
    theme_bw(base_rect_size = 1) +
    theme(plot.margin = margin(t = 0, b = 5, l = 5, r = 5))
  
  # --- Combine Panels using Patchwork ---
  combined <- p1 / p2 + 
    plot_layout(heights = panel_height)
  
  # --- Save to PDF ---
  ggsave(
    filename = output_file,
    plot = combined,
    width = width,
    height = height,
    device = "pdf"
  )
  
  # Return plot object silently
  invisible(combined)
}



