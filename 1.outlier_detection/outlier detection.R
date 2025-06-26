# Author: cxf
# Date: 2024.06.13
# Purpose: Outlier detection using Local Outlier Factor (LOF)

# Clean environment and load libraries
rm(list = ls())
options(stringsAsFactors = FALSE)
library(dplyr)        # Data manipulation
library(RIdeogram)    # Visualization
library(solitude)     # Isolation Forest (alternative to LOF)
library(DMwR2)        # LOF implementation
library(ggplot2)      # Data visualization
library(kneedle)
# Load methylation data (CpG sites)
load("../0.raw_data/cpg_info.Rdata")
# ------------------------------------------------------------
# LOF Outlier Detection
# ------------------------------------------------------------
# Calculate LOF scores for two subsets of the data
lof_scores_1 <- lofactor(combined_cpg[1:48, ], k = 10)  
lof_scores_2 <- lofactor(combined_cpg[49:96, ], k = 10) 
# Combine results into a data frame for plotting
plot_data <- data.frame(
  group = rep(c("Dataset 1", "Dataset 2"), each = 48),
  LOF_score = c(lof_scores_1, lof_scores_2)
)
# Calculate boxplot statistics (Q1, Q3, IQR, fences)
box_stats <- plot_data %>%
  group_by(group) %>%
  summarize(
    q1 = quantile(LOF_score, 0.25),
    q3 = quantile(LOF_score, 0.75),
    iqr = q3 - q1,
    lower_fence = q1 - 1.5 * iqr,
    upper_fence = q3 + 1.5 * iqr
  )
# Identify outliers (points beyond fences)
outliers <- plot_data %>%
  left_join(box_stats, by = "group") %>%
  filter(LOF_score < lower_fence | LOF_score > upper_fence)
# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------
plot <- ggplot(plot_data, aes(x = group, y = LOF_score, color = group)) +
  geom_boxplot(outlier.shape = NA) 
# Enhanced plot with:
# - Custom-sized outlier points (jittered)
# - Color scheme and theme adjustments
final_plot <- plot +
  # Add outlier points with size proportional to LOF score
  geom_point(
    data = outliers,
    aes(size = abs(LOF_score)),
    shape = 21,
    stroke = 1.1,
    position = position_jitter(width = 0.1, height = 0)
  ) +
  # Customize scaling and colors
  scale_size(range = c(1, 5)) +  # Control point size range
  scale_color_manual(values = c("#ab4a2a", "#003366")) +  # Dataset colors
  
  # Axes and reference line
  labs(x = "", y = "LOF Score") +
  geom_hline(
    yintercept = 1.5,
    linetype = "dashed",
    color = "black"
  ) +
  
  # Theme customization
  theme_minimal() +
  theme(
    title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.key.size = unit(25, "pt"),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# Save plot (90mm x 70mm is a common journal figure size)
ggsave(
  "outlier_samples.pdf",
  width = 120,
  height = 110,
  units = "mm"
)

# output outliers' sample index
plot_data[plot_data$LOF_score>1.5,]
# group LOF_score
# 91 Dataset 2  1.741449
# 93 Dataset 2  3.001136