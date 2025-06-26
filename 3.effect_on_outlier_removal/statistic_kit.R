#' Perform Paired MAE Comparison Between Feature Selection Strategies
#'
#' This function conducts paired t-tests to compare Mean Absolute Error (MAE) between matched
#' "non-removal" and "removal" feature selection strategies across different models. Results
#' include confidence intervals, effect sizes, and statistical significance indicators.
#'
#' @param mae_summary A data frame containing the following columns:
#' \itemize{
#'   \item Model: Character or factor - model identifiers
#'   \item MAE: Numeric - mean absolute error values
#'   \item type: Character or factor with two levels: "removal" and "non-removal"
#' }
#'
#' @return
#' A tibble containing the following statistical results for each model:
#' \itemize{
#'   \item Model: Model identifiers (retained from input)
#'   \item Comparison: Description of paired comparison ("non-removal vs removal")
#'   \item Mean_Diff: Mean difference (non-removal MAE - removal MAE) rounded to 3 decimals
#'   \item CI_Lower/CI_Upper: 95% confidence interval bounds (rounded to 3 decimals)
#'   \item t_statistic: t-value (rounded to 3 decimals)
#'   \item df: Degrees of freedom (integer)
#'   \item p_value: Two-tailed p-value (rounded to 3 decimals)
#'   \item Significance: Asterisk notation ("***" = p < 0.001, "**" = p < 0.01, "*" = p < 0.05, "ns" = not significant)

paired_mae_test <- function(mae_summary) {
  mae_summary %>%
    group_by(Model) %>%
    summarise(
      t_test_res = list(
        t.test(
          x = MAE[type == "non-removal"], 
          y = MAE[type == "removal"], 
          paired = TRUE
        )
      ),
      Comparison = "non-removal vs removal",  
      Mean_Diff = round(t_test_res[[1]]$estimate, 3),        
      CI_Lower = round(t_test_res[[1]]$conf.int[1], 3),      
      CI_Upper = round(t_test_res[[1]]$conf.int[2], 3),     
      t_statistic = round(t_test_res[[1]]$statistic, 3),    
      df = t_test_res[[1]]$parameter,                       
      p_value = round(t_test_res[[1]]$p.value, 3),          
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    select(-t_test_res)  
}

