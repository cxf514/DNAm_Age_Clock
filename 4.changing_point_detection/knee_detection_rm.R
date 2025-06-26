# Author: cxf
# Date:2024.6.24
# Aim: break point detection fro age-related DNAm changes 
rm(list=ls())
options(stringAsFactors = F)
library(dplyr)
library(kneedle)
library(tidyr)
library(ggplot2)
library(changepoint)
library(segmented)
library(patchwork)
load('../0.raw_data/cpg_info.Rdata')
source('zscore_transform.R')
exclude_samples <- c(91, 93)
retain_id<- setdiff(1:nrow(combined_cpg),exclude_samples)
combined_info<-combined_info[retain_id,]
combined_cpg<-as.data.frame(combined_cpg[retain_id,])
# load data ---------------------------------------------------------------
correlation_df <-cpg_age_correlation(combined_cpg, combined_info)
pos_sites <- correlation_df[correlation_df$cor >= 0, "sites"] 
neg_sites <- correlation_df[correlation_df$cor < 0, "sites"] 


# mode=all ----------------------------------------------------------------
result_all <- analyze_cpg_inflection(
  combined_cpg = combined_cpg,
  combined_info = combined_info,
  neg_sites = neg_sites,
  pos_sites = pos_sites,
  option = "all"
)
agg_result<-result_all$aggregated
cpg_result<-result_all$by_site

## kneedle -----------------------------------------------------------------
knee=kneedle(as.vector(agg_result$age),as.vector(agg_result$value), 
             sensitivity = 1,decreasing = FALSE,concave=FALSE)
kneedle_point=knee[1]
# [1] 50.000000  7.291073

## CPM ---------------------------------------------------------------------
cpt_var <- cpt.var(agg_result$value, method = "PELT")
cpt_index <- cpts(cpt_var)
cpt_ages <- agg_result$age[cpt_index]
# [1] 29 50

## Segmented ---------------------------------------------------------------
lm_model <- lm(value ~ age, data = agg_result)
seg_model <- segmented(
  lm_model,
  seg.Z = ~ age,
  psi = 50
)
seg_psi=seg_model$psi[,'Est.']
# [1] 50.61535
plot_change_points(
  data = agg_result,
  kneedle_point = kneedle_point,
  cpt_ages = cpt_ages,
  seg_psi = seg_model$psi[,'Est.'],
  output_file = "break_points_all.pdf"
)

plot_dnam_results(
  cpg_result = cpg_result,
  agg_result = agg_result,
  kneedle_point = kneedle_point, 
  output_file = "pos_neg_inversed.pdf")


# Pos ---------------------------------------------------------------------
pos_all <- analyze_cpg_inflection(
  combined_cpg = combined_cpg,
  combined_info = combined_info,
  neg_sites = neg_sites,
  pos_sites = pos_sites,
  option = "pos"
)
agg_result<-pos_all$aggregated
cpg_result<-pos_all$by_site


## kneedle -----------------------------------------------------------------
knee=kneedle(as.vector(agg_result$age),as.vector(agg_result$value), 
             sensitivity = 1,decreasing = FALSE,concave=FALSE)
kneedle_point=knee[1]
#[1] 53.00000  5.28407

## CPM ---------------------------------------------------------------------
cpt_var <- cpt.var(agg_result$value, method = "PELT")
cpt_index <- cpts(cpt_var)
cpt_ages <- agg_result$age[cpt_index]
#[1] 30 51

## segmented ---------------------------------------------------------------
lm_model <- lm(value ~ age, data = agg_result)
seg_model <- segmented(
  lm_model,
  seg.Z = ~ age,
  psi = 50
)
seg_psi=seg_model$psi[,'Est.']
#[1] 54.58786

plot_change_points(
  data = agg_result,
  kneedle_point = kneedle_point,
  cpt_ages = cpt_ages,
  seg_psi = seg_model$psi[,'Est.'],
  output_file = "break_points_pos.pdf"
)

plot_dnam_results(
  cpg_result = cpg_result,
  agg_result = agg_result,
  kneedle_point = kneedle_point, 
  output_file = "pos.pdf")

# Neg ---------------------------------------------------------------------
neg_all <- analyze_cpg_inflection(
  combined_cpg = combined_cpg,
  combined_info = combined_info,
  neg_sites = neg_sites,
  pos_sites = pos_sites,
  option = "neg"
)
agg_result<-neg_all$aggregated
cpg_result<-neg_all$by_site


## kneedle -----------------------------------------------------------------
knee=kneedle(as.vector(agg_result$age),as.vector(agg_result$value), 
             sensitivity = 1,decreasing = TRUE,concave=TRUE)
kneedle_point=knee[1]
#[1] 46.000000 -2.197368

## CPM ---------------------------------------------------------------------
cpt_var <- cpt.var(agg_result$value, method = "PELT")
cpt_index <- cpts(cpt_var)
cpt_ages <- agg_result$age[cpt_index]
#[1] 28 48

## segmented ---------------------------------------------------------------
lm_model <- lm(value ~ age, data = agg_result)
seg_model <- segmented(
  lm_model,
  seg.Z = ~ age,
  psi = 50
)
seg_psi=seg_model$psi[,'Est.']
#[1] 46.61851

plot_change_points(
  data = agg_result,
  kneedle_point = kneedle_point,
  cpt_ages = cpt_ages,
  seg_psi = seg_model$psi[,'Est.'],
  output_file = "break_points_neg.pdf"
)

plot_dnam_results(
  cpg_result = cpg_result,
  agg_result = agg_result,
  kneedle_point = kneedle_point, 
  output_file = "neg.pdf")