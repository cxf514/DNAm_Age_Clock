# Author: cxf  
# Date: 2025-06-14  
# Objective: MAE comparision of models bewtween different strategies 
rm(list=ls())  
options(stringAsFactors = F)  
.libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")  
library(dplyr)  
library(tidyverse)  
library(broom)  
library(rstatix)  
source('statistic_kit.R')
# Load non-removal data  
load('../2.models_comparision/2.1.methods_comparision-Data_Without_outlier_Removal/cpg_formal_nonrm.Rdata')  
feature_stats_raw = feature_stats  
lasso_summary_raw = lasso_summary  
rlasso_summary_raw = rlasso_summary  
mae_data_raw = mae_data  
mae_data_raw$type = 'non-removal'  

# Load removal data  
load('../2.models_comparision/2.2.methods_comparision-Cleaned_Data/cpg_formal_removal.Rdata')  
feature_stats_rm = feature_stats  
lasso_summary_rm = lasso_summary  
rlasso_summary_rm = rlasso_summary  
mae_data_rm = mae_data  
mae_data_rm$type = 'removal'  

# Combine MAE data  
mae_summary <- rbind(mae_data_raw, mae_data_rm)  

# Perform Paired MAE Comparison Between Two Strategies
result_df<-paired_mae_test(mae_summary)
write.table(result_df, file="strategy_comparision_summary.txt", sep="\t", row.names=F, col.names=T)  
