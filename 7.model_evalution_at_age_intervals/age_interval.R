# Author: cxf
# Date: 2025-06-24
# Objective: model performance at different age intervals
# Environment Setup -------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
.libPaths("/public/home/chuxufeng/R/x86_64-redhat-linux-gnu-library/3.6")
load("../0.raw_data/cpg_info.Rdata")
load(file = "../4.SVR_singleton/SVR.RData")
load('Datasets_for_reported_models.RData') # beta_matrix of CpGs in Liu's and Huang's model in dataset1 or dataset2
source("function_CpG_age.r")
library(dplyr)
library(caret)
library(glmnet)
library(e1071)
library(reshape2)

exclude_samples <- c(91, 93)
retained_id<- setdiff(1:nrow(combined_cpg), exclude_samples)
combined_cpg <- combined_cpg[retained_id,]
combined_info<- combined_info[retained_id,]
age_info<-combined_info$Age

age_dataset1<-age_info[1:48]
age_dataset2<-age_info[49:94]
pred_dataset1<-predict(model_global, combined_cpg[1:48,])
pred_dataset2<-predict(model_global, combined_cpg[49:94,])
svr.dat1<-data_frame('Age'=age_dataset1,'Pre_Age'=pred_dataset1,'model'='SVR')
svr.dat2<-data_frame('Age'=age_dataset2,'Pre_Age'=pred_dataset2,'model'='SVR')
svr.dat1$diff<-abs(svr.dat1$Age-svr.dat1$Pre_Age)
svr.dat2$diff<-abs(svr.dat2$Age-svr.dat2$Pre_Age)



# dataset1 ----------------------------------------------------------------


CpG_list = colnames(dataset1)
CpG_L = CpG_list[1:9]
CpG_H = CpG_list[10:15]

pre_age_m1m2 = t(sapply(1:length(rownames(dataset1)),function(n){
  age_m1_d1 = age.predict.m1.f(dataset1, n, CpG_L)
  age_m2_d1 = age.predict.m2.f(dataset1, n, CpG_H)
  c(age_m1_d1, age_m2_d1)
}))
colnames(pre_age_m1m2) = c("model1", "model2")
pre_age_m1m2<-as.data.frame(pre_age_m1m2)
pre_age_m1m2$Age<-age_dataset1
m1m2.dat <- melt(pre_age_m1m2, id.vars = c("Age"), measure.vars = c("model1", "model2"))
m1m2.dat$diff<-abs(m1m2.dat$Age-m1m2.dat$value)
colnames(m1m2.dat)<-c('Age','model','Pre_Age','diff')

m1m2m3.dat<-rbind(m1m2.dat,svr.dat1)

breaks <- c(-Inf, 18, 30, 40, 50, 60, 70, Inf)
labels <- c("~18", "18~30", "30~40", "40~50", "50~60", "60~70", "70~")
m1m2m3.dat$age_group <- cut(m1m2m3.dat$Age, breaks = breaks, labels = labels, include.lowest = TRUE)



# pic_dict ----------------------------------------------------------------
ggplot(data = m1m2m3.dat, 
            aes(x = age_group, y = abs(diff), colour = model)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "darkgray", lty = 2, linewidth = 0.8) +
  stat_summary(aes(group = model, shape = model), fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", position = position_dodge(width = 0.7),
               colour = "#264653") +
  ylim(-20,75)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 11),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position.inside  = c(0.1, 0.85),
        legend.title = element_blank(),
        plot.title = element_text(colour = "black", size = 13, 
                                  face = "bold", hjust = 0.5),
        legend.text = element_text(colour = "black", size = 10)) +
  scale_shape_manual(values = c(16, 15,18), labels = c("Model 1", "Model 2",'Model 3')) +
  scale_colour_manual(values = c("orange", "brown","#8d93fe"),labels = c("Model 1", "Model 2",'Model 3')) +
  labs(title = "Dataset1", x = "Age Interval", y = "|Chronological Age - Predicted Age|")
ggsave("d1_interval.tiff", dpi=500,width = 10, height = 5, units = "in")
ggsave("d1_interval.pdf",device = "pdf",width = 10,height = 5,units = "in")


# dataset2 ----------------------------------------------------------------

CpG_list = colnames(dataset2);
CpG_H = CpG_list[1:6]

pre_age_m2 = t(sapply(1:length(rownames(dataset2)),function(n){
  age_m2_d2 = age.predict.m2.f(dataset2, n, CpG_H)
}))

m2.dat<-data_frame('Age'=age_dataset2,'Pre_Age'=as.vector(pre_age_m2))

m2.dat$diff<-abs(m2.dat$Age-m2.dat$Pre_Age)
m2.dat$model<-'model2'
m2m3.dat<-rbind(m2.dat,svr.dat2)

breaks <- c(-Inf, 18, 30, 40, 50, 60, 70, Inf)
labels <- c("~18", "18~30", "30~40", "40~50", "50~60", "60~70", "70~")
m2m3.dat$age_group <- cut(m2m3.dat$Age, breaks = breaks, labels = labels, include.lowest = TRUE)

ggplot(data = m2m3.dat, 
       aes(x = age_group, y = abs(diff), 
           shape = model, colour = model)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "darkgray", lty = 2, linewidth = 0.8) +
  stat_summary(aes(group = model, shape = model), fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "pointrange", position = position_dodge(width = 0.7),
               colour = "#264653") +
  ylim(-20,75)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 11),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position.inside  = c(0.1, 0.85),
        legend.title = element_blank(),
        plot.title = element_text(colour = "black", size = 13, 
                                  face = "bold", hjust = 0.5),
        legend.text = element_text(colour = "black", size = 10)) +
  scale_shape_manual(values = c(15,18), labels = c("Model 2",'SVR')) +
  scale_colour_manual(values = c("brown","#8d93fe"),labels = c("Model 2",'SVR')) +
  labs(title = "Dataset2", x = "Age Interval", y = "|Chronological Age - Predicted Age|")
ggsave("d2_interval.tiff", dpi=500,width = 10, height = 5, units = "in")
ggsave("d2_interval.pdf",device = "pdf",width = 10,height = 5,units = "in")


# statistical_analysis ----------------------------------------------------
# dataset1 ----------------------------------------------------------------
m1m2m3.tdat<-data_frame(
                      'diff_model1'=m1m2.dat[m1m2.dat$model=='model1',]$diff,
                      'diff_model2'=m1m2.dat[m1m2.dat$model=='model2',]$diff,
                      'diff_SVR'=svr.dat1$diff,
                      Age=svr.dat1$Age)

m1m2m3.tdat <- m1m2m3.tdat %>% 
  mutate(
    paired_diff = diff_model2 - diff_SVR,  # model2与SVR的预测差异
    age_group = cut(Age, 
                    breaks = breaks,
                    labels = labels))

m1m2m3.tdat <- m1m2m3.tdat %>%
  mutate(
    model1_vs_svr = diff_model1 - diff_SVR,  # Model1与SVR的差值对比
    model2_vs_svr = diff_model2 - diff_SVR   # Model2与SVR的差值对比
  )

model1_results <- m1m2m3.tdat %>%
  group_by(age_group) %>%
  summarise(
    n = n(),
    mean_diff = mean(model1_vs_svr),
    t_test = list(t.test(diff_model1,diff_SVR,alternative = 'greater',  paired = TRUE)),
    p_value = t_test[[1]]$p.value,
    ci_lower = t_test[[1]]$conf.int[1],
    ci_upper = Inf,
    .groups = "drop"
  ) %>%
  select(-t_test)

model2_results <- m1m2m3.tdat %>%
  group_by(age_group) %>%
  summarise(
    n = n(),
    mean_diff = mean(model2_vs_svr),
    t_test = list(t.test(diff_model2,diff_SVR,alternative = 'greater', paired = TRUE)),
    p_value = t_test[[1]]$p.value,
    ci_lower = t_test[[1]]$conf.int[1],
    ci_upper = Inf,
    .groups = "drop"
  ) %>%
  select(-t_test)


final_results_d1 <- bind_rows(
  model1_results %>% mutate(comparison = "Model1 vs SVR"),
  model2_results %>% mutate(comparison = "Model2 vs SVR")
) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    ci_display = sprintf("%.2f [%.2f, +∞)", mean_diff, ci_lower)
  )
write.csv(final_results_d1, "model_performance_comparison_age_interval_dataset1.csv", row.names = FALSE)
# dataset2

m2m3.tdat<-data_frame('diff_model2'=m2.dat$diff,
                      'diff_SVR'=svr.dat2$diff,
                      Age=svr.dat2$Age)

m2m3.tdat <- m2m3.tdat %>% 
  mutate(
    paired_diff = diff_model2 - diff_SVR, 
    age_group = cut(Age, 
                    breaks = breaks,
                    labels = labels)
  )
results <- m2m3.tdat %>% 
  group_by(age_group) %>% 
  summarise(
    n = n(),
    mean_diff = mean(paired_diff),
    t_test = list(t.test(diff_model2,diff_SVR,alternative = 'greater',  paired = TRUE)),
    p_value = t_test[[1]]$p.value,
    ci_lower = t_test[[1]]$conf.int[1],
    ci_upper = Inf,
    .groups = "drop"
  ) %>%
  select(-t_test)


final_results_d2 <- results%>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    ci_display = sprintf("%.2f [%.2f, +∞)", mean_diff, ci_lower)
  )
write.csv(final_results_d2, "model_performance_comparison_age_interval_dataset2.csv", row.names = FALSE)
save(svr.dat1,svr.dat2,m1m2.dat,m2.dat,file = 'model_performance.RData')
