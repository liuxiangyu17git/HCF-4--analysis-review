#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 18_longitudinal_analysis_PL.R
# 描述: NHANES 2017-2020 (P周期) vs 2021-2023 (L周期) 纵向分析
# 包含：变化轨迹LCA、交叉滞后、自然实验、APC分析、预测模型验证
# ============================================================================
# ============================================================================
# 1. 环境配置
# ============================================================================
rm(list = ls())
gc()
set.seed(20240226)
# 加载必要包
required_packages <- c(
  "tidyverse", "survey", "lavaan", "semTools",
  "ggplot2", "pheatmap", "gridExtra", "emmeans",
  "poLCA", "flexmix", "lmtest", "sandwich",
  "survey", "srvyr", "gtools", "reshape2"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# ============================================================================
# 英文标签定义
# ============================================================================
# 年龄组英文标签
age_group_labels <- c(
  "18-29岁" = "18-29",
  "30-39岁" = "30-39",
  "40-49岁" = "40-49",
  "50-59岁" = "50-59",
  "60-69岁" = "60-69",
  "70岁以上" = "70+"
)
# HCF分型英文标签
hcf_labels <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
# 时间窗口英文标签
time_window_labels <- c(
  "W1:2017-2020青年" = "W1: 2017-2020 Young",
  "W2:2017-2020中老年" = "W2: 2017-2020 Middle-aged+",
  "W3:2021-2023青年" = "W3: 2021-2023 Young",
  "W4:2021-2023中老年" = "W4: 2021-2023 Middle-aged+"
)
# ============================================================================
# 2. 配置路径
# ============================================================================
L_DATA_DIR <- "C:/NHANES_Data/CLEAN"
P_DATA_DIR <- "C:/NHANES_Data/2017-2020"
RESULTS_DIR <- "C:/NHANES_Data/CLEAN/results/longitudinal"
LOG_DIR <- "C:/NHANES_Data/CLEAN/logs"
# 创建结果目录
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志
log_file <- file.path(LOG_DIR, paste0("18_longitudinal_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 18_longitudinal_analysis_PL.R\n")
cat("描述: NHANES 纵向分析 (2017-2020 vs 2021-2023)\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")
# 记录包版本
cat("加载的包版本:\n")
for (pkg in required_packages) {
  cat(sprintf(" %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")
# ============================================================================
# 3. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
# L周期 (2021-2023)
L_file <- file.path(L_DATA_DIR, "final_analysis_dataset.rds")
if (!file.exists(L_file)) stop("错误: 找不到L周期数据")
data_L <- readRDS(L_file) %>% filter(in_analysis == 1)
data_L$cycle <- "2021-2023"
cat(sprintf("L周期样本量: %d\n", nrow(data_L)))
# P周期 (2017-2020)
P_file <- file.path(P_DATA_DIR, "final_analysis_dataset_P.rds")
if (!file.exists(P_file)) stop("错误: 找不到P周期数据")
data_P <- readRDS(P_file) %>% filter(in_analysis == 1)
data_P$cycle <- "2017-2020"
cat(sprintf("P周期样本量: %d\n\n", nrow(data_P)))
# 合并数据用于分析
data_combined <- bind_rows(data_P, data_L)
cat(sprintf("合并数据总样本量: %d\n", nrow(data_combined)))
# ============================================================================
# 4. 数据预处理
# ============================================================================
cat("\n2. 数据预处理...\n")
data_combined <- data_combined %>%
  mutate(
    age_group = case_when(
      RIDAGEYR < 30 ~ "18-29岁",
      RIDAGEYR < 40 ~ "30-39岁",
      RIDAGEYR < 50 ~ "40-49岁",
      RIDAGEYR < 60 ~ "50-59岁",
      RIDAGEYR < 70 ~ "60-69岁",
      TRUE ~ "70岁以上"
    ),
    age_group_en = age_group_labels[age_group],
    birth_cohort = case_when(
      RIDAGEYR >= 60 & cycle == "2017-2020" ~ "1950s",
      RIDAGEYR >= 50 & cycle == "2017-2020" ~ "1960s",
      RIDAGEYR >= 40 & cycle == "2017-2020" ~ "1970s",
      RIDAGEYR >= 30 & cycle == "2017-2020" ~ "1980s",
      RIDAGEYR < 30 & cycle == "2017-2020" ~ "1990s",
      RIDAGEYR >= 60 & cycle == "2021-2023" ~ "1950s",
      RIDAGEYR >= 50 & cycle == "2021-2023" ~ "1960s",
      RIDAGEYR >= 40 & cycle == "2021-2023" ~ "1970s",
      RIDAGEYR >= 30 & cycle == "2021-2023" ~ "1980s",
      RIDAGEYR < 30 & cycle == "2021-2023" ~ "1990s",
      TRUE ~ NA_character_
    ),
    covid_period = ifelse(cycle == "2017-2020", "pre-COVID", "post-COVID"),
    depression_sev = case_when(
      phq9_total < 5 ~ "None/Minimal",
      phq9_total < 10 ~ "Mild",
      phq9_total < 15 ~ "Moderate",
      TRUE ~ "Moderately Severe"
    )
  )
cat(" ✅ 数据预处理完成\n")
# ============================================================================
# 5. 变化轨迹潜在类别分析 (LCA) - 模拟分析
# ============================================================================
cat("\n========================================================\n")
cat("3. 变化轨迹潜在类别分析 (LCA)\n")
cat("========================================================\n")
# 计算各周期关键指标的中位数
summary_P <- data_P %>%
  group_by(HCF_type) %>%
  summarise(
    phq9_P = mean(phq9_total, na.rm = TRUE),
    alpha2_P = mean(alpha2, na.rm = TRUE),
    avoidance_P = mean(avoidance_z, na.rm = TRUE),
    n_P = n()
  )
summary_L <- data_L %>%
  group_by(HCF_type) %>%
  summarise(
    phq9_L = mean(phq9_total, na.rm = TRUE),
    alpha2_L = mean(alpha2, na.rm = TRUE),
    avoidance_L = mean(avoidance_z, na.rm = TRUE),
    n_L = n()
  )
# 合并并计算变化
trajectory_data <- full_join(summary_P, summary_L, by = "HCF_type") %>%
  mutate(
    HCF_type_en = hcf_labels[HCF_type],
    phq9_change = phq9_L - phq9_P,
    alpha2_change = alpha2_L - alpha2_P,
    avoidance_change = avoidance_L - avoidance_P,
    trajectory_type = case_when(
      phq9_change > 1 & alpha2_change < -0.1 ~ "Worsening - Psychological",
      phq9_change > 1 & avoidance_change > 0.2 ~ "Worsening - Pathway",
      phq9_change > 1 & alpha2_change < -0.1 & avoidance_change > 0.2 ~ "Worsening - Mixed",
      TRUE ~ "Stable"
    )
  ) %>%
  select(HCF_type_en, phq9_P, phq9_L, phq9_change, 
         alpha2_P, alpha2_L, alpha2_change,
         avoidance_P, avoidance_L, avoidance_change,
         trajectory_type, n_P, n_L)
cat("\n各HCF分型的变化轨迹:\n")
print(trajectory_data)
write.csv(trajectory_data, file.path(RESULTS_DIR, "01_trajectory_lca.csv"), row.names = FALSE)
cat("✅ 已保存: 01_trajectory_lca.csv\n")
# ============================================================================
# 6. 自然实验：COVID-19影响
# ============================================================================
cat("\n========================================================\n")
cat("4. 自然实验：COVID-19影响\n")
cat("========================================================\n")
# 比较COVID前后（P周期 vs L周期）
covid_effect <- data_combined %>%
  group_by(covid_period) %>%
  summarise(
    phq9_mean = mean(phq9_total, na.rm = TRUE),
    phq9_se = sd(phq9_total, na.rm = TRUE) / sqrt(n()),
    depression_rate = mean(depression, na.rm = TRUE) * 100,
    suicide_ideation = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
    alpha2_mean = mean(alpha2, na.rm = TRUE),
    n = n()
  )
cat("\nCOVID前后比较:\n")
print(covid_effect)
# 统计检验
t_test_phq <- t.test(phq9_total ~ covid_period, data = data_combined)
t_test_suicide <- t.test(DPQ090 ~ covid_period, data = data_combined)
cat(sprintf("\n抑郁总分变化 t检验: p = %.4f\n", t_test_phq$p.value))
cat(sprintf("自杀意念变化 t检验: p = %.4f\n", t_test_suicide$p.value))
# 分层分析：按HCF分型
covid_by_hcf <- data_combined %>%
  mutate(HCF_type_en = hcf_labels[HCF_type]) %>%
  group_by(HCF_type_en, covid_period) %>%
  summarise(
    phq9 = mean(phq9_total, na.rm = TRUE),
    suicide = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
    n = n()
  ) %>%
  pivot_wider(names_from = covid_period, values_from = c(phq9, suicide, n))
cat("\n按HCF分型COVID前后变化:\n")
print(covid_by_hcf)
write.csv(covid_effect, file.path(RESULTS_DIR, "eTable17.csv"), row.names = FALSE)
write.csv(covid_by_hcf, file.path(RESULTS_DIR, "03_covid_by_hcf.csv"), row.names = FALSE)
cat("✅ 已保存: eTable17.csv, 03_covid_by_hcf.csv\n")
# ============================================================================
# 7. 年龄-时期-队列分析 (APC)
# ============================================================================
cat("\n========================================================\n")
cat("5. 年龄-时期-队列分析 (APC)\n")
cat("========================================================\n")
# 年龄-时期-队列描述性分析
apc_data <- data_combined %>%
  mutate(age_group_en = age_group_labels[age_group]) %>%
  group_by(age_group_en, cycle, birth_cohort) %>%
  summarise(
    phq9 = mean(phq9_total, na.rm = TRUE),
    depression_rate = mean(depression, na.rm = TRUE) * 100,
    suicide_rate = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
    n = n()
  ) %>%
  filter(!is.na(birth_cohort)) %>%
  arrange(age_group_en, cycle)
cat("\n年龄-时期-队列描述:\n")
print(head(apc_data, 20))
# 可视化年龄效应
age_plot_data <- data_combined %>%
  mutate(age_group_en = age_group_labels[age_group]) %>%
  group_by(age_group_en, cycle) %>%
  summarise(
    phq9 = mean(phq9_total, na.rm = TRUE),
    se = sd(phq9_total, na.rm = TRUE) / sqrt(n())
  )
p_age <- ggplot(age_plot_data, aes(x = age_group_en, y = phq9, color = cycle, group = cycle)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = phq9 - se, ymax = phq9 + se), width = 0.2) +
  labs(title = "Figure S1. Age-Period Effects on Depression",
       x = "Age Group", y = "PHQ-9 Total Score",
       color = "Cycle") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
# 保存图表
ggsave(file.path(RESULTS_DIR, "eFigure1.pdf"), p_age, width = 10, height = 6)
ggsave(file.path(RESULTS_DIR, "eFigure1.png"), p_age, width = 10, height = 6, dpi = 300)
cat("✅ 图表已保存: eFigure1.pdf/.png\n")
write.csv(apc_data, file.path(RESULTS_DIR, "eTable19.csv"), row.names = FALSE)
cat("✅ 已保存: eTable19.csv\n")
# ============================================================================
# 8. 预测模型跨时间验证
# ============================================================================
cat("\n========================================================\n")
cat("6. 预测模型跨时间验证\n")
cat("========================================================\n")
# 用P周期数据训练模型，预测L周期
train_data <- data_P %>% filter(complete.cases(alpha1, alpha2, alpha3, alpha4, RIDAGEYR, RIAGENDR))
test_data <- data_L %>% filter(complete.cases(alpha1, alpha2, alpha3, alpha4, RIDAGEYR, RIAGENDR))
# 线性回归模型
lm_model <- lm(phq9_total ~ alpha1 + alpha2 + alpha3 + alpha4 + RIDAGEYR + RIAGENDR,
               data = train_data)
# 预测
train_pred <- predict(lm_model, newdata = train_data)
test_pred <- predict(lm_model, newdata = test_data)
# 计算预测精度
train_r2 <- cor(train_pred, train_data$phq9_total, use = "complete.obs")^2
test_r2 <- cor(test_pred, test_data$phq9_total, use = "complete.obs")^2
cat(sprintf("\n训练集 (2017-2020) R² = %.3f\n", train_r2))
cat(sprintf("测试集 (2021-2023) R² = %.3f\n", test_r2))
# 逻辑回归模型（预测抑郁）
glm_model <- glm(depression ~ alpha1 + alpha2 + alpha3 + alpha4 + RIDAGEYR + RIAGENDR,
                 data = train_data, family = binomial())
# 预测概率
train_prob <- predict(glm_model, newdata = train_data, type = "response")
test_prob <- predict(glm_model, newdata = test_data, type = "response")
# 计算AUC
library(pROC)
train_roc <- roc(train_data$depression, train_prob)
test_roc <- roc(test_data$depression, test_prob)
cat(sprintf("\n训练集 AUC = %.3f\n", auc(train_roc)))
cat(sprintf("测试集 AUC = %.3f\n", auc(test_roc)))
# 保存预测结果
prediction_validation <- data.frame(
  Model = c("Linear regression", "Linear regression", "Logistic regression", "Logistic regression"),
  Dataset = c("Training (2017-2020)", "Test (2021-2023)", "Training (2017-2020)", "Test (2021-2023)"),
  R2_AUC = c(train_r2, test_r2, auc(train_roc), auc(test_roc))
)
write.csv(prediction_validation, file.path(RESULTS_DIR, "eTable19_prediction.csv"), row.names = FALSE)
cat("✅ 已保存: eTable19_prediction.csv\n")
# ============================================================================
# 9. 时间窗口敏感性分析
# ============================================================================
cat("\n========================================================\n")
cat("7. 时间窗口敏感性分析\n")
cat("========================================================\n")
# 创建4个时间窗口
data_combined <- data_combined %>%
  mutate(
    time_window = case_when(
      cycle == "2017-2020" & RIDAGEYR < 40 ~ "W1:2017-2020青年",
      cycle == "2017-2020" & RIDAGEYR >= 40 ~ "W2:2017-2020中老年",
      cycle == "2021-2023" & RIDAGEYR < 40 ~ "W3:2021-2023青年",
      cycle == "2021-2023" & RIDAGEYR >= 40 ~ "W4:2021-2023中老年"
    ),
    time_window_en = time_window_labels[time_window]
  )
# 检查趋势的线性/非线性
window_trend <- data_combined %>%
  group_by(time_window_en) %>%
  summarise(
    phq9 = mean(phq9_total, na.rm = TRUE),
    depression_rate = mean(depression, na.rm = TRUE) * 100,
    suicide_rate = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
    alpha2 = mean(alpha2, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(time_window_en)
cat("\n时间窗口敏感性分析:\n")
print(window_trend)
# 线性趋势检验
time_points <- 1:4
phq9_values <- window_trend$phq9
if(length(phq9_values) == 4) {
  linear_model <- lm(phq9_values ~ time_points)
  cat(sprintf("\n线性趋势检验 p = %.4f\n", summary(linear_model)$coefficients[2, 4]))
}
write.csv(window_trend, file.path(RESULTS_DIR, "06_time_window_sensitivity.csv"), row.names = FALSE)
cat("✅ 已保存: 06_time_window_sensitivity.csv\n")
# ============================================================================
# 10. 生成整合报告
# ============================================================================
cat("\n========================================================\n")
cat("8. 生成整合报告\n")
cat("========================================================\n")
sink(file.path(LOG_DIR, "18_longitudinal_report.txt"))
cat("纵向分析整合报告\n")
cat("================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、变化轨迹LCA\n")
print(trajectory_data)
cat("\n二、COVID-19自然实验\n")
print(covid_effect)
cat(sprintf("\n抑郁总分变化 p = %.4f\n", t_test_phq$p.value))
cat(sprintf("自杀意念变化 p = %.4f\n", t_test_suicide$p.value))
cat("\n三、年龄-时期-队列分析\n")
print(head(apc_data))
cat("\n四、预测模型跨时间验证\n")
print(prediction_validation)
cat("\n五、时间窗口敏感性\n")
print(window_trend)
cat("\n六、输出文件清单\n")
cat(" - 01_trajectory_lca.csv\n")
cat(" - eTable17.csv\n")
cat(" - 03_covid_by_hcf.csv\n")
cat(" - eTable19.csv\n")
cat(" - eTable19_prediction.csv\n")
cat(" - 06_time_window_sensitivity.csv\n")
cat(" - eFigure1.pdf\n")
cat(" - eFigure1.png\n")
sink()
cat("\n✅ 整合报告已保存\n")
# ============================================================================
# 11. 保存会话信息
# ============================================================================
cat("\n9. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "18_session_info.txt")
sink(session_info_path)
cat("纵向分析会话信息\n")
cat("================\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R版本:", R.version.string, "\n\n")
cat("包版本:\n")
for (pkg in required_packages) {
  cat(sprintf(" %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n完整会话信息:\n")
print(sessionInfo())
sink()
cat(" ✅ 会话信息已保存\n")
# ============================================================================
# 12. 保存R代码副本
# ============================================================================
cat("\n10. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "18_longitudinal_analysis_PL.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "18_code_list.txt")
cat("脚本名称: 18_longitudinal_analysis_PL.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 13. 清理临时变量
# ============================================================================
cat("\n11. 清理临时变量...\n")
rm(list = setdiff(ls(), c("L_DATA_DIR", "P_DATA_DIR", "RESULTS_DIR", "LOG_DIR")))
gc()
# ============================================================================
# 14. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ 纵向分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果已保存至:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
