#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 12_paper2_pathway_analysis_L.R
# 描述: NHANES 2021-2023 (L周期) 论文2 - 四通路K-means聚类分析
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q2, 假设 H4-H5
# 对应研究计划: 第五部分 图2、表2
# 对应变量详表: 第五部分 5.2 pathway_proxy_vars.rds, 5.3 pathway_clustering.rds
#
# 伦理声明: 使用NHANES公开数据，无需IRB批准
# 数据来源: https://www.cdc.gov/nchs/nhanes/index.htm
# 数据周期: 2021-2023 (L系列)
#
# 随机种子: 20240226 (固定以确保可重复性)
# 最后修改: 2026-02-20
# ============================================================================
# ============================================================================
# 1. 环境配置
# ============================================================================
rm(list = ls())
gc()
# 设置随机种子（期刊要求）
set.seed(20240226)
# 加载必要包
required_packages <- c(
  "tidyverse", "survey", "cluster", "ggplot2", 
  "pheatmap", "scales", "gridExtra"
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
# 四通路聚类英文标签
pathway_cluster_labels <- c(
  "1" = "Low-hyperactivation Healthy",
  "2" = "Aversion-Exhaustion",
  "3" = "Perseveration-dominant",
  "4" = "Pure Hyperactivation"
)
# 四通路维度英文标签
dimension_labels <- c(
  "avoidance" = "Avoidance",
  "perseveration" = "Perseveration",
  "hyperactivation" = "Hyperactivation",
  "exhaustion" = "Exhaustion"
)
# ============================================================================
# 2. 配置路径
# ============================================================================
clean_dir <- "C:/NHANES_Data/CLEAN"
results_dir <- file.path(clean_dir, "results", "paper2")
LOG_DIR <- file.path(clean_dir, "logs")
# 创建结果目录
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("12_paper2_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 12_paper2_pathway_analysis_L.R\n")
cat("描述: NHANES 2021-2023 (L周期) 论文2 - 四通路聚类分析\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R版本:", R.version.string, "\n")
cat("随机种子: 20240226\n")
cat("========================================================\n\n")
# 记录包版本（期刊要求）
cat("加载的包版本:\n")
for (pkg in required_packages) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")
# ============================================================================
# 3. 数据加载
# ============================================================================
cat("1. 加载数据...\n")
data_path <- file.path(clean_dir, "analysis_dataset_subset.rds")
if(!file.exists(data_path)) {
  stop("错误: analysis_dataset_subset.rds不存在!")
}
final_data <- readRDS(data_path)
cat(sprintf("数据加载成功！\n数据集维度: %d行 x %d列\n\n", 
            nrow(final_data), ncol(final_data)))
# ============================================================================
# 4. 定义高危人群（身心混合型 + 纯心理型）
# ============================================================================
cat("2. 定义高危人群...\n")
cat("HCF_type原始分布:\n")
print(table(final_data$HCF_type, useNA = "ifany"))
# 根据您的数据选择高危人群
high_risk_data <- final_data %>%
  filter(
    !is.na(WTINT2YR) & WTINT2YR > 0,
    HCF_type %in% c("身心混合型", "纯心理型")
  )
cat(sprintf("\n高危人群样本量: %d\n", nrow(high_risk_data)))
cat(sprintf("身心混合型: %d\n", sum(high_risk_data$HCF_type == "身心混合型")))
cat(sprintf("纯心理型: %d\n\n", sum(high_risk_data$HCF_type == "纯心理型")))
# ============================================================================
# 5. 准备四通路变量
# ============================================================================
cat("3. 准备四通路变量...\n")
pathway_vars <- c(
  "avoidance_z",      # 回避倾向
  "perseveration_z",  # 固着倾向
  "hyperactivation_z", # 过激状态
  "exhaustion_z"      # 枯竭状态
)
# 检查变量是否存在
for (var in pathway_vars) {
  if (!var %in% names(high_risk_data)) {
    stop("变量不存在: ", var)
  }
}
# 检查缺失值
missing_summary <- data.frame(
  Variable = pathway_vars,
  Missing = sapply(pathway_vars, function(v) sum(is.na(high_risk_data[[v]]))),
  Percent = NA
)
missing_summary$Percent <- round(missing_summary$Missing / nrow(high_risk_data) * 100, 2)
cat("\n四通路变量缺失情况:\n")
print(missing_summary)
# 处理缺失值 - 用中位数填充
for (var in pathway_vars) {
  if (missing_summary$Missing[missing_summary$Variable == var] > 0) {
    fill_val <- median(high_risk_data[[var]], na.rm = TRUE)
    high_risk_data[[var]][is.na(high_risk_data[[var]])] <- fill_val
    cat(sprintf("填充 %s 缺失值: %d个，使用中位数: %.3f\n", 
                var, missing_summary$Missing[missing_summary$Variable == var], fill_val))
  }
}
# 创建聚类数据矩阵
cluster_matrix <- as.matrix(high_risk_data[, pathway_vars])
cat(" ✅ 四通路变量准备完成\n\n")
# ============================================================================
# 6. K-means聚类分析
# ============================================================================
cat("4. K-means聚类分析 (k=4)...\n")
set.seed(20240226)
km_result <- kmeans(cluster_matrix, centers = 4, nstart = 50, iter.max = 100)
# 添加聚类结果
high_risk_data$cluster_raw <- km_result$cluster
# ============================================================================
# 7. 计算聚类特征并命名
# ============================================================================
cat("\n5. 计算聚类特征并命名...\n")
cluster_profiles <- high_risk_data %>%
  group_by(cluster_raw) %>%
  summarise(
    avoidance = mean(avoidance_z, na.rm = TRUE),
    perseveration = mean(perseveration_z, na.rm = TRUE),
    hyperactivation = mean(hyperactivation_z, na.rm = TRUE),
    exhaustion = mean(exhaustion_z, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(cluster_raw)
cat("\n【聚类原始均值】\n")
print(cluster_profiles)
# 根据实际数据特征命名（中文，用于内部处理）
cluster_names <- c(
  "1" = "低过激健康型",
  "2" = "对抗-枯竭混合型",
  "3" = "痴固着主导型",
  "4" = "纯生理过激型"
)
# 应用命名（中文，用于内部处理）
high_risk_data$pathway_cluster <- factor(
  high_risk_data$cluster_raw,
  levels = 1:4,
  labels = cluster_names[as.character(1:4)]
)
# 创建英文因子（用于图表输出）
high_risk_data$pathway_cluster_en <- factor(
  high_risk_data$cluster_raw,
  levels = 1:4,
  labels = pathway_cluster_labels[as.character(1:4)]
)
# 保存命名映射（英文版）
name_mapping <- data.frame(
  Cluster = 1:4,
  Cluster_Name = pathway_cluster_labels[as.character(1:4)],
  Avoidance_Mean = cluster_profiles$avoidance,
  Perseveration_Mean = cluster_profiles$perseveration,
  Hyperactivation_Mean = cluster_profiles$hyperactivation,
  Exhaustion_Mean = cluster_profiles$exhaustion,
  N = cluster_profiles$n
)
write.csv(name_mapping, file.path(results_dir, "cluster_naming.csv"), row.names = FALSE)
cat(" ✅ 聚类命名已保存\n\n")
# ============================================================================
# 8. 可视化1：四通路剖面图 (Figure 1)
# ============================================================================
cat("6. 生成可视化结果...\n")
# 准备绘图数据
plot_data <- high_risk_data %>%
  group_by(pathway_cluster) %>%
  summarise(
    avoidance = mean(avoidance_z, na.rm = TRUE),
    perseveration = mean(perseveration_z, na.rm = TRUE),
    hyperactivation = mean(hyperactivation_z, na.rm = TRUE),
    exhaustion = mean(exhaustion_z, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = -pathway_cluster,
               names_to = "dimension",
               values_to = "value") %>%
  mutate(
    dimension_en = dimension_labels[dimension]
  )
# 合并英文聚类标签
plot_data <- plot_data %>%
  left_join(
    high_risk_data %>% 
      select(pathway_cluster, pathway_cluster_en) %>% 
      distinct(),
    by = "pathway_cluster"
  )
p1 <- ggplot(plot_data, aes(x = dimension_en, y = value, fill = pathway_cluster_en)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Figure 1. Four Psychophysiological Pathways",
       x = "Dimension", y = "Z-score",
       fill = "Pathway Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  guides(fill = guide_legend(nrow = 2))
# 保存图表
ggsave(file.path(results_dir, "Figure1.pdf"), p1, width = 10, height = 6)
ggsave(file.path(results_dir, "Figure1.png"), p1, width = 10, height = 6, dpi = 300)
cat(" ✅ Figure 1. Four Psychophysiological Pathways saved\n")
# ============================================================================
# 9. 可视化2：特征热图 (Figure 2)
# ============================================================================
heatmap_data <- high_risk_data %>%
  group_by(pathway_cluster) %>%
  summarise(
    Age = mean(RIDAGEYR, na.rm = TRUE),
    Male_Pct = mean(RIAGENDR == 1, na.rm = TRUE) * 100,
    College_Pct = mean(DMDEDUC2 >= 4, na.rm = TRUE) * 100,
    Poverty_Pct = mean(INDFMPIR < 1, na.rm = TRUE) * 100,
    PHQ9 = mean(phq9_total, na.rm = TRUE),
    Depression_Pct = mean(phq9_total >= 10, na.rm = TRUE) * 100,
    BMI = mean(BMXBMI, na.rm = TRUE),
    CRP = mean(hs_crp_mgl, na.rm = TRUE),
    HR = mean(BPXOPLS1, na.rm = TRUE),
    Avoidance = mean(avoidance_z, na.rm = TRUE),
    Perseveration = mean(perseveration_z, na.rm = TRUE),
    Hyperactivation = mean(hyperactivation_z, na.rm = TRUE),
    Exhaustion = mean(exhaustion_z, na.rm = TRUE),
    n = n(),
    cluster_raw = first(cluster_raw)
  )
heat_matrix <- as.matrix(heatmap_data[, c("Age", "Male_Pct", "College_Pct", "Poverty_Pct",
                                          "PHQ9", "Depression_Pct", "BMI", "CRP", "HR",
                                          "Avoidance", "Perseveration", "Hyperactivation", "Exhaustion")])
# 用英文标签作为行名
rownames(heat_matrix) <- pathway_cluster_labels[as.character(heatmap_data$cluster_raw)]
heat_matrix_scaled <- scale(heat_matrix)
# 保存热图
pdf(file.path(results_dir, "figure2_cluster_heatmap.pdf"), width = 14, height = 7)
pheatmap(
  heat_matrix_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Figure 2. Pathway Cluster Characteristics Heatmap",
  fontsize_row = 11,
  fontsize_col = 10,
  angle_col = 45,
  display_numbers = TRUE,
  number_format = "%.1f",
  border_color = NA,
  clustering_method = "ward.D2"
)
dev.off()
cat(" ✅ Figure 2. Pathway Cluster Characteristics Heatmap saved\n\n")
# ============================================================================
# 10. 聚类特征表（加权）- Table 1
# ============================================================================
cat("7. 生成聚类特征表...\n")
design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT2YR,
  nest = TRUE,
  data = high_risk_data
)
# 获取英文聚类标签
cluster_names_levels <- levels(high_risk_data$pathway_cluster_en)
n_clusters <- length(cluster_names_levels)
# 定义变量名称（英文）
var_names <- c(
  "Variable",
  pathway_cluster_labels[as.character(1:n_clusters)]
)
# 创建结果列表
result_list <- list()
# 添加变量名行
result_list[["Variable"]] <- c(
  "N",
  "Age (years), Mean (SE)",
  "Male, % (SE)",
  "College or above, % (SE)",
  "Poverty (income<1.0), % (SE)",
  "PHQ-9 total, Mean (SE)",
  "Depression (PHQ-9≥10), % (SE)",
  "BMI (kg/m²), Mean (SE)",
  "hs-CRP (mg/L), Mean (SE)",
  "Resting heart rate (bpm), Mean (SE)"
)
# 对每个聚类计算统计量
for (i in 1:n_clusters) {
  cluster_num <- i  # 直接用数字
  cluster_name_en <- pathway_cluster_labels[as.character(cluster_num)]
  # 用 cluster_raw 进行子集选择
  design_sub <- subset(design, cluster_raw == cluster_num)
  cluster_values <- character(10)
  # 1. 样本量
  cluster_values[1] <- as.character(nrow(design_sub$variables))
  # 2. 年龄
  age_mean <- svymean(~RIDAGEYR, design_sub, na.rm = TRUE)
  cluster_values[2] <- sprintf("%.1f (%.2f)", coef(age_mean), SE(age_mean))
  # 3. 男性比例
  male_prop <- svymean(~I(RIAGENDR==1), design_sub, na.rm = TRUE)
  cluster_values[3] <- sprintf("%.1f (%.2f)", 100 * coef(male_prop), 100 * SE(male_prop))
  # 4. 大学及以上
  college_prop <- svymean(~I(DMDEDUC2>=4), design_sub, na.rm = TRUE)
  cluster_values[4] <- sprintf("%.1f (%.2f)", 100 * coef(college_prop), 100 * SE(college_prop))
  # 5. 贫困
  poverty_prop <- svymean(~I(INDFMPIR<1), design_sub, na.rm = TRUE)
  cluster_values[5] <- sprintf("%.1f (%.2f)", 100 * coef(poverty_prop), 100 * SE(poverty_prop))
  # 6. PHQ-9总分
  phq9_mean <- svymean(~phq9_total, design_sub, na.rm = TRUE)
  cluster_values[6] <- sprintf("%.1f (%.2f)", coef(phq9_mean), SE(phq9_mean))
  # 7. 抑郁比例
  dep_prop <- svymean(~I(phq9_total>=10), design_sub, na.rm = TRUE)
  cluster_values[7] <- sprintf("%.1f (%.2f)", 100 * coef(dep_prop), 100 * SE(dep_prop))
  # 8. BMI
  bmi_mean <- svymean(~BMXBMI, design_sub, na.rm = TRUE)
  cluster_values[8] <- sprintf("%.1f (%.2f)", coef(bmi_mean), SE(bmi_mean))
  # 9. hs-CRP
  crp_mean <- svymean(~hs_crp_mgl, design_sub, na.rm = TRUE)
  cluster_values[9] <- sprintf("%.1f (%.2f)", coef(crp_mean), SE(crp_mean))
  # 10. 静息心率
  hr_mean <- svymean(~BPXOPLS1, design_sub, na.rm = TRUE)
  cluster_values[10] <- sprintf("%.1f (%.2f)", coef(hr_mean), SE(hr_mean))
  result_list[[cluster_name_en]] <- cluster_values
}
# 转换为数据框
characteristics <- as.data.frame(result_list, stringsAsFactors = FALSE)
# 保存
write.csv(characteristics, file.path(results_dir, "Table1.csv"), row.names = FALSE)
cat(" ✅ 聚类特征表已保存\n\n")
# ============================================================================
# 11. 四通路扩展分析
# ============================================================================
cat("8. 四通路扩展分析...\n")
# 11.1 通路演进风险评分
if ("alpha3" %in% names(high_risk_data)) {
  high_risk_data <- high_risk_data %>%
    mutate(
      alpha3_scaled = (alpha3 - min(alpha3, na.rm = TRUE)) /
                      (max(alpha3, na.rm = TRUE) - min(alpha3, na.rm = TRUE)),
      risk_progression_score = 0.3 * avoidance_z +
                                0.3 * hyperactivation_z +
                                0.4 * (1 - alpha3_scaled),
      risk_category = case_when(
        risk_progression_score > quantile(risk_progression_score, 0.75, na.rm = TRUE) ~ "High",
        risk_progression_score > quantile(risk_progression_score, 0.25, na.rm = TRUE) ~ "Medium",
        TRUE ~ "Low"
      )
    )
} else {
  high_risk_data <- high_risk_data %>%
    mutate(
      risk_progression_score = 0.5 * avoidance_z + 0.5 * hyperactivation_z,
      risk_category = case_when(
        risk_progression_score > quantile(risk_progression_score, 0.75, na.rm = TRUE) ~ "High",
        risk_progression_score > quantile(risk_progression_score, 0.25, na.rm = TRUE) ~ "Medium",
        TRUE ~ "Low"
      )
    )
}
risk_by_cluster <- high_risk_data %>%
  group_by(pathway_cluster_en) %>%
  summarise(
    Risk_Mean = mean(risk_progression_score, na.rm = TRUE),
    Risk_SD = sd(risk_progression_score, na.rm = TRUE),
    N = n()
  )
cat("\n各聚类风险评分:\n")
print(risk_by_cluster)
write.csv(risk_by_cluster, file.path(results_dir, "risk_score_by_cluster.csv"), row.names = FALSE)
# 11.2 通路与临床结局的关联
outcome_vars <- c("phq9_total", "BMXBMI", "hs_crp_mgl", "BPXOPLS1")
outcome_names <- c("PHQ-9", "BMI", "hs-CRP", "HR")
outcome_analysis <- data.frame()
for (j in 1:length(outcome_vars)) {
  outcome <- outcome_vars[j]
  outcome_name <- outcome_names[j]
  for (i in 1:n_clusters) {
    cluster_num <- i
    cluster_name_en <- pathway_cluster_labels[as.character(cluster_num)]
    # 用 cluster_raw 进行子集选择
    cluster_data <- subset(high_risk_data, cluster_raw == cluster_num)
    mean_val <- mean(cluster_data[[outcome]], na.rm = TRUE)
    sd_val <- sd(cluster_data[[outcome]], na.rm = TRUE)
    outcome_analysis <- rbind(outcome_analysis, data.frame(
      Outcome = outcome_name,
      Pathway = cluster_name_en,
      Mean = round(mean_val, 2),
      SD = round(sd_val, 2),
      N = nrow(cluster_data)
    ))
  }
}
write.csv(outcome_analysis, file.path(results_dir, "outcome_by_cluster.csv"), row.names = FALSE)
cat(" ✅ 扩展分析完成\n\n")
# ============================================================================
# 12. 保存所有结果
# ============================================================================
cat("9. 保存结果...\n")
save_vars <- c(
  "SEQN", "WTINT2YR", "SDMVPSU", "SDMVSTRA", "HCF_type",
  "pathway_cluster", "pathway_cluster_en", "cluster_raw", pathway_vars,
  "risk_progression_score", "risk_category",
  "RIDAGEYR", "RIAGENDR", "DMDEDUC2", "INDFMPIR",
  "phq9_total", "BMXBMI", "hs_crp_mgl", "BPXOPLS1"
)
save_vars <- save_vars[save_vars %in% names(high_risk_data)]
high_risk_final <- high_risk_data[, save_vars]
saveRDS(high_risk_final, file.path(results_dir, "paper2_analysis_data.rds"))
cat(" ✅ 分析数据已保存: paper2_analysis_data.rds\n")
save(
  km_result,
  cluster_profiles,
  name_mapping,
  risk_by_cluster,
  outcome_analysis,
  characteristics,
  file = file.path(results_dir, "paper2_results.RData")
)
cat(" ✅ 聚类结果已保存: paper2_results.RData\n\n")
# ============================================================================
# 13. 生成分析报告
# ============================================================================
cat("10. 生成分析报告...\n")
report_file <- file.path(LOG_DIR, "12_paper2_report_L.txt")
sink(report_file)
cat("Paper 2: Pathway Clustering Analysis Report\n")
cat("==========================================\n\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、Sample Information\n")
cat(sprintf("  High-risk sample size: %d\n", nrow(high_risk_data)))
cat(sprintf("  Psychosomatic Mixed: %d\n", sum(high_risk_data$HCF_type == "身心混合型")))
cat(sprintf("  Psychological: %d\n\n", sum(high_risk_data$HCF_type == "纯心理型")))
cat("二、Cluster Distribution\n")
print(table(high_risk_data$pathway_cluster_en))
cat("\n三、Cluster Quality Metrics\n")
cat(sprintf("  Total SS: %.2f\n", km_result$totss))
cat(sprintf("  Within SS: %.2f\n", km_result$tot.withinss))
cat(sprintf("  Between SS: %.2f\n", km_result$betweenss))
cat(sprintf("  Variance explained: %.2f%%\n", km_result$betweenss / km_result$totss * 100))
cat("\n四、Output Files\n")
cat("  📊 Figure 1: Figure1.pdf\n")
cat("  📊 Figure 2: figure2_cluster_heatmap.pdf\n")
cat("  📋 Table 1: Table1.csv\n")
cat("  📋 Table 2: cluster_naming.csv\n")
cat("  📋 Table 3: risk_score_by_cluster.csv\n")
cat("  📋 Table 4: outcome_by_cluster.csv\n")
cat("  💾 Data: paper2_analysis_data.rds\n")
cat("  💾 Results: paper2_results.RData\n")
sink()
cat(" ✅ 分析报告已保存\n\n")
# ============================================================================
# 14. 保存会话信息（期刊要求）
# ============================================================================
cat("11. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "12_session_info_L.txt")
sink(session_info_path)
cat("NHANES Paper 2 Analysis Session Information\n")
cat("===========================================\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n\n")
cat("Package versions:\n")
for (pkg in required_packages) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\nComplete session information:\n")
print(sessionInfo())
sink()
cat(" ✅ 会话信息已保存\n")
# ============================================================================
# 15. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n12. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "12_paper2_pathway_analysis_L.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "12_code_list_L.txt")
cat("脚本名称: 12_paper2_pathway_analysis_L.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 16. 输出摘要
# ============================================================================
cat("\n========================================================\n")
cat("【Final Cluster Distribution】\n")
print(table(high_risk_data$pathway_cluster_en))
cat("\n【Cluster Means】\n")
print(cluster_profiles %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("\n========================================================\n")
cat("✅ Paper 2 analysis complete!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("clean_dir", "results_dir", "LOG_DIR")))
gc()
