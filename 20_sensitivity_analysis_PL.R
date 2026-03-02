#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 20_sensitivity_analysis.R
# 描述: 敏感性分析合集 - 排除服药者、非线性关系、青年亚组分析
# ============================================================================
# ============================================================================
# 1. 环境配置
# ============================================================================
rm(list = ls())
gc()
set.seed(20240226)
# 加载必要包
required_packages <- c(
  "tidyverse", "survey", "splines", "rms",
  "ggplot2", "gridExtra", "pROC"
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
hcf_labels <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
pathway_labels <- c(
  "低过激健康型" = "Low-hyperactivation",
  "对抗-枯竭混合型" = "Aversion-Exhaustion",
  "痴固着主导型" = "Perseveration",
  "纯生理过激型" = "Hyperactivation"
)
disease_labels <- c(
  "糖尿病" = "Diabetes",
  "高血压" = "Hypertension",
  "抑郁" = "Depression",
  "心血管疾病" = "CVD"
)
# ============================================================================
# 2. 配置路径
# ============================================================================
L_DATA_DIR <- "C:/NHANES_Data/CLEAN"
RESULTS_DIR <- "C:/NHANES_Data/CLEAN/results/sensitivity"
LOG_DIR <- "C:/NHANES_Data/CLEAN/logs"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志
log_file <- file.path(LOG_DIR, paste0("20_sensitivity_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 20_sensitivity_analysis.R\n")
cat("描述: 敏感性分析合集 - 排除服药者、非线性关系、青年亚组\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")
# 记录包版本
cat("加载的包版本:\n")
for (pkg in required_packages) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")
# ============================================================================
# 3. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
L_file <- file.path(L_DATA_DIR, "analysis_dataset_subset.rds")
if (!file.exists(L_file)) stop("错误: 找不到L周期数据")
data <- readRDS(L_file)
cat(sprintf("样本量: %d\n", nrow(data)))
# ============================================================================
# 3.1 合并通路聚类数据
# ============================================================================
cat("\n=== 合并通路聚类数据 ===\n")
pathway_file <- file.path(L_DATA_DIR, "results", "paper2", "paper2_analysis_data.rds")
if(file.exists(pathway_file)) {
  pathway_data <- readRDS(pathway_file)
  cat("通路数据加载成功，维度:", dim(pathway_data), "\n")
  cat("通路数据列名:", paste(names(pathway_data), collapse = ", "), "\n")
  # 合并前先重命名，避免冲突
  pathway_data <- pathway_data %>%
    rename(
      pathway_cluster_raw = pathway_cluster,
      pathway_cluster_en_raw = pathway_cluster_en
    )
  # 合并到主数据
  data <- data %>%
    left_join(pathway_data %>% select(SEQN, pathway_cluster_raw, pathway_cluster_en_raw), 
              by = "SEQN")
  # 重命名为标准名称
  data$pathway_cluster <- data$pathway_cluster_raw
  data$pathway_cluster_en <- data$pathway_cluster_en_raw
  cat("合并后数据维度:", dim(data), "\n")
  cat("合并后 pathway_cluster 缺失值:", sum(is.na(data$pathway_cluster)), "\n")
  cat("合并后 pathway_cluster_en 缺失值:", sum(is.na(data$pathway_cluster_en)), "\n")
} else {
  cat("⚠️ 通路数据文件不存在:", pathway_file, "\n")
}
# ============================================================================
# 3.2 创建设计对象（原有代码）
# ============================================================================
design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = data
)
cat(" ✅ 设计对象创建成功\n\n")
# ============================================================================
# 4. 敏感性分析1：排除服药人群 (Table S8)
# ============================================================================
cat("\n========================================================\n")
cat("2. 敏感性分析1：排除服药人群 (Table S8)\n")
cat("========================================================\n")
# 定义服药变量
data$on_medication <- with(data,
  (antihypertensive_med == 1 & !is.na(antihypertensive_med)) |
  (antidiabetic_med == 1 & !is.na(antidiabetic_med)) |
  (phq9_total >= 10 & !is.na(phq9_total))
)
cat(sprintf("服药人数: %d (%.1f%%)\n",
            sum(data$on_medication, na.rm = TRUE),
            mean(data$on_medication, na.rm = TRUE) * 100))
# 创建未服药子集
data_no_med <- data %>% filter(!on_medication)
cat(sprintf("未服药样本量: %d\n", nrow(data_no_med)))
# 创建设计对象
design_no_med <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = data_no_med
)
# 重新运行主要分析
diseases <- c("diabetes_doctor", "hypertension", "depression", "any_cvd")
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
# 疾病英文名映射
disease_names <- c(
  "diabetes_doctor" = "Diabetes",
  "hypertension" = "Hypertension", 
  "depression" = "Depression",
  "any_cvd" = "CVD"
)
sensitivity_results <- data.frame()
for(disease in diseases) {
  if(!disease %in% names(data_no_med)) next
  cat(sprintf("\n分析疾病: %s\n", disease))
  formula <- as.formula(paste0(disease, " ~ ", 
                                paste(alpha_vars, collapse = " + "),
                                " + RIDAGEYR + RIAGENDR"))
  model <- tryCatch({
    svyglm(formula, design = design_no_med, family = quasibinomial())
  }, error = function(e) {
    cat("  模型失败:", e$message, "\n")
    return(NULL)
  })
  if(!is.null(model)) {
    coefs <- summary(model)$coefficients
    for(alpha in alpha_vars) {
      if(alpha %in% rownames(coefs)) {
        sensitivity_results <- rbind(sensitivity_results, data.frame(
          Disease = disease_names[disease],
          Alpha_Factor = alpha,
          OR_No_Med = exp(coefs[alpha, "Estimate"]),
          CI_lower = exp(coefs[alpha, "Estimate"] - 1.96 * coefs[alpha, "Std. Error"]),
          CI_upper = exp(coefs[alpha, "Estimate"] + 1.96 * coefs[alpha, "Std. Error"]),
          P_value = coefs[alpha, "Pr(>|t|)"]
        ))
        cat(sprintf("  ✅ %s: OR=%.2f\n", alpha, exp(coefs[alpha, "Estimate"])))
      }
    }
  }
}
# 读取原始结果进行比较
original_file <- file.path("C:/NHANES_Data/CLEAN/results/paper3", "alpha_disease_all.csv")
if(file.exists(original_file) && nrow(sensitivity_results) > 0) {
  original <- read.csv(original_file)
  # 创建中文到英文的疾病映射
  disease_map <- c(
    "糖尿病" = "Diabetes",
    "高血压" = "Hypertension",
    "抑郁" = "Depression",
    "心血管疾病" = "CVD"
  )
  # 创建α因子映射
  alpha_map <- c(
    "alpha1" = "alpha1",
    "alpha2" = "alpha2", 
    "alpha3" = "alpha3",
    "alpha4" = "alpha4"
  )
original_subset <- original %>%
  filter(疾病 %in% names(disease_map)) %>%
  mutate(
    Disease = disease_map[as.character(疾病)],  
    Alpha_Factor = as.character(α因子)          
  ) %>%
  select(Disease, Alpha_Factor, OR) %>%
  rename(OR_Original = OR)
  sensitivity_results <- sensitivity_results %>%
    left_join(original_subset, by = c("Disease", "Alpha_Factor")) %>%
    mutate(
      Percent_Change = (OR_No_Med - OR_Original) / OR_Original * 100
    )
  cat("\n排除服药者敏感性分析结果:\n")
  print(sensitivity_results)
} else {
  if(nrow(sensitivity_results) == 0) {
    cat("\n⚠️ 没有生成敏感性分析结果\n")
  } else {
    cat("\n⚠️ 原始结果文件不存在，跳过比较\n")
  }
}
write.csv(sensitivity_results, file.path(RESULTS_DIR, "eTable12.csv"), row.names = FALSE)
cat("✅ Table S8 saved: eTable12.csv\n")
# ============================================================================
# 5. 敏感性分析2：α因子与炎症的非线性关系 (Figure S4 & Table S9)
# ============================================================================
cat("\n========================================================\n")
cat("3. 敏感性分析2：α因子与炎症的非线性关系 (Figure S4 & Table S9)\n")
cat("========================================================\n")
nonlinear_data <- data %>% filter(complete.cases(alpha2, hs_crp_mgl))
cat(sprintf("分析样本量: %d\n", nrow(nonlinear_data)))
cat("\nalpha2分布:\n")
print(summary(nonlinear_data$alpha2))
library(splines)
tryCatch({
  model_rcs <- lm(log(hs_crp_mgl + 0.1) ~ ns(alpha2, df = 3) + RIDAGEYR + RIAGENDR,
                  data = nonlinear_data)
  alpha_seq <- seq(min(nonlinear_data$alpha2, na.rm = TRUE),
                   max(nonlinear_data$alpha2, na.rm = TRUE), length.out = 100)
  pred <- predict(model_rcs, newdata = data.frame(
    alpha2 = alpha_seq,
    RIDAGEYR = mean(nonlinear_data$RIDAGEYR),
    RIAGENDR = 1
  ), se.fit = TRUE)
  nonlinear_results <- data.frame(
    alpha2 = alpha_seq,
    log_crp_pred = pred$fit,
    log_crp_lower = pred$fit - 1.96 * pred$se.fit,
    log_crp_upper = pred$fit + 1.96 * pred$se.fit,
    crp_pred = exp(pred$fit) - 0.1
  )
  model_linear <- lm(log(hs_crp_mgl + 0.1) ~ alpha2 + RIDAGEYR + RIAGENDR, data = nonlinear_data)
  anova_result <- anova(model_linear, model_rcs)
  cat(sprintf("\n非线性检验 p值: %.4f\n", anova_result$`Pr(>F)`[2]))
  # 可视化
  p <- ggplot(nonlinear_results, aes(x = alpha2, y = crp_pred)) +
    geom_line(color = "blue", size = 1.2) +
    geom_ribbon(aes(ymin = exp(log_crp_lower) - 0.1,
                    ymax = exp(log_crp_upper) - 0.1), alpha = 0.2) +
    labs(title = "Figure S4. Nonlinear Relationship Between α₂ and CRP",
         x = "α₂ (Emotional Regulation)", y = "Predicted CRP (mg/L)") +
    theme_minimal() +
    geom_vline(xintercept = c(-0.5, 0, 0.5), linetype = "dashed", color = "gray50") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(RESULTS_DIR, "eFigure4.pdf"), p, width = 8, height = 6)
  ggsave(file.path(RESULTS_DIR, "eFigure4.png"), p, width = 8, height = 6, dpi = 300)
  cat("✅ Figure S4 saved: eFigure4.pdf/.png\n")
  write.csv(nonlinear_results, file.path(RESULTS_DIR, "eTable13.csv"), row.names = FALSE)
  cat("✅ Table S9 saved: eTable13.csv\n")
}, error = function(e) {
  cat("⚠️ 非线性分析失败:", e$message, "\n")
})
# ============================================================================
# 6. 敏感性分析3：青年亚组分析 (18-29岁) - eTables 14-16
# ============================================================================
cat("\n========================================================\n")
cat("4. 敏感性分析3：青年亚组分析 (18-29岁)\n")
cat("========================================================\n")
# 筛选青年样本（18-29岁）
youth_data <- data %>% filter(RIDAGEYR >= 18 & RIDAGEYR <= 29)
cat(sprintf("青年样本量: %d\n", nrow(youth_data)))
# 创建青年亚组设计对象（survey包要求）
design_youth <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = youth_data
)
# --------------------------------------------------------------------------
# eTable 14: 青年亚组HCF分型分布
# --------------------------------------------------------------------------
cat("\n--- 生成 eTable 14: 青年亚组HCF分型分布 ---\n")
youth_hcf <- youth_data %>%
  mutate(HCF_type_en = hcf_labels[HCF_type]) %>%
  group_by(HCF_type_en) %>%
  summarise(
    N = n(),
    Pct = n() / nrow(youth_data) * 100,
    PHQ9_Mean = mean(phq9_total, na.rm = TRUE),
    Suicide_Rate = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(PHQ9_Mean))
cat("\n青年亚组HCF分型分布 (eTable 14):\n")
print(youth_hcf)
write.csv(youth_hcf, file.path(RESULTS_DIR, "eTable14.csv"), row.names = FALSE)
cat("✅ eTable14.csv 已保存\n")
# --------------------------------------------------------------------------
# eTable 15: 青年亚组α因子效应（加权）
# --------------------------------------------------------------------------
cat("\n--- 生成 eTable 15: 青年亚组α因子效应 ---\n")
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
youth_alpha_effects <- data.frame()
for(alpha in alpha_vars) {
  # 确保变量存在
  if(!alpha %in% names(youth_data)) {
    cat(sprintf("⚠️ 变量 %s 不存在，跳过\n", alpha))
    next
  }
  # 构建公式
  formula <- as.formula(paste0("phq9_total ~ ", alpha, " + RIDAGEYR + RIAGENDR"))
  # 拟合加权回归模型
  model <- tryCatch({
    svyglm(formula, design = design_youth)
  }, error = function(e) {
    cat(sprintf("⚠️ 模型拟合失败 (%s): %s\n", alpha, e$message))
    return(NULL)
  })
  if(!is.null(model)) {
    coefs <- summary(model)$coefficients
    if(alpha %in% rownames(coefs)) {
      youth_alpha_effects <- rbind(youth_alpha_effects, data.frame(
        Alpha_Factor = alpha,
        Beta = coefs[alpha, "Estimate"],
        SE = coefs[alpha, "Std. Error"],
        P_value = coefs[alpha, "Pr(>|t|)"],
        CI_lower = coefs[alpha, "Estimate"] - 1.96 * coefs[alpha, "Std. Error"],
        CI_upper = coefs[alpha, "Estimate"] + 1.96 * coefs[alpha, "Std. Error"],
        stringsAsFactors = FALSE
      ))
      cat(sprintf("  ✅ %s: Beta = %.3f, P = %.4f\n", 
                  alpha, coefs[alpha, "Estimate"], coefs[alpha, "Pr(>|t|)"]))
    }
  }
}
cat("\n青年亚组α因子效应 (eTable 15):\n")
print(youth_alpha_effects)
write.csv(youth_alpha_effects, file.path(RESULTS_DIR, "eTable15.csv"), row.names = FALSE)
cat("✅ eTable15.csv 已保存\n")
# --------------------------------------------------------------------------
# eTable 16: 青年中最脆弱的群体（HCF分型 × 四通路）
# --------------------------------------------------------------------------
cat("\n--- 生成 eTable 16: 青年中最脆弱的群体 ---\n")
# 检查通路变量
cat("通路变量检查:\n")
cat("  pathway_cluster 是否存在:", "pathway_cluster" %in% names(youth_data), "\n")
if("pathway_cluster" %in% names(youth_data)) {
  cat("  pathway_cluster 非缺失值:", sum(!is.na(youth_data$pathway_cluster)), "\n")
  cat("  pathway_cluster 分布:\n")
  print(table(youth_data$pathway_cluster, useNA = "ifany"))
}
# 只保留有通路分类的样本
youth_with_pathway <- youth_data %>% 
  filter(!is.na(pathway_cluster) & !is.na(HCF_type))
cat(sprintf("\n有完整HCF和通路分类的青年样本: %d (%.1f%%)\n", 
            nrow(youth_with_pathway), 
            nrow(youth_with_pathway)/nrow(youth_data)*100))
if(nrow(youth_with_pathway) >= 10) {
  youth_severe <- youth_with_pathway %>%
    mutate(
      HCF_type_en = hcf_labels[HCF_type],
      Pathway_en = pathway_labels[as.character(pathway_cluster)]
    ) %>%
    group_by(HCF_type_en, Pathway_en) %>%
    summarise(
      N = n(),
      PHQ9_Mean = mean(phq9_total, na.rm = TRUE),
      Depression_Rate = mean(depression, na.rm = TRUE) * 100,
      Suicide_Rate = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
      Alpha2_Mean = mean(alpha2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(PHQ9_Mean))
  cat("\n青年中最脆弱的群体 (eTable 16):\n")
  print(youth_severe)
  write.csv(youth_severe, file.path(RESULTS_DIR, "eTable16.csv"), row.names = FALSE)
  cat("✅ eTable16.csv 已保存\n")
} else {
  cat("\n⚠️ 有通路分类的样本太少，无法生成可靠的 eTable16\n")
  # 生成简化版本（仅HCF分型）
  youth_severe_simple <- youth_data %>%
    mutate(HCF_type_en = hcf_labels[HCF_type]) %>%
    group_by(HCF_type_en) %>%
    summarise(
      N = n(),
      PHQ9_Mean = mean(phq9_total, na.rm = TRUE),
      Depression_Rate = mean(depression, na.rm = TRUE) * 100,
      Suicide_Rate = mean(DPQ090 >= 1, na.rm = TRUE) * 100,
      Alpha2_Mean = mean(alpha2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(PHQ9_Mean))
  cat("\n简化版本 (仅HCF分型):\n")
  print(youth_severe_simple)
  write.csv(youth_severe_simple, file.path(RESULTS_DIR, "eTable16_simple.csv"), row.names = FALSE)
  cat("✅ eTable16_simple.csv 已保存\n")
}
# --------------------------------------------------------------------------
# 总结
# --------------------------------------------------------------------------
cat("\n✅ 青年亚组分析完成\n")
cat("   生成文件:\n")
cat("   - eTable14.csv: 青年亚组HCF分型分布\n")
cat("   - eTable15.csv: 青年亚组α因子效应\n")
if(exists("youth_severe")) {
  cat("   - eTable16.csv: 青年中最脆弱的群体（HCF×通路）\n")
} else {
  cat("   - eTable16_simple.csv: 青年中最脆弱的群体（简化版）\n")
}
# ============================================================================
# 7. 生成敏感性分析报告
# ============================================================================
cat("\n========================================================\n")
cat("5. 生成敏感性分析报告\n")
cat("========================================================\n")
sink(file.path(LOG_DIR, "20_sensitivity_report.txt"))
cat("Sensitivity Analysis Report\n")
cat("===========================\n\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、Excluding Medicated Participants (Table S8)\n")
if(exists("sensitivity_results")) print(head(sensitivity_results))
cat("\n二、Nonlinear Relationship Between α₂ and CRP (Figure S4, Table S9)\n")
if(exists("anova_result")) {
  cat(sprintf("Nonlinearity test p-value: %.4f\n", anova_result$`Pr(>F)`[2]))
}
cat("\n三、Youth Subgroup Analysis (Tables S10-12)\n")
cat("\nHCF Type Distribution in Youth (Table S10):\n")
print(youth_hcf)
cat("\nAlpha Factor Effects in Youth (Table S11):\n")
print(youth_alpha_effects)
cat("\n四、Output Files\n")
cat(" - eTable12.csv\n")
cat(" - eFigure4.pdf\n")
cat(" - eFigure4.png\n")
cat(" - eTable13.csv\n")
cat(" - eTable14.csv\n")
cat(" - eTable15.csv\n")
cat(" - eTable16.csv\n")
sink()
cat("\n✅ Sensitivity analysis report saved\n")
# ============================================================================
# 8. 保存会话信息
# ============================================================================
cat("\n6. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "20_session_info.txt")
sink(session_info_path)
cat("Sensitivity Analysis Session Information\n")
cat("=======================================\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n\n")
cat("Package versions:\n")
for (pkg in required_packages) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\nComplete session information:\n")
print(sessionInfo())
sink()
cat(" ✅ Session information saved\n")
# ============================================================================
# 9. 保存R代码副本
# ============================================================================
cat("\n7. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) dir.create(scripts_dir, recursive = TRUE)
code_save_path <- file.path(scripts_dir, "20_sensitivity_analysis.R")
cat("\n⚠️  Please manually save the current script to:\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   This is a JAMA Psychiatry requirement: all analysis code must be saved and made public.\n")
code_list_path <- file.path(LOG_DIR, "20_code_list.txt")
cat("Script name: 20_sensitivity_analysis.R\n", file = code_list_path)
cat("Generation time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("Suggested save location:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ Code list saved\n")
# ============================================================================
# 10. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ Sensitivity analysis complete!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results saved to:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
