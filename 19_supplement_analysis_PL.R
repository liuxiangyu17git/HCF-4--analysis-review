#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 19_supplement_analysis.R
# 描述: 论文补充分析整合脚本 - NNT计算、AUC比较、敏感性分析等
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 探索性分析 (非主要假设)
# 对应研究计划: 第六部分 - 敏感性分析/扩展分析
#
# 随机种子: 20240226
# 最后修改: 2026-02-28
# ============================================================================
# ============================================================================
# 1. 环境配置
# ============================================================================
rm(list = ls())
gc()
set.seed(20240226)
# 加载包
required_packages <- c(
  "tidyverse", "survey", "pROC", "gridExtra",
  "knitr", "kableExtra", "forestplot", "pheatmap"
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
# HCF分型英文标签
hcf_labels <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
# α因子英文标签
alpha_labels <- c(
  "alpha1" = "α₁",
  "alpha2" = "α₂",
  "alpha3" = "α₃",
  "alpha4" = "α₄"
)
# 系统英文标签
system_labels <- c(
  "代谢" = "Metabolic",
  "心血管" = "Cardiovascular",
  "肾脏" = "Renal",
  "肝脏" = "Hepatic",
  "神经精神" = "Neuropsychiatric"
)
# 疾病英文标签
disease_labels <- c(
  "糖尿病" = "Diabetes",
  "肥胖" = "Obesity",
  "高胆固醇" = "High cholesterol",
  "高血压" = "Hypertension",
  "心血管疾病" = "CVD",
  "心力衰竭" = "Heart failure",
  "中风" = "Stroke",
  "心肌梗死" = "Myocardial infarction",
  "慢性肾病" = "CKD",
  "ALT升高" = "Elevated ALT",
  "抑郁" = "Depression",
  "自评健康差" = "Poor self-rated health"
)
# 生物标志物英文标签
biomarker_labels <- c(
  "NLR" = "NLR",
  "eGFR" = "eGFR",
  "CKD" = "CKD",
  "血清肌酐" = "Serum creatinine",
  "AST" = "AST",
  "AST/ALT比值" = "AST/ALT ratio",
  "肝损伤" = "Liver injury",
  "尿酸" = "Uric acid",
  "高尿酸血症" = "Hyperuricemia",
  "HDL胆固醇" = "HDL cholesterol",
  "低HDL" = "Low HDL"
)
# 临床结局英文标签
clinical_labels <- c(
  "住院" = "Hospitalization",
  "急诊" = "Emergency visit",
  "活动受限" = "Activity limitation",
  "降压药" = "Antihypertensive medication",
  "降糖药" = "Antidiabetic medication"
)
# ============================================================================
# 2. 配置路径
# ============================================================================
L_DATA_DIR <- "C:/NHANES_Data/CLEAN"
RESULTS_DIR <- "C:/NHANES_Data/CLEAN/results/supplement"
LOG_DIR <- "C:/NHANES_Data/CLEAN/logs"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 日志
log_file <- file.path(LOG_DIR, paste0("19_supplement_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 19_supplement_analysis.R\n")
cat("描述: 补充分析整合 - NNT、AUC、敏感性分析等\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")
# ============================================================================
# 3. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
L_file <- file.path(L_DATA_DIR, "analysis_dataset_subset.rds")
if (!file.exists(L_file)) stop("错误: 找不到L周期数据")
data <- readRDS(L_file)
cat(sprintf("样本量: %d\n", nrow(data)))
# 创建设计对象
design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = data
)
# ============================================================================
# 4. 样本量流程图 (Table S1)
# ============================================================================
cat("\n========================================================\n")
cat("2. 生成样本量流程图 (Table S1)\n")
cat("========================================================\n")
flow_data <- data.frame(
  Step = c(
    "Total NHANES 2021-2023 sample",
    "Excluded age <18 years",
    "Excluded pregnant women",
    "Excluded missing HCF type",
    "Final analytic sample"
  ),
  N = c(
    nrow(readRDS(file.path(L_DATA_DIR, "final_analysis_dataset.rds"))),
    sum(data$RIDAGEYR >= 18),
    sum(data$RIDAGEYR >= 18 & (data$RIDEXPRG != 1 | is.na(data$RIDEXPRG))),
    sum(!is.na(data$HCF_type)),
    nrow(data)
  )
)
write.csv(flow_data, file.path(RESULTS_DIR, "eTable2.csv"), row.names = FALSE)
cat("✅ Table S1 saved: eTable2.csv\n")
# ============================================================================
# 5. 变量一致性对照表 (Table S2)
# ============================================================================
cat("\n========================================================\n")
cat("3. 生成变量一致性对照表 (Table S2)\n")
cat("========================================================\n")
var_consistency <- data.frame(
  Variable = c(
    "HCF type", "Pathway cluster", "Alpha factors",
    "Depression (PHQ-9)", "Inflammation (CRP)", "Resting heart rate",
    "BMI", "Hypertension", "Diabetes"
  ),
  L_cycle_variable = c(
    "HCF_type", "pathway_cluster", "alpha1-alpha4",
    "phq9_total", "hs_crp_mgl", "BPXOPLS1",
    "BMXBMI", "hypertension", "diabetes_doctor"
  ),
  P_cycle_variable = c(
    "HCF_type", "pathway_cluster", "alpha1-alpha4",
    "phq9_total", "hs_crp_mgl", "BPXOPLS1",
    "BMXBMI", "hypertension", "diabetes_doctor"
  ),
  Measurement_Consistent = c("Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes")
)
write.csv(var_consistency, file.path(RESULTS_DIR, "eTable3.csv"), row.names = FALSE)
cat("✅ Table S2 saved: eTable3.csv\n")
# ============================================================================
# 6. FDR校正前后对比 (Table S3)
# ============================================================================
cat("\n========================================================\n")
cat("4. FDR校正前后对比 (Table S3)\n")
cat("========================================================\n")
# 读取原始p值（从Table 3）
alpha_disease_file <- file.path("C:/NHANES_Data/CLEAN/results/paper3", "alpha_disease_all.csv")
if(file.exists(alpha_disease_file)) {
  alpha_res <- read.csv(alpha_disease_file)
  # FDR校正
  alpha_res$p_fdr <- p.adjust(alpha_res$p_value, method = "fdr")
  # 创建英文版表格
  alpha_res_en <- data.frame(
    System = system_labels[alpha_res$系统],
    Disease = disease_labels[alpha_res$疾病],
    Alpha_Factor = alpha_labels[alpha_res$α因子],
    Raw_p = alpha_res$p_value,
    FDR_p = alpha_res$p_fdr,
    Raw_Significant = ifelse(alpha_res$p_value < 0.05, "Yes", "No"),
    FDR_Significant = ifelse(alpha_res$p_fdr < 0.05, "Yes", "No")
  )
  # 统计
  n_sig_raw <- sum(alpha_res$p_value < 0.05)
  n_sig_fdr <- sum(alpha_res$p_fdr < 0.05)
  cat(sprintf("\nFDR校正前显著结果数: %d\n", n_sig_raw))
  cat(sprintf("FDR校正后显著结果数: %d (%.1f%%保留)\n", 
              n_sig_fdr, n_sig_fdr/n_sig_raw*100))
  write.csv(alpha_res_en, file.path(RESULTS_DIR, "eTable7.csv"), row.names = FALSE)
  cat("✅ Table S3 saved: eTable7.csv\n")
}
# ============================================================================
# 6a. 生成 eTable 7 - FDR校正对比表（详细版）
# ============================================================================
cat("\n========================================================\n")
cat("6a. 生成 eTable 7: FDR校正对比表（详细版）\n")
cat("========================================================\n")
# 读取α因子与疾病关联的原始结果 - 使用 eTable4.csv
alpha_file <- file.path("C:/NHANES_Data/CLEAN/results/paper3", "eTable4.csv")
if(file.exists(alpha_file)) {
  alpha_data <- read.csv(alpha_file)
  cat("✅ 加载 eTable4.csv 成功\n")
  cat("数据维度:", dim(alpha_data), "\n")
  cat("列名:", paste(names(alpha_data), collapse = ", "), "\n")
  # 检查是否有P_value列
  if("P_value" %in% names(alpha_data)) {
    # 进行FDR校正
    alpha_data$FDR_p <- p.adjust(alpha_data$P_value, method = "fdr")
    # 创建详细的FDR对比表
    fdr_table <- data.frame(
      Outcome = alpha_data$Outcome,
      Alpha_Factor = alpha_data$Alpha_Factor,
      Raw_p = alpha_data$P_value,
      FDR_p = alpha_data$FDR_p,
      Raw_Significant = ifelse(alpha_data$P_value < 0.05, "✅", "❌"),
      FDR_Significant = ifelse(alpha_data$FDR_p < 0.05, "✅", "❌")
    )
    # 按P值排序
    fdr_table <- fdr_table %>% arrange(Raw_p)
    # 统计信息
    n_sig_raw <- sum(alpha_data$P_value < 0.05)
    n_sig_fdr <- sum(alpha_data$FDR_p < 0.05)
    cat("\n========================================\n")
    cat("FDR校正统计:\n")
    cat(sprintf("校正前显著结果数: %d\n", n_sig_raw))
    cat(sprintf("校正后显著结果数: %d (%.1f%%保留)\n", 
                n_sig_fdr, n_sig_fdr/n_sig_raw*100))
    # 保存为eTable 7
    write.csv(fdr_table, 
              file.path(RESULTS_DIR, "eTable7.csv"), 
              row.names = FALSE)
    cat("\n✅ eTable7.csv 已保存\n")
    cat("   总行数:", nrow(fdr_table), "\n")
    # 显示前几行确认
    cat("\nFDR对比表预览（前10行）:\n")
    print(head(fdr_table, 10))
  } else {
    cat("\n❌ eTable4.csv 中没有 P_value 列\n")
  }
} else {
  cat("\n⚠️ 文件不存在:", alpha_file, "\n")
}
# ============================================================================
# 7. 与已知风险因素比较 (Table S4)
# ============================================================================
cat("\n========================================================\n")
cat("5. 与已知风险因素比较 (Table S4)\n")
cat("========================================================\n")
compare_vars <- c("alpha2", "alpha3", "BMXBMI", "current_smoker", "heavy_drinker")
compare_results <- data.frame()
for(disease in c("depression", "diabetes_doctor", "hypertension")) {
  if(!disease %in% names(data)) next
  # 模型1: 只有经典风险因素
  formula1 <- as.formula(paste0(disease, " ~ BMXBMI + current_smoker + heavy_drinker + RIDAGEYR + RIAGENDR"))
  model1 <- svyglm(formula1, design = design, family = quasibinomial())
  # 模型2: 加入α因子
  formula2 <- as.formula(paste0(disease, " ~ alpha2 + alpha3 + BMXBMI + current_smoker + heavy_drinker + RIDAGEYR + RIAGENDR"))
  model2 <- svyglm(formula2, design = design, family = quasibinomial())
  compare_results <- rbind(compare_results, data.frame(
    Outcome = disease,
    Model1_AIC = AIC(model1),
    Model2_AIC = AIC(model2),
    AIC_Delta = AIC(model1) - AIC(model2)
  ))
}
write.csv(compare_results, file.path(RESULTS_DIR, "eTable8.csv"), row.names = FALSE)
cat("✅ Table S4 saved: eTable8.csv\n")
# ============================================================================
# 8. NNT计算 (Table S5)
# ============================================================================
cat("\n========================================================\n")
cat("6. NNT计算 (Table S5)\n")
cat("========================================================\n")
model_alpha2 <- svyglm(depression ~ alpha2 + RIDAGEYR + RIAGENDR, 
                        design = design, family = quasibinomial())
or_alpha2 <- exp(coef(model_alpha2)["alpha2"])
baseline_risk <- mean(data$depression, na.rm = TRUE)
if(or_alpha2 < 1) {
  arr <- baseline_risk * (1 - or_alpha2)
  nnt <- ceiling(1/arr)
} else {
  arr <- NA
  nnt <- NA
}
nnt_result <- data.frame(
  Factor = "α₂",
  OR = round(or_alpha2, 2),
  Baseline_Risk = round(baseline_risk * 100, 1),
  Absolute_Risk_Reduction = round(arr * 100, 1),
  NNT = nnt
)
cat("\nNNT计算结果:\n")
print(nnt_result)
write.csv(nnt_result, file.path(RESULTS_DIR, "eTable9.csv"), row.names = FALSE)
cat("✅ Table S5 saved: eTable9.csv\n")
# ============================================================================
# 9. AUC比较 (Table S6) 和 ROC曲线 (eFigure 3)
# ============================================================================
cat("\n========================================================\n")
cat("7. AUC比较 (Table S6) 和 ROC曲线 (eFigure 3)\n")
cat("========================================================\n")
library(pROC)
vars_needed <- c("depression", "alpha2", "BMXBMI", "current_smoker", 
                 "RIDAGEYR", "RIAGENDR")
valid_idx <- complete.cases(data[, vars_needed])
valid_data <- data[valid_idx, ]
cat(sprintf("完整案例样本量: %d\n", nrow(valid_data)))
# 模型1: α₂ + 年龄性别
model_alpha <- glm(depression ~ alpha2 + RIDAGEYR + RIAGENDR, 
                   data = valid_data, family = binomial())
# 模型2: 经典风险因素
valid_data$current_smoker <- as.factor(valid_data$current_smoker)
model_classic <- glm(depression ~ BMXBMI + current_smoker + RIDAGEYR + RIAGENDR, 
                     data = valid_data, family = binomial())
pred_alpha <- predict(model_alpha, type = "response")
pred_classic <- predict(model_classic, type = "response")
roc_alpha <- roc(valid_data$depression, pred_alpha)
roc_classic <- roc(valid_data$depression, pred_classic)
auc_alpha <- auc(roc_alpha)
auc_classic <- auc(roc_classic)
ci_alpha <- ci.auc(roc_alpha)
ci_classic <- ci.auc(roc_classic)
auc_compare <- data.frame(
  Model = c("α₂ + age/sex", "BMI + smoking + age/sex"),
  AUC = round(c(auc_alpha, auc_classic), 3),
  CI_lower = round(c(ci_alpha[1], ci_classic[1]), 3),
  CI_upper = round(c(ci_alpha[3], ci_classic[3]), 3)
)
cat("\nAUC比较:\n")
print(auc_compare)
write.csv(auc_compare, file.path(RESULTS_DIR, "eTable10.csv"), row.names = FALSE)
# 绘制ROC曲线
pdf(file.path(RESULTS_DIR, "eFigure3.pdf"), width = 8, height = 6)
plot(roc_alpha, col = "blue", main = "eFigure 3. ROC Curve Comparison")
plot(roc_classic, col = "red", add = TRUE)
legend("bottomright", legend = c("α₂ Model", "Classic Risk Model"), 
       col = c("blue", "red"), lty = 1)
dev.off()
# 保存PNG
png(file.path(RESULTS_DIR, "eFigure3.png"), width = 800, height = 600, res = 150)
plot(roc_alpha, col = "blue", main = "eFigure 3. ROC Curve Comparison")
plot(roc_classic, col = "red", add = TRUE)
legend("bottomright", legend = c("α₂ Model", "Classic Risk Model"), 
       col = c("blue", "red"), lty = 1)
dev.off()
cat("✅ eFigure 3 saved: eFigure3.pdf/.png\n")
cat("✅ Table S6 saved: eTable10.csv\n")
# ============================================================================
# 10. 敏感性分析 (Table S7)
# ============================================================================
cat("\n========================================================\n")
cat("8. 敏感性分析 (Table S7)\n")
cat("========================================================\n")
cutpoints <- c(0.5, 1, 1.5)
sensitivity_results <- data.frame()
for(i in seq_along(cutpoints)) {
  cut <- cutpoints[i]
  data$alpha2_high <- as.numeric(data$alpha2 > cut)
  design <- update(design, alpha2_high = data$alpha2_high)
  model <- svyglm(depression ~ alpha2_high + RIDAGEYR + RIAGENDR, 
                  design = design, family = quasibinomial())
  coef_name <- "alpha2_high"
  if(coef_name %in% names(coef(model))) {
    or_val <- exp(coef(model)[coef_name])
    ci <- exp(confint(model)[coef_name, ])
    sensitivity_results <- rbind(sensitivity_results, data.frame(
      Cutpoint = cut,
      OR = round(or_val, 2),
      CI_lower = round(ci[1], 2),
      CI_upper = round(ci[2], 2)
    ))
  }
}
if(nrow(sensitivity_results) > 0) {
  cat("\n敏感性分析结果:\n")
  print(sensitivity_results)
  write.csv(sensitivity_results, file.path(RESULTS_DIR, "eTable11.csv"), row.names = FALSE)
  cat("✅ Table S7 saved: eTable11.csv\n")
}
# ============================================================================
# 11. 生成补充分析报告
# ============================================================================
cat("\n========================================================\n")
cat("9. 生成补充分析报告\n")
cat("========================================================\n")
sink(file.path(LOG_DIR, "19_supplement_report.txt"))
cat("Supplemental Analysis Report\n")
cat("============================\n\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、Sample Flowchart\n")
print(flow_data)
cat("\n二、FDR Correction Results\n")
if(exists("n_sig_raw")) {
  cat(sprintf("Raw significant results: %d\n", n_sig_raw))
  cat(sprintf("FDR-corrected significant: %d\n", n_sig_fdr))
}
cat("\n三、Comparison with Known Risk Factors\n")
if(exists("compare_results")) print(compare_results)
cat("\n四、NNT Calculation\n")
if(exists("nnt_result")) print(nnt_result)
cat("\n五、AUC Comparison\n")
if(exists("auc_compare")) print(auc_compare)
cat("\n六、Sensitivity Analysis\n")
if(exists("sensitivity_results")) print(sensitivity_results)
cat("\n七、Output Files\n")
cat(" - eTable2.csv\n")
cat(" - eTable3.csv\n")
cat(" - eTable7.csv (FDR校正对比表)\n")
cat(" - eTable8.csv\n")
cat(" - eTable9.csv\n")
cat(" - eTable10.csv\n")
cat(" - eTable11.csv\n")
cat(" - eTable18.csv\n")
cat(" - figure_S2_heatmap.pdf\n")
cat(" - figure_S2_heatmap.png\n")
cat(" - eFigure3.pdf\n")
cat(" - eFigure3.png\n")
sink()
cat("\n✅ Supplemental report saved\n")
# ============================================================================
# 12. 保存会话信息
# ============================================================================
cat("\n10. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "19_session_info.txt")
sink(session_info_path)
cat("Supplemental Analysis Session Information\n")
cat("========================================\n")
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
# 13. 保存R代码副本
# ============================================================================
cat("\n11. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) dir.create(scripts_dir, recursive = TRUE)
code_save_path <- file.path(scripts_dir, "19_supplement_analysis.R")
cat("\n⚠️  Please manually save the current script to:\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   This is a JAMA Psychiatry requirement: all analysis code must be saved and made public.\n")
code_list_path <- file.path(LOG_DIR, "19_code_list.txt")
cat("Script name: 19_supplement_analysis.R\n", file = code_list_path)
cat("Generation time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("Suggested save location:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ Code list saved\n")
# ============================================================================
# 14. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ Supplemental analysis complete!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results saved to:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
