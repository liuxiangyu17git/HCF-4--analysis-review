#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 10_paper1_supplement_P.R
# 描述: NHANES 2017-2020 (P周期) 论文1补充分析
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 假设 H3 (C层中介), H5 (行为交互)
# 对应研究计划: 第六部分 - 敏感性分析
# 对应变量详表: 第四部分 HCF分型
#
# 伦理声明: 使用NHANES公开数据，无需IRB批准
# 数据来源: https://www.cdc.gov/nchs/nhanes/index.htm
# 数据周期: 2017-2020 (P系列)
#
# 随机种子: 20240226 (固定以确保可重复性)
# 最后修改: 2026-02-21
# ============================================================================
# ============================================================================
# 1. 环境配置
# ============================================================================
rm(list = ls())
gc()
set.seed(20240226)
# 加载必要包
required_packages <- c("survey", "dplyr", "emmeans", "ggplot2",
                       "tidyr", "pheatmap", "RColorBrewer")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# ============================================================================
# 英文标签定义
# ============================================================================
# HCF分型英文映射
hcf_labels <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
# 临床结局英文标签
clinical_labels <- c(
  "住院" = "Hospitalization",
  "急诊" = "Emergency Visit",
  "活动受限" = "Activity Limitation",
  "降压药" = "Antihypertensive Medication",
  "降糖药" = "Antidiabetic Medication"
)
# 生物标志物英文标签
biomarker_labels <- c(
  "NLR" = "NLR",
  "eGFR" = "eGFR",
  "CKD" = "CKD",
  "血清肌酐" = "Serum Creatinine",
  "AST" = "AST",
  "AST/ALT比值" = "AST/ALT Ratio",
  "肝损伤" = "Liver Injury",
  "尿酸" = "Uric Acid",
  "高尿酸血症" = "Hyperuricemia",
  "HDL胆固醇" = "HDL Cholesterol",
  "低HDL" = "Low HDL"
)
# ============================================================================
# 2. 配置路径 - P周期独立目录
# ============================================================================
PROJECT_ROOT <- "C:/NHANES_Data"
CLEAN_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")
RESULTS_DIR <- file.path(CLEAN_DATA_DIR, "results")
LOG_DIR <- file.path(CLEAN_DATA_DIR, "logs")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 日志
log_file <- file.path(LOG_DIR, paste0("10_paper1_supp_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
# ============================================================================
# 3. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
df <- readRDS(file.path(CLEAN_DATA_DIR, "final_analysis_dataset_P.rds"))
df <- df %>% filter(in_analysis == 1)
# 加载健康行为变量
health <- readRDS(file.path(CLEAN_DATA_DIR, "healthbehavior_vars_P.rds"))
# 使用 match() 直接赋值
df$pa_meets_guideline <- health$pa_meets_guideline[match(df$SEQN, health$SEQN)]
df$sleep_adequate <- health$sleep_adequate[match(df$SEQN, health$SEQN)]
cat(sprintf("分析样本: %d人\n", nrow(df)))
cat("\npa_meets_guideline 分布:\n")
print(table(df$pa_meets_guideline, useNA = "ifany"))
# ============================================================================
# 4. 创建设计对象
# ============================================================================
cat("\n2. 创建设计对象...\n")
mec_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMECPRP,
  nest = TRUE,
  data = df
)
cat(" ✅ 设计对象创建成功\n\n")
# ============================================================================
# 5. 分析1：C层中介模型
# ============================================================================
cat("========================================================\n")
cat("分析1：C层中介模型\n")
cat("========================================================\n")
model_total <- svyglm(phq9_total ~ HCF_A + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                      design = mec_design)
model_c <- svyglm(HCF_C ~ HCF_A + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                  family = quasibinomial(), design = mec_design)
model_outcome <- svyglm(phq9_total ~ HCF_A + HCF_C + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                        design = mec_design)
total_effect <- coef(model_total)["HCF_A"]
a_effect <- coef(model_c)["HCF_A"]
b_effect <- coef(model_outcome)["HCF_C"]
direct_effect <- coef(model_outcome)["HCF_A"]
indirect_effect <- a_effect * b_effect
prop_mediated <- indirect_effect / total_effect * 100
cat("\n=== C层中介分析结果 ===\n")
cat(sprintf("总效应: %.3f\n", total_effect))
cat(sprintf("直接效应: %.3f\n", direct_effect))
cat(sprintf("间接效应: %.3f\n", indirect_effect))
cat(sprintf("中介比例: %.1f%%\n", prop_mediated))
# 保存中介结果（英文列名）
mediation_results <- data.frame(
  Effect = c("Total", "Direct", "Indirect", "Prop_mediated"),
  Estimate = c(total_effect, direct_effect, indirect_effect, prop_mediated)
)
write.csv(mediation_results, file.path(RESULTS_DIR, "paper1_supp_table1_mediation_P.csv"), row.names = FALSE)
cat("✅ Table 1 (Mediation results) saved\n\n")
# ============================================================================
# 6. 分析2：分型特异性行为效应
# ============================================================================
cat("\n========================================================\n")
cat("分析2：分型特异性行为效应\n")
cat("========================================================\n")
# 只分析三种类型（排除健康型）
hcf_types <- c("纯心理型", "纯生理型", "身心混合型")
pa_results <- data.frame()
for(type in hcf_types) {
  cat(sprintf("\n分析 %s...\n", type))
  subset_design <- tryCatch({
    subset(mec_design, HCF_type == type)
  }, error = function(e) {
    cat("  子集创建失败:", e$message, "\n")
    return(NULL)
  })
  if(is.null(subset_design)) next
  # 检查该亚组中 pa_meets_guideline 的分布
  vals <- df$pa_meets_guideline[df$HCF_type == type]
  cat(sprintf("  样本量: %d\n", length(vals)))
  cat(sprintf("  0: %d, 1: %d, NA: %d\n", 
              sum(vals == 0, na.rm = TRUE),
              sum(vals == 1, na.rm = TRUE),
              sum(is.na(vals))))
  model <- tryCatch({
    svyglm(phq9_total ~ pa_meets_guideline + RIDAGEYR + RIAGENDR,
           design = subset_design)
  }, error = function(e) {
    cat("  模型失败:", e$message, "\n")
    return(NULL)
  })
  if(!is.null(model)) {
    coefs <- summary(model)$coefficients
    if("pa_meets_guideline" %in% rownames(coefs)) {
      pa_results <- rbind(pa_results, data.frame(
        HCF_type = type,
        Beta = coefs["pa_meets_guideline", "Estimate"],
        SE = coefs["pa_meets_guideline", "Std. Error"],
        P_value = coefs["pa_meets_guideline", "Pr(>|t|)"],
        N = sum(df$HCF_type == type, na.rm = TRUE)
      ))
      cat(sprintf("  ✅ beta=%.3f, p=%.4f\n", 
                  coefs["pa_meets_guideline", "Estimate"],
                  coefs["pa_meets_guideline", "Pr(>|t|)"]))
    }
  }
}
# 保存行为效应结果（英文列名）
if(nrow(pa_results) > 0) {
  # 将HCF_type列的值映射为英文
  pa_results$HCF_type_en <- hcf_labels[pa_results$HCF_type]
  pa_results_out <- data.frame(
    HCF_Type = pa_results$HCF_type_en,
    Beta = pa_results$Beta,
    SE = pa_results$SE,
    P_value = pa_results$P_value,
    N = pa_results$N
  )
  write.csv(pa_results_out, file.path(RESULTS_DIR, "paper1_supp_table2_pa_effects_P.csv"), row.names = FALSE)
  cat("\n✅ Table 2 (Behavior effects) saved\n")
  # 绘图
  plot_data <- data.frame(
    HCF_type = pa_results$HCF_type_en,
    Beta = pa_results$Beta,
    CI_lower = pa_results$Beta - 1.96 * pa_results$SE,
    CI_upper = pa_results$Beta + 1.96 * pa_results$SE
  )
  # 颜色映射
  color_map <- c(
    "Psychological" = "#FF7F00",
    "Pure Physiological" = "#377EB8",
    "Psychosomatic Mixed" = "#E41A1C"
  )
  p <- ggplot(plot_data, aes(x = HCF_type, y = Beta, fill = HCF_type)) +
    geom_col() +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = color_map) +
    labs(title = "Figure S1. Physical Activity Effects on Depression by HCF Type (P Cycle)",
         x = "", y = "β Coefficient (95% CI)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  # 同时保存PDF和PNG
  ggsave(file.path(RESULTS_DIR, "paper1_figure_s1_pa_effects_P.pdf"), p, width = 8, height = 6)
  ggsave(file.path(RESULTS_DIR, "paper1_figure_s1_pa_effects_P.png"), p, width = 8, height = 6, dpi = 300)
  cat(" ✅ Figure S1 saved: PDF and PNG formats (P Cycle)\n")
}
# ============================================================================
# 7. 分析3：生物标志物梯度分析
# ============================================================================
cat("\n========================================================\n")
cat("分析3：生物标志物梯度分析\n")
cat("========================================================\n")
biomarkers <- c(
  "nlr_ratio", "egfr_ckdepi_2021", "ckd_flag", "serum_creatinine_mgdl",
  "ast_ui_l", "ast_alt_ratio", "liver_injury", "uric_acid_mgdl",
  "hyperuricemia", "hdl_cholesterol_mgdl", "low_hdl"
)
names(biomarkers) <- c("NLR", "eGFR", "CKD", "Serum creatinine",
                        "AST", "AST/ALT ratio", "Liver injury", "Uric acid",
                        "Hyperuricemia", "HDL cholesterol", "Low HDL")
biomarker_results <- data.frame()
for(i in 1:length(biomarkers)) {
  bio_name <- names(biomarkers)[i]
  bio_var <- biomarkers[i]
  if(!bio_var %in% names(df)) {
    cat(sprintf("\n跳过 %s: 变量不存在\n", bio_name))
    next
  }
  cat(sprintf("\n分析: %s (%s)\n", bio_name, bio_var))
  # 连续变量
  if(is.numeric(df[[bio_var]])) {
    formula <- as.formula(paste0(bio_var, " ~ HCF_type + RIDAGEYR + RIAGENDR"))
    model <- tryCatch(svyglm(formula, design = mec_design), error = function(e) NULL)
    if(!is.null(model)) {
      coefs <- summary(model)$coefficients
      for(type in c("纯生理型", "纯心理型", "身心混合型")) {
        coef_name <- paste0("HCF_type", type)
        if(coef_name %in% rownames(coefs)) {
          biomarker_results <- rbind(biomarker_results, data.frame(
            Biomarker = bio_name,
            HCF_type = type,
            Beta = coefs[coef_name, "Estimate"],
            SE = coefs[coef_name, "Std. Error"],
            P_value = coefs[coef_name, "Pr(>|t|)"],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}
if(nrow(biomarker_results) > 0) {
  # 将HCF_type列的值映射为英文
  biomarker_results$HCF_type_en <- hcf_labels[biomarker_results$HCF_type]
  biomarker_results_out <- data.frame(
    Biomarker = biomarker_results$Biomarker,
    HCF_Type = biomarker_results$HCF_type_en,
    Beta = biomarker_results$Beta,
    SE = biomarker_results$SE,
    P_value = biomarker_results$P_value
  )
  write.csv(biomarker_results_out, file.path(RESULTS_DIR, "paper1_supp_table3_biomarker_P.csv"), row.names = FALSE)
  cat("\n✅ Table 3 (Biomarker gradient) saved\n")
}
# ============================================================================
# 8. 分析4：临床意义分析
# ============================================================================
cat("\n========================================================\n")
cat("分析4：临床意义分析\n")
cat("========================================================\n")
# 创建临床结局
if("HUQ090" %in% names(df)) {
  df$hospitalization <- ifelse(df$HUQ090 %in% c(1,2,3), 1, 0)
  df$emergency <- ifelse(df$HUQ090 == 2, 1, 0)
  df$function_limitation <- ifelse(df$HUQ090 %in% c(1,2), 1, 0)
  mec_design <- update(mec_design,
                       hospitalization = df$hospitalization,
                       emergency = df$emergency,
                       function_limitation = df$function_limitation)
}
clinical_vars <- c(
  "Hospitalization" = "hospitalization",
  "Emergency Visit" = "emergency",
  "Activity Limitation" = "function_limitation",
  "Antihypertensive Medication" = "antihypertensive_med",
  "Antidiabetic Medication" = "antidiabetic_med"
)
clinical_results <- data.frame()
for(i in 1:length(clinical_vars)) {
  outcome_name <- names(clinical_vars)[i]
  outcome_var <- clinical_vars[i]
  if(!outcome_var %in% names(df)) next
  event_count <- sum(df[[outcome_var]] == 1, na.rm = TRUE)
  if(event_count < 10) next
  cat(sprintf("\n分析: %s, 事件数: %d\n", outcome_name, event_count))
  formula <- as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR"))
  model <- tryCatch(svyglm(formula, design = mec_design, family = quasibinomial()),
                    error = function(e) NULL)
  if(!is.null(model)) {
    coefs <- summary(model)$coefficients
    for(type in c("纯生理型", "纯心理型", "身心混合型")) {
      coef_name <- paste0("HCF_type", type)
      if(coef_name %in% rownames(coefs)) {
        est <- coefs[coef_name, "Estimate"]
        se <- coefs[coef_name, "Std. Error"]
        clinical_results <- rbind(clinical_results, data.frame(
          Clinical_Outcome = outcome_name,
          HCF_type = type,
          OR = exp(est),
          CI_lower = exp(est - 1.96 * se),
          CI_upper = exp(est + 1.96 * se),
          P_value = coefs[coef_name, "Pr(>|t|)"],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}
if(nrow(clinical_results) > 0) {
  # 将HCF_type列的值映射为英文
  clinical_results$HCF_type_en <- hcf_labels[clinical_results$HCF_type]
  clinical_results_out <- data.frame(
    Clinical_Outcome = clinical_results$Clinical_Outcome,
    HCF_Type = clinical_results$HCF_type_en,
    OR = clinical_results$OR,
    CI_lower = clinical_results$CI_lower,
    CI_upper = clinical_results$CI_upper,
    P_value = clinical_results$P_value
  )
  write.csv(clinical_results_out, file.path(RESULTS_DIR, "paper1_supp_table4_clinical_P.csv"), row.names = FALSE)
  cat("\n✅ Table 4 (Clinical outcomes) saved\n")
}
# ============================================================================
# 9. 生成补充分析报告
# ============================================================================
cat("\n========================================================\n")
cat("9. 生成补充分析报告\n")
cat("========================================================\n")
report_file <- file.path(LOG_DIR, "10_paper1_supp_report_P.txt")
sink(report_file)
cat("论文1补充分析报告 (P周期)\n")
cat("==========================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、中介分析结果 (Table 1)\n")
print(mediation_results)
cat("\n二、行为效应分层结果 (Table 2)\n")
if(nrow(pa_results) > 0) print(pa_results_out)
cat("\n三、生物标志物梯度 (Table 3)\n")
if(nrow(biomarker_results) > 0) print(head(biomarker_results_out))
cat("\n四、临床意义 (Table 4)\n")
if(nrow(clinical_results) > 0) print(head(clinical_results_out))
cat("\n五、输出文件清单\n")
cat(" - paper1_supp_table1_mediation_P.csv\n")
cat(" - paper1_supp_table2_pa_effects_P.csv\n")
cat(" - paper1_supp_table3_biomarker_P.csv\n")
cat(" - paper1_supp_table4_clinical_P.csv\n")
cat(" - paper1_figure_s1_pa_effects_P.pdf\n")
cat(" - paper1_figure_s1_pa_effects_P.png\n")
sink()
cat("\n ✅ 补充分析报告已保存\n\n")
# ============================================================================
# 10. 保存会话信息
# ============================================================================
cat("10. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "10_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期论文1补充分析会话信息\n")
cat("==================================\n")
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
# 11. 保存R代码副本
# ============================================================================
cat("\n11. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "10_paper1_supplement_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "10_code_list_P.txt")
cat("脚本名称: 10_paper1_supplement_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 12. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期论文1补充分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果已保存至:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "RESULTS_DIR", "LOG_DIR")))
gc()
