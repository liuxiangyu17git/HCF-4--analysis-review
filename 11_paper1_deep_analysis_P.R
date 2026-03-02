#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 11_paper1_deep_analysis_P.R
# 描述: NHANES 2017-2020 (P周期) paper1深入分析 - C层中介、D层深化、行为交互
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
# 设置随机种子（期刊要求）
set.seed(20240226)
# 加载必要包
required_packages <- c("survey", "dplyr", "emmeans", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# HCF分型英文标签
hcf_labels <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
# 配置路径 - P周期独立目录
PROJECT_ROOT <- "C:/NHANES_Data"
CLEAN_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")
RESULTS_DIR <- file.path(CLEAN_DATA_DIR, "results")
LOG_DIR <- file.path(CLEAN_DATA_DIR, "logs")
# 创建结果目录
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("11_paper1_deep_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 11_paper1_deep_analysis_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 论文1深入分析\n")
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
# 2. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
data_file <- file.path(CLEAN_DATA_DIR, "final_analysis_dataset_P.rds")
if(!file.exists(data_file)) {
  stop("错误: final_analysis_dataset_P.rds 不存在")
}
data_full <- readRDS(data_file)
data <- data_full %>% filter(in_analysis == 1)
cat(sprintf(" 分析样本: %d人\n\n", nrow(data)))
# ============================================================================
# 3. 创建设计对象（P周期使用WTMECPRP）
# ============================================================================
cat("2. 创建设计对象...\n")
mec_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMECPRP,
  nest = TRUE,
  data = data
)
cat(" ✅ 设计对象创建成功\n\n")
# ============================================================================
# 4. 分析1：C层中介模型（重复验证）
# ============================================================================
cat("\n========================================================\n")
cat("分析1：C层中介模型（重复验证）\n")
cat("========================================================\n")
# 步骤1：A层对抑郁的总效应
model_total <- svyglm(phq9_total ~ HCF_A + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                      design = mec_design)
# 步骤2：A层对C层的影响
model_c <- svyglm(HCF_C ~ HCF_A + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                  family = quasibinomial(), design = mec_design)
# 步骤3：A层 + C层对抑郁的影响
model_outcome <- svyglm(phq9_total ~ HCF_A + HCF_C + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                        design = mec_design)
# 提取系数
total_effect <- coef(model_total)["HCF_A"]
a_effect <- coef(model_c)["HCF_A"]
b_effect <- coef(model_outcome)["HCF_C"]
direct_effect <- coef(model_outcome)["HCF_A"]
# 计算中介效应
indirect_effect <- a_effect * b_effect
prop_mediated <- indirect_effect / total_effect * 100
cat("\n=== C层中介分析结果 ===\n")
cat(sprintf("总效应: %.3f\n", total_effect))
cat(sprintf("直接效应: %.3f\n", direct_effect))
cat(sprintf("间接效应: %.3f\n", indirect_effect))
cat(sprintf("中介比例: %.1f%%\n", prop_mediated))
# 保存结果
mediation_results <- data.frame(
  effect = c("Total", "Direct", "Indirect", "Prop_mediated"),
  estimate = c(total_effect, direct_effect, indirect_effect, prop_mediated)
)
write.csv(mediation_results, file.path(RESULTS_DIR, "paper1_deep_mediation_P.csv"), row.names = FALSE)
cat("✅ 结果已保存: paper1_deep_mediation_P.csv\n\n")
# ============================================================================
# 5. 分析2：D层深化分析
# ============================================================================
cat("\n========================================================\n")
cat("分析2：D层深化分析\n")
cat("========================================================\n")
# 检查D层相关变量
if("HCF_D" %in% names(data) && "DMDEDUC2" %in% names(data)) {
  # D层与教育水平的关联
  model_d1 <- svyglm(HCF_D ~ DMDEDUC2 + RIDAGEYR + RIAGENDR,
                     family = quasibinomial(), design = mec_design)
  cat("\nD层与教育水平的关联:\n")
  coefs_d1 <- summary(model_d1)$coefficients
  print(coefs_d1)
  # D层与静息心率的关联
  if("BPXOPLS1" %in% names(data)) {
    model_d2 <- svyglm(HCF_D ~ BPXOPLS1 + RIDAGEYR + RIAGENDR,
                       family = quasibinomial(), design = mec_design)
    cat("\nD层与静息心率的关联:\n")
    coefs_d2 <- summary(model_d2)$coefficients
    print(coefs_d2)
  }
  # D层预测抑郁
  model_d3 <- svyglm(phq9_total ~ HCF_D + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                     design = mec_design)
  cat("\nD层预测抑郁:\n")
  coefs_d3 <- summary(model_d3)$coefficients
  print(coefs_d3["HCF_D", ])
  # 保存D层分析结果
  d_layer_results <- list(
    education = as.data.frame(coefs_d1),
    heart_rate = if(exists("coefs_d2")) as.data.frame(coefs_d2) else NULL,
    depression = as.data.frame(t(coefs_d3["HCF_D", ]))
  )
  saveRDS(d_layer_results, file.path(RESULTS_DIR, "paper1_deep_Dlayer_P.rds"))
  cat("✅ D层分析结果已保存: paper1_deep_Dlayer_P.rds\n\n")
}
# ============================================================================
# 6. 分析3：分型特异性行为效应深入分析
# ============================================================================
cat("\n========================================================\n")
cat("分析3：分型特异性行为效应深入分析\n")
cat("========================================================\n")
# 定义行为因素
behaviors <- list(
  "体力活动达标" = "pa_meets_guideline",
  "睡眠充足" = "sleep_adequate",
  "从不吸烟" = "never_smoker",
  "非重度饮酒" = "not_heavy_drinker"
)
# 加载健康行为变量
health_file <- file.path(CLEAN_DATA_DIR, "healthbehavior_vars_P.rds")
if(file.exists(health_file)) {
  health <- readRDS(health_file)
  data$pa_meets_guideline <- health$pa_meets_guideline[match(data$SEQN, health$SEQN)]
  data$sleep_adequate <- health$sleep_adequate[match(data$SEQN, health$SEQN)]
}
# 创建非重度饮酒变量
if("heavy_drinker" %in% names(data)) {
  data$not_heavy_drinker <- ifelse(data$heavy_drinker == 0, 1, 0)
  mec_design <- update(mec_design, not_heavy_drinker = data$not_heavy_drinker)
}
behavior_results <- data.frame()
for(behavior_name in names(behaviors)) {
  behavior_var <- behaviors[[behavior_name]]
  if(!behavior_var %in% names(data)) {
    next
  }
  cat(sprintf("\n=== 分析: %s ===\n", behavior_name))
  hcf_types <- c("健康型", "纯生理型", "纯心理型", "身心混合型")
  for(hcf_type in hcf_types) {
    subset_design <- tryCatch({
      subset(mec_design, HCF_type == hcf_type)
    }, error = function(e) NULL)
    if(is.null(subset_design)) next
    type_data <- data[data$HCF_type == hcf_type, ]
    var_data <- type_data[[behavior_var]]
    n_yes <- sum(var_data == 1, na.rm = TRUE)
    n_no <- sum(var_data == 0, na.rm = TRUE)
    if(n_yes < 5 || n_no < 5) {
      cat(sprintf(" 跳过 %s: 样本量不足 (1=%d, 0=%d)\n", hcf_type, n_yes, n_no))
      next
    }
    formula <- as.formula(paste0("phq9_total ~ ", behavior_var, " + RIDAGEYR + RIAGENDR"))
    model <- tryCatch({
      svyglm(formula, design = subset_design)
    }, error = function(e) NULL)
    if(!is.null(model)) {
      coef_matrix <- summary(model)$coefficients
      if(behavior_var %in% rownames(coef_matrix)) {
        behavior_results <- rbind(behavior_results, data.frame(
          行为因素 = behavior_name,
          HCF_type = hcf_type,
          beta = coef_matrix[behavior_var, "Estimate"],
          SE = coef_matrix[behavior_var, "Std. Error"],
          p_value = coef_matrix[behavior_var, "Pr(>|t|)"],
          N_yes = n_yes,
          N_no = n_no,
          stringsAsFactors = FALSE
        ))
        cat(sprintf(" ✅ %s: beta=%.3f, p=%.4f\n",
                    hcf_type,
                    coef_matrix[behavior_var, "Estimate"],
                    coef_matrix[behavior_var, "Pr(>|t|)"]))
      }
    }
  }
}
if(nrow(behavior_results) > 0) {
  write.csv(behavior_results, file.path(RESULTS_DIR, "paper1_deep_behavior_P.csv"), row.names = FALSE)
  cat("\n✅ 行为因素深入分析结果已保存: paper1_deep_behavior_P.csv\n")
}
# ============================================================================
# 7. 分析4：复合中介模型（A层 → C层 → D层 → 抑郁）
# ============================================================================
cat("\n========================================================\n")
cat("分析4：复合中介模型（A层 → C层 → D层 → 抑郁）\n")
cat("========================================================\n")
if(all(c("HCF_A", "HCF_C", "HCF_D", "phq9_total") %in% names(data))) {
  # 模型1：A层 → 抑郁 (总效应)
  model_total2 <- svyglm(phq9_total ~ HCF_A + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                         design = mec_design)
  # 模型2：A层 → C层
  model_a_c <- svyglm(HCF_C ~ HCF_A + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                      family = quasibinomial(), design = mec_design)
  # 模型3：A层 + C层 → D层
  model_ac_d <- svyglm(HCF_D ~ HCF_A + HCF_C + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                       family = quasibinomial(), design = mec_design)
  # 模型4：A层 + C层 + D层 → 抑郁
  model_acd_dep <- svyglm(phq9_total ~ HCF_A + HCF_C + HCF_D + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                          design = mec_design)
  # 提取系数
  a_to_c <- coef(model_a_c)["HCF_A"]
  c_to_d <- coef(model_ac_d)["HCF_C"]
  d_to_dep <- coef(model_acd_dep)["HCF_D"]
  # 计算间接效应
  indirect_c <- a_to_c * coef(model_acd_dep)["HCF_C"]
  indirect_d <- a_to_c * c_to_d * d_to_dep
  total_indirect <- indirect_c + indirect_d
  total_effect2 <- coef(model_total2)["HCF_A"]
  direct_effect2 <- coef(model_acd_dep)["HCF_A"]
  cat("\n=== 复合中介模型结果 ===\n")
  cat(sprintf("总效应: %.3f\n", total_effect2))
  cat(sprintf("直接效应: %.3f\n", direct_effect2))
  cat(sprintf("间接效应 (C层): %.3f\n", indirect_c))
  cat(sprintf("间接效应 (C→D链): %.3f\n", indirect_d))
  cat(sprintf("总间接效应: %.3f\n", total_indirect))
  cat(sprintf("中介比例: %.1f%%\n", total_indirect / total_effect2 * 100))
  # 保存复合中介结果
  composite_mediation <- data.frame(
    path = c("Total", "Direct", "Indirect_C", "Indirect_CD", "Total_Indirect", "Proportion"),
    estimate = c(total_effect2, direct_effect2, indirect_c, indirect_d, total_indirect,
                 total_indirect / total_effect2 * 100)
  )
  write.csv(composite_mediation, file.path(RESULTS_DIR, "paper1_deep_composite_mediation_P.csv"), row.names = FALSE)
  cat("✅ 复合中介结果已保存: paper1_deep_composite_mediation_P.csv\n")
}
# ============================================================================
# 8. 分析5：行为因素×HCF分型交互作用可视化
# ============================================================================
cat("\n========================================================\n")
cat("分析5：行为因素×HCF分型交互作用可视化\n")
cat("========================================================\n")
if(exists("behavior_results") && nrow(behavior_results) > 0) {
# 筛选体力活动结果
pa_plot_data <- behavior_results %>%
  filter(行为因素 == "体力活动达标") %>%
  mutate(
    CI_lower = beta - 1.96 * SE,
    CI_upper = beta + 1.96 * SE
  )
if(nrow(pa_plot_data) > 0) {
  p1 <- ggplot(pa_plot_data, aes(x = HCF_type, y = beta, fill = HCF_type)) +
    geom_col() +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_x_discrete(labels = hcf_labels) +  # ← 添加
    scale_fill_manual(values = c(
      "健康型" = "#4DAF4A",
      "纯生理型" = "#377EB8",
      "纯心理型" = "#FF7F00",
      "身心混合型" = "#E41A1C"
    ), labels = hcf_labels) +  # ← 添加labels
    labs(title = "Figure S2. Physical Activity Effects on Depression by HCF Type",
         x = "", y = "β Coefficient (95% CI)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  # 同时保存PDF和PNG
ggsave(file.path(RESULTS_DIR, "paper1_figure_deep_pa_effects_P.pdf"), p1, width = 8, height = 6)
ggsave(file.path(RESULTS_DIR, "paper1_figure_deep_pa_effects_P.png"), p1, width = 8, height = 6, dpi = 300)
  cat("✅ 体力活动效应图已保存: PDF和PNG格式\n")
}
# 筛选睡眠结果
sleep_plot_data <- behavior_results %>%
  filter(行为因素 == "睡眠充足") %>%
  mutate(
    CI_lower = beta - 1.96 * SE,
    CI_upper = beta + 1.96 * SE
  )
if(nrow(sleep_plot_data) > 0) {
  p2 <- ggplot(sleep_plot_data, aes(x = HCF_type, y = beta, fill = HCF_type)) +
    geom_col() +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_x_discrete(labels = hcf_labels) +  # ← 添加
    scale_fill_manual(values = c(
      "健康型" = "#4DAF4A",
      "纯生理型" = "#377EB8",
      "纯心理型" = "#FF7F00",
      "身心混合型" = "#E41A1C"
    ), labels = hcf_labels) +  # ← 添加labels
    labs(title = "Figure S3. Sleep Effects on Depression by HCF Type",
         x = "", y = "β Coefficient (95% CI)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  # 同时保存PDF和PNG
  ggsave(file.path(RESULTS_DIR, "paper1_figure_deep_sleep_effects_P.pdf"), p2, width = 8, height = 6)
ggsave(file.path(RESULTS_DIR, "paper1_figure_deep_sleep_effects_P.png"), p2, width = 8, height = 6, dpi = 300)
  cat("✅ 睡眠效应图已保存: PDF和PNG格式\n")
}
}
# ============================================================================
# 9. 生成深入分析报告
# ============================================================================
cat("\n========================================================\n")
cat("9. 生成深入分析报告\n")
cat("========================================================\n")
report_file <- file.path(LOG_DIR, "11_paper1_deep_report_P.txt")
sink(report_file)
cat("论文1深入分析报告 (P周期)\n")
cat("==========================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、C层中介分析结果\n")
print(mediation_results)
cat("\n二、行为因素深入分析结果（摘要）\n")
if(exists("behavior_results") && nrow(behavior_results) > 0) {
  print(behavior_results %>% select(行为因素, HCF_type, beta, p_value))
}
cat("\n三、复合中介模型结果\n")
if(exists("composite_mediation")) {
  print(composite_mediation)
}
cat("\n四、输出文件清单\n")
cat(" - paper1_deep_mediation_P.csv\n")
cat(" - paper1_deep_Dlayer_P.rds\n")
cat(" - paper1_deep_behavior_P.csv\n")
cat(" - paper1_deep_composite_mediation_P.csv\n")
cat(" - paper1_figure_deep_pa_effects_P.pdf\n")
cat(" - paper1_figure_deep_sleep_effects_P.pdf\n")
sink()
cat(" ✅ 深入分析报告已保存\n\n")
# ============================================================================
# 10. 保存会话信息（期刊要求）
# ============================================================================
cat("10. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "11_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期论文1深入分析会话信息\n")
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
# 11. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n11. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "11_paper1_deep_analysis_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "11_code_list_P.txt")
cat("脚本名称: 11_paper1_deep_analysis_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 12. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期论文1深入分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果已保存至:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "RESULTS_DIR", "LOG_DIR")))
gc()
