#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 14_paper2_deep_analysis_P.R
# 描述: NHANES 2017-2020 (P周期) 论文2深入分析 - 通路深化分析合集
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q2, 假设 H4-H5 (深入验证)
# 对应研究计划: 第六部分 - 敏感性分析
# 对应变量详表: 第五部分 5.2 pathway_proxy_vars.rds
#
# 伦理声明: 使用NHANES公开数据，无需IRB批准
# 数据来源: https://www.cdc.gov/nchs/nhanes/index.htm
# 数据周期: 2017-2020 (P系列)
#
# 随机种子: 20240226 (固定以确保可重复性)
# 最后修改: 2026-02-22
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
  "tidyverse", "survey", "ggplot2", "pROC",
  "reshape2", "gridExtra", "scales", "pheatmap"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# 通路英文标签
pathway_labels <- c(
  "低过激健康型" = "Low-hyperactivation",
  "对抗-枯竭混合型" = "Aversion-Exhaustion",
  "痴固着主导型" = "Perseveration",
  "纯生理过激型" = "Hyperactivation"
)
# 年龄组英文标签
age_labels <- c(
  "18-39岁" = "18-39",
  "40-59岁" = "40-59",
  "60岁以上" = "60+"
)
# 配置路径 - P周期独立目录
PROJECT_ROOT <- "C:/NHANES_Data"
CLEAN_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")
RESULTS_DIR <- file.path(CLEAN_DATA_DIR, "results", "paper2")
LOG_DIR <- file.path(CLEAN_DATA_DIR, "logs")
# 创建结果目录
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
if (!dir.exists(file.path(RESULTS_DIR, "figures"))) {
  dir.create(file.path(RESULTS_DIR, "figures"), recursive = TRUE)
}
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("14_paper2_deep_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 14_paper2_deep_analysis_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 论文2深入分析\n")
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
data_file <- file.path(CLEAN_DATA_DIR, "analysis_dataset_subset_P.rds")
if (!file.exists(data_file)) {
  stop("错误: 找不到 analysis_dataset_subset_P.rds")
}
data_raw <- readRDS(data_file)
cat(sprintf("数据加载成功！\n原始数据维度: %d行 x %d列\n\n", 
            nrow(data_raw), ncol(data_raw)))
# ============================================================================
# 3. 准备四通路变量
# ============================================================================
cat("2. 准备四通路变量...\n")
pathway_vars <- c("avoidance_z", "perseveration_z", "hyperactivation_z", "exhaustion_z")
# 检查变量
missing_vars <- pathway_vars[!pathway_vars %in% names(data_raw)]
if (length(missing_vars) > 0) {
  stop("缺少四通路变量: ", paste(missing_vars, collapse = ", "))
}
# 处理缺失值
for (var in pathway_vars) {
  if (sum(is.na(data_raw[[var]])) > 0) {
    fill_val <- median(data_raw[[var]], na.rm = TRUE)
    data_raw[[var]][is.na(data_raw[[var]])] <- fill_val
    cat(sprintf(" 填充 %s: %d → 中位数 %.3f\n", 
                var, sum(is.na(data_raw[[var]])), fill_val))
  }
}
# ============================================================================
# 4. 重新创建四通路聚类
# ============================================================================
cat("\n3. 重新创建四通路聚类...\n")
set.seed(20240226)
cluster_matrix <- as.matrix(data_raw[, pathway_vars])
km_result <- kmeans(cluster_matrix, centers = 4, nstart = 50, iter.max = 100)
data_raw$cluster_raw <- km_result$cluster
cluster_names <- c(
  "1" = "低过激健康型",
  "2" = "对抗-枯竭混合型",
  "3" = "痴固着主导型",
  "4" = "纯生理过激型"
)
data_raw$pathway_cluster <- factor(
  data_raw$cluster_raw,
  levels = 1:4,
  labels = cluster_names[as.character(1:4)]
)
cat("聚类分布:\n")
print(table(data_raw$pathway_cluster))
# ============================================================================
# 5. 创建行为变量
# ============================================================================
cat("\n4. 创建行为变量...\n")
# 加载健康行为变量
health_file <- file.path(CLEAN_DATA_DIR, "healthbehavior_vars_P.rds")
if(file.exists(health_file)) {
  health <- readRDS(health_file)
  data_raw$pa_meets_guideline <- health$pa_meets_guideline[match(data_raw$SEQN, health$SEQN)]
  data_raw$sleep_adequate <- health$sleep_adequate[match(data_raw$SEQN, health$SEQN)]
  cat(" ✅ 健康行为变量合并完成\n")
}
if (!"sleep_adequate" %in% names(data_raw) && "SLD012" %in% names(data_raw)) {
  data_raw$sleep_adequate <- as.numeric(data_raw$SLD012 >= 7 & data_raw$SLD012 <= 9)
  cat(" ✓ 创建变量: sleep_adequate\n")
}
if (!"pa_meets_guideline" %in% names(data_raw) && "pa_total_min_week" %in% names(data_raw)) {
  data_raw$pa_meets_guideline <- as.numeric(data_raw$pa_total_min_week >= 600)
  cat(" ✓ 创建变量: pa_meets_guideline\n")
}
data_raw <- data_raw %>%
  mutate(
    age_group = case_when(
      RIDAGEYR < 40 ~ "18-39岁",
      RIDAGEYR < 60 ~ "40-59岁",
      TRUE ~ "60岁以上"
    ),
    age_group = factor(age_group, levels = c("18-39岁", "40-59岁", "60岁以上"))
  )
# ============================================================================
# 6. 深化分析1：通路 × α因子 交互
# ============================================================================
cat("\n========================================================\n")
cat("深化分析1：通路 × α因子 交互作用\n")
cat("========================================================\n")
design1 <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINTPRP,
  nest = TRUE,
  data = data_raw
)
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
alpha_names <- c("元认知", "情感调节", "系统协调", "目标效能")
interaction_results <- data.frame()
for (i in seq_along(alpha_vars)) {
  alpha <- alpha_vars[i]
  alpha_name <- alpha_names[i]
  cat("\n【", alpha_name, "(", alpha, ")×通路交互】\n")
  model_int <- svyglm(as.formula(paste0("phq9_total ~ ", alpha, " * pathway_cluster + RIDAGEYR + RIAGENDR")),
                      design = design1)
  int_test <- regTermTest(model_int, as.formula(paste0("~", alpha, ":pathway_cluster")), method = "Wald")
  cat(" 交互项整体检验: F =", round(int_test$Ftest, 3), 
      ", p =", format.pval(int_test$p, digits = 3), "\n")
  coefs <- summary(model_int)$coefficients
  alpha_rows <- grep(paste0("^", alpha, ":"), rownames(coefs))
  for (row in alpha_rows) {
    pathway_name <- gsub(paste0(alpha, ":pathway_cluster"), "", rownames(coefs)[row])
    interaction_results <- rbind(interaction_results, data.frame(
      alpha因子 = alpha_name,
      通路 = pathway_name,
      beta = coefs[row, 1],
      SE = coefs[row, 2],
      p_value = coefs[row, 4],
      CI_lower = coefs[row, 1] - 1.96 * coefs[row, 2],
      CI_upper = coefs[row, 1] + 1.96 * coefs[row, 2],
      stringsAsFactors = FALSE
    ))
    cat(" ", pathway_name, ": beta =", round(coefs[row, 1], 3),
        ", p =", format.pval(coefs[row, 4], digits = 3), "\n")
  }
}
write.csv(interaction_results, file.path(RESULTS_DIR, "deep1_interaction_P.csv"), row.names = FALSE)
cat("\n✅ 已保存: deep1_interaction_P.csv\n")
# ============================================================================
# 7. 深化分析2：通路特异性炎症谱
# ============================================================================
cat("\n========================================================\n")
cat("深化分析2：通路特异性炎症谱\n")
cat("========================================================\n")
design2 <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMECPRP,
  nest = TRUE,
  data = data_raw
)
inflammation_vars <- c(
  "hs_crp_mgl" = "hs-CRP (mg/L)",
  "nlr_ratio" = "中性粒-淋巴比值",
  "white_blood_cells_x10e9_l" = "白细胞计数"
)
inflammation_vars <- inflammation_vars[names(inflammation_vars) %in% names(data_raw)]
if (length(inflammation_vars) > 0) {
  inflam_by_pathway <- svyby(
    make.formula(paste0("~", paste(names(inflammation_vars), collapse = "+"))),
    ~pathway_cluster,
    design2,
    svymean,
    na.rm = TRUE,
    keep.var = TRUE
  )
  cat("\n【各通路炎症指标均值】\n")
  print(inflam_by_pathway)
  inflam_pvals <- data.frame()
  for (var in names(inflammation_vars)) {
    model <- svyglm(as.formula(paste0(var, " ~ pathway_cluster")), design = design2)
    ft <- regTermTest(model, ~pathway_cluster, method = "Wald")
    inflam_pvals <- rbind(inflam_pvals, data.frame(
      指标 = inflammation_vars[var],
      F = ft$Ftest,
      p = ft$p
    ))
  }
  cat("\n【通路间差异检验】\n")
  print(inflam_pvals)
  write.csv(inflam_by_pathway, file.path(RESULTS_DIR, "deep2_inflammation_means_P.csv"), row.names = FALSE)
  write.csv(inflam_pvals, file.path(RESULTS_DIR, "deep2_inflammation_pvals_P.csv"), row.names = FALSE)
  cat("\n✅ 已保存炎症分析结果\n")
}
# ============================================================================
# 8. 深化分析3：通路介导HCF对抑郁的影响
# ============================================================================
cat("\n========================================================\n")
cat("深化分析3：通路介导HCF对抑郁的影响\n")
cat("========================================================\n")
if ("HCF_type" %in% names(data_raw) && "phq9_total" %in% names(data_raw)) {
  data_analysis3 <- data_raw
  hcf_levels <- levels(factor(data_analysis3$HCF_type))
  cat("HCF_type水平:", paste(hcf_levels, collapse = ", "), "\n")
  hcf_map <- c(
    "健康型" = "Healthy",
    "纯生理型" = "Physiological",
    "纯心理型" = "Psychological",
    "身心混合型" = "Mixed"
  )
  hcf_ref <- "健康型"
  hcf_dummies <- c()
  for (level in hcf_levels) {
    if (level != hcf_ref) {
      eng_name <- hcf_map[level]
      var_name <- paste0("HCF_", eng_name)
      data_analysis3[[var_name]] <- as.numeric(data_analysis3$HCF_type == level)
      hcf_dummies <- c(hcf_dummies, var_name)
    }
  }
  cat("参照组: 健康型 (Healthy)\n")
  cat("HCF虚拟变量:", paste(hcf_dummies, collapse = ", "), "\n")
  pathway_levels <- levels(data_analysis3$pathway_cluster)
  cat("通路水平:", paste(pathway_levels, collapse = ", "), "\n")
  pathway_map <- c(
    "低过激健康型" = "Low",
    "对抗-枯竭混合型" = "Avoidant",
    "痴固着主导型" = "Perseverative",
    "纯生理过激型" = "Hyper"
  )
  pathway_ref <- "低过激健康型"
  pathway_dummies <- c()
  for (level in pathway_levels) {
    if (level != pathway_ref) {
      eng_name <- pathway_map[level]
      var_name <- paste0("P_", eng_name)
      data_analysis3[[var_name]] <- as.numeric(data_analysis3$pathway_cluster == level)
      pathway_dummies <- c(pathway_dummies, var_name)
    }
  }
  cat("通路参照组: 低过激健康型 (Low)\n")
  cat("通路虚拟变量:", paste(pathway_dummies, collapse = ", "), "\n")
  design3 <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTINTPRP,
    nest = TRUE,
    data = data_analysis3
  )
  formula_total <- as.formula(paste0("phq9_total ~ ", 
                                      paste(hcf_dummies, collapse = " + "),
                                      " + RIDAGEYR + RIAGENDR"))
  model_total <- svyglm(formula_total, design = design3)
  formula_direct <- as.formula(paste0("phq9_total ~ ",
                                       paste(hcf_dummies, collapse = " + "), " + ",
                                       paste(pathway_dummies, collapse = " + "),
                                       " + RIDAGEYR + RIAGENDR"))
  model_direct <- svyglm(formula_direct, design = design3)
  mediation_results <- data.frame()
  for (hcf in hcf_dummies) {
    total <- coef(model_total)[hcf]
    direct <- coef(model_direct)[hcf]
    indirect <- total - direct
    hcf_name <- gsub("^HCF_", "", hcf)
    mediation_results <- rbind(mediation_results, data.frame(
      HCF类型 = hcf_name,
      总效应 = total,
      直接效应 = direct,
      间接效应 = indirect,
      中介比例 = indirect / total * 100,
      stringsAsFactors = FALSE
    ))
  }
  cat("\n【中介分析结果】\n")
  print(mediation_results)
  write.csv(mediation_results, file.path(RESULTS_DIR, "deep3_mediation_P.csv"), row.names = FALSE)
  cat("\n✅ 已保存: deep3_mediation_P.csv\n")
  sink(file.path(RESULTS_DIR, "deep3_model_details_P.txt"))
  cat("【总效应模型】\n")
  print(summary(model_total))
  cat("\n【直接效应模型】\n")
  print(summary(model_direct))
  sink()
  cat("✅ 已保存: deep3_model_details_P.txt\n")
}
# ============================================================================
# 9. 深化分析4：通路演进风险验证
# ============================================================================
cat("\n========================================================\n")
cat("深化分析4：通路演进风险验证\n")
cat("========================================================\n")
if ("alpha3" %in% names(data_raw)) {
  data_analysis4 <- data_raw %>%
    mutate(
      alpha3_scaled = (alpha3 - min(alpha3, na.rm = TRUE)) /
                      (max(alpha3, na.rm = TRUE) - min(alpha3, na.rm = TRUE)),
      risk_score = 0.3 * avoidance_z + 0.3 * hyperactivation_z + 0.4 * (1 - alpha3_scaled)
    )
  design4_int <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTINTPRP,
    nest = TRUE,
    data = data_analysis4
  )
  design4_mec <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMECPRP,
    nest = TRUE,
    data = data_analysis4
  )
  model_exh <- svyglm(exhaustion_z ~ risk_score + RIDAGEYR + RIAGENDR, design = design4_int)
  cat("\n【风险评分预测枯竭】\n")
  print(coef(summary(model_exh))["risk_score", ])
  if ("hs_crp_mgl" %in% names(data_analysis4)) {
    model_crp <- svyglm(log(hs_crp_mgl + 0.1) ~ risk_score + RIDAGEYR + RIAGENDR,
                        design = design4_mec)
    cat("\n【风险评分预测log(CRP)】\n")
    print(coef(summary(model_crp))["risk_score", ])
  }
  risk_by_path <- svyby(~risk_score, ~pathway_cluster, design4_int, svymean)
  cat("\n【各通路风险评分】\n")
  print(risk_by_path)
  write.csv(risk_by_path, file.path(RESULTS_DIR, "deep4_risk_by_pathway_P.csv"), row.names = FALSE)
  cat("\n✅ 已保存风险验证结果\n")
}
# ============================================================================
# 10. 深化分析5：多周期验证（如果有）
# ============================================================================
cat("\n========================================================\n")
cat("深化分析5：多周期验证\n")
cat("========================================================\n")
if ("SDDSRVYR" %in% names(data_raw) && length(unique(data_raw$SDDSRVYR)) > 1) {
  design5 <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTINTPRP,
    nest = TRUE,
    data = data_raw
  )
  cycle_by_path <- svyby(~pathway_cluster, ~SDDSRVYR, design5, svymean)
  cat("\n【各周期通路分布】\n")
  print(cycle_by_path)
  cycle_test <- svychisq(~pathway_cluster + SDDSRVYR, design5)
  cat("\n【周期×通路独立性检验】\n")
  print(cycle_test)
  write.csv(cycle_by_path, file.path(RESULTS_DIR, "deep5_pathway_by_cycle_P.csv"), row.names = FALSE)
  cat("\n✅ 已保存多周期验证结果\n")
} else {
  cat("\n⚠️ 无多周期数据，跳过深化分析5\n")
}
# ============================================================================
# 12. 炎症轨迹图 (Figure S5)
# ============================================================================
cat("\n========================================================\n")
cat("12. 生成炎症轨迹图: Figure S5\n")
cat("========================================================\n")
if ("hs_crp_mgl" %in% names(data_raw) && "age_group" %in% names(data_raw)) {
  # 过滤掉 pathway_cluster 为 NA 的样本
  data_filtered <- data_raw %>%
    filter(!is.na(pathway_cluster), !is.na(age_group), !is.na(hs_crp_mgl))
  cat(sprintf("过滤前样本量: %d\n", nrow(data_raw)))
  cat(sprintf("过滤后样本量: %d\n", nrow(data_filtered)))
  # 创建设计对象（用过滤后的数据）
  design_mec <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMECPRP,
    nest = TRUE,
    data = data_filtered
  )
  # 计算各通路各年龄组的炎症均值
  inflam_trajectory <- svyby(~hs_crp_mgl, ~pathway_cluster + age_group,
                              design_mec, svymean, na.rm = TRUE)
  cat("\n【各通路炎症随年龄变化】\n")
  # 转换为数据框
  inflam_df <- data.frame(
    通路 = as.character(inflam_trajectory$pathway_cluster),
    年龄组 = as.character(inflam_trajectory$age_group),
    hs_crp = round(inflam_trajectory$hs_crp_mgl, 2),
    SE = round(inflam_trajectory$se.hs_crp_mgl, 2),
    stringsAsFactors = FALSE
  )
  print(inflam_df)
  # 保存数据
  write.csv(inflam_df, file.path(RESULTS_DIR, "table_deep_inflammation_trajectory_P.csv"), row.names = FALSE)
  if (nrow(inflam_df) > 0 && !all(is.na(inflam_df$SE))) {
    # 创建英文标签
    inflam_df$通路_en <- factor(inflam_df$通路,
                                 levels = names(pathway_labels),
                                 labels = pathway_labels)
    inflam_df$年龄组_en <- factor(inflam_df$年龄组,
                                   levels = names(age_labels),
                                   labels = age_labels)
    # 绘制轨迹图
    p <- ggplot(inflam_df, aes(x = 年龄组_en, y = hs_crp, 
                                color = 通路_en, group = 通路_en)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = hs_crp - 1.96 * SE, 
                        ymax = hs_crp + 1.96 * SE), width = 0.2) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "Figure S5. Inflammation Trajectories by Pathway (P Cycle)",
           x = "Age Group", y = "hs-CRP (mg/L)",
           color = "Pathway Type") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
      ) +
      guides(color = guide_legend(nrow = 2, byrow = TRUE))
    # 显示图形
    print(p)
    # 保存PDF
    ggsave(file.path(RESULTS_DIR, "figure_S5_inflammation_trajectory_P.pdf"), 
           p, width = 10, height = 6)
    # 保存PNG
    ggsave(file.path(RESULTS_DIR, "figure_S5_inflammation_trajectory_P.png"), 
           p, width = 10, height = 6, dpi = 300)
    cat(" ✅ Figure S5 saved: figure_S5_inflammation_trajectory_P.pdf and .png\n")
  }
}
# ============================================================================
# 13. 保存会话信息
# ============================================================================
cat("\n13. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "14_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期论文2深入分析会话信息\n")
cat("====================================\n")
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
# 14. 保存R代码副本
# ============================================================================
cat("\n14. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "14_paper2_deep_analysis_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "14_code_list_P.txt")
cat("脚本名称: 14_paper2_deep_analysis_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 15. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期论文2深入分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果已保存至:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
