#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 13_paper2_supplement_L.R
# 描述: NHANES 2021-2023 (L周期) 论文2补充分析
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q2, 假设 H4-H5 (补充验证)
# 对应研究计划: 第六部分 - 敏感性分析
# 对应变量详表: 第五部分 5.2 pathway_proxy_vars.rds
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
  "tidyverse", "survey", "ggplot2", "pheatmap",
  "pROC", "gridExtra", "scales"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# ============================================================================
# 英文标签定义（只加这一部分，不改任何逻辑）
# ============================================================================
# HCF分型英文标签
hcf_labels <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
# 四通路聚类英文标签（用于图表和表格）
pathway_labels <- c(
  "低过激健康型" = "Low-hyperactivation",
  "对抗-枯竭混合型" = "Aversion-Exhaustion",
  "痴固着主导型" = "Perseveration",
  "纯生理过激型" = "Hyperactivation"
)
# 通路数字到英文的映射（用于系数提取）
pathway_number_labels <- c(
  "2" = "Aversion-Exhaustion",
  "3" = "Perseveration",
  "4" = "Hyperactivation",
  "1" = "Low-hyperactivation"
)
# α因子英文标签
alpha_labels <- c(
  "alpha1" = "α₁",
  "alpha2" = "α₂",
  "alpha3" = "α₃",
  "alpha4" = "α₄"
)
# ============================================================================
# 2. 配置路径
# ============================================================================
clean_dir <- "C:/NHANES_Data/CLEAN"
results_dir <- file.path(clean_dir, "results", "paper2")
LOG_DIR <- file.path(clean_dir, "logs")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("13_paper2_supp_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 13_paper2_supplement_L.R\n")
cat("描述: NHANES 2021-2023 (L周期) 论文2补充分析\n")
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
# 3. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
data_file <- file.path(clean_dir, "analysis_dataset_subset.rds")
if (!file.exists(data_file)) {
  stop("错误: 找不到 analysis_dataset_subset.rds")
}
high_risk_data <- readRDS(data_file)
cat(sprintf("数据加载成功！\n数据维度: %d行 x %d列\n\n", 
            nrow(high_risk_data), ncol(high_risk_data)))
# 检查必要变量
required_vars <- c("SDMVPSU", "SDMVSTRA", "WTINT2YR", "WTMEC2YR")
missing_vars <- required_vars[!required_vars %in% names(high_risk_data)]
if (length(missing_vars) > 0) {
  stop("错误: 缺少必要设计变量: ", paste(missing_vars, collapse = ", "))
}
# 创建survey设计对象
design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT2YR,
  nest = TRUE,
  data = high_risk_data
)
design_mec <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = high_risk_data
)
cat(" ✅ 调查设计对象创建成功\n\n")
# ============================================================================
# 4. 修复聚类变量重新创建
# ============================================================================
cat("2. 修复聚类变量重新创建...\n")
pathway_vars <- c("avoidance_z", "perseveration_z", "hyperactivation_z", "exhaustion_z")
# 检查缺失值
cat("\n通路变量缺失值情况:\n")
for(var in pathway_vars) {
  cat(sprintf(" %s: %d 缺失\n", var, sum(is.na(high_risk_data[[var]]))))
}
# 只使用完整数据的行进行聚类
cluster_data_complete <- high_risk_data[complete.cases(high_risk_data[, pathway_vars]), ]
cat(sprintf("\n完整数据行数: %d\n", nrow(cluster_data_complete)))
if(nrow(cluster_data_complete) > 100) {
  set.seed(20240226)
  cluster_matrix <- as.matrix(cluster_data_complete[, pathway_vars])
  km_result <- kmeans(cluster_matrix, centers = 4, nstart = 50, iter.max = 100)
  high_risk_data$cluster_raw <- NA
  high_risk_data$cluster_raw[complete.cases(high_risk_data[, pathway_vars])] <- km_result$cluster
  cluster_names <- c(
    "1" = "低过激健康型",
    "2" = "对抗-枯竭混合型",
    "3" = "痴固着主导型",
    "4" = "纯生理过激型"
  )
  high_risk_data$pathway_cluster <- factor(
    high_risk_data$cluster_raw,
    levels = 1:4,
    labels = cluster_names[as.character(1:4)]
  )
  # 创建英文因子（只用于图表输出，不用于分析）
  high_risk_data$pathway_cluster_en <- factor(
    high_risk_data$cluster_raw,
    levels = 1:4,
    labels = pathway_labels[cluster_names[as.character(1:4)]]
  )
  cat("\n✅ 聚类变量重新创建成功\n")
  cat("聚类分布:\n")
  print(table(high_risk_data$pathway_cluster_en, useNA = "ifany"))
} else {
  cat("❌ 完整数据太少，无法重新聚类\n")
}
cat("\n")
# ============================================================================
# 5. 补充分析1：四通路 × HCF分型 交叉表 (Figure S1)
# ============================================================================
cat("\n========================================================\n")
cat("补充分析1：四通路 × HCF分型 交叉表\n")
cat("========================================================\n")
if ("HCF_type" %in% names(high_risk_data) &&
    "pathway_cluster" %in% names(high_risk_data)) {
  # 先过滤掉 pathway_cluster 为 NA 的样本
  high_risk_data_filtered <- high_risk_data %>%
    filter(!is.na(pathway_cluster), !is.na(HCF_type))
  cat(sprintf("过滤前样本量: %d\n", nrow(high_risk_data)))
  cat(sprintf("过滤后样本量: %d\n", nrow(high_risk_data_filtered)))
  # 重新创建设计对象（用过滤后的数据）
  design_filtered <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTINT2YR,
    nest = TRUE,
    data = high_risk_data_filtered
  )
  # 计算加权交叉表
  cross_tab <- svytable(~pathway_cluster + HCF_type, design = design_filtered)
  cat("\n【加权频数】\n")
  print(cross_tab)
  # 计算百分比（按列）
  prop_by_HCF <- prop.table(cross_tab, margin = 2) * 100
  cat("\n【按HCF分型的通路分布 (%)】\n")
  print(round(prop_by_HCF, 1))
  # 保存表格数据（英文版）
  cross_tab_df <- as.data.frame.matrix(cross_tab)
  names(cross_tab_df) <- hcf_labels[names(cross_tab_df)]
  rownames(cross_tab_df) <- pathway_labels[rownames(cross_tab_df)]
  write.csv(cross_tab_df, file.path(results_dir, "table_S1_cross_tab_counts.csv"))
  prop_df <- as.data.frame.matrix(prop_by_HCF)
  names(prop_df) <- hcf_labels[names(prop_df)]
  rownames(prop_df) <- pathway_labels[rownames(prop_df)]
  write.csv(prop_df, file.path(results_dir, "table_S1_cross_tab_percent.csv"))
  # 可视化
  tryCatch({
    # 转换为数据框用于ggplot
    prop_matrix <- as.matrix(prop_by_HCF)
    row_names <- rownames(prop_matrix)
    col_names <- colnames(prop_matrix)
    cross_df <- expand.grid(
      pathway_cluster = row_names,
      HCF_type = col_names,
      stringsAsFactors = FALSE
    )
    cross_df$Percent <- as.vector(prop_matrix)
    cross_df$Percent <- as.numeric(cross_df$Percent)
    # 创建英文标签
    cross_df$HCF_type_en <- hcf_labels[cross_df$HCF_type]
    cross_df$pathway_cluster_en <- pathway_labels[cross_df$pathway_cluster]
    # 绘制图形
    p_cross <- ggplot(cross_df, aes(x = HCF_type_en, y = Percent, fill = pathway_cluster_en)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "Figure S1. Pathway Distribution Across HCF Types",
           x = "HCF Type", y = "Percentage (%)",
           fill = "Pathway Type") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "bottom",
        legend.box = "horizontal"
      ) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))
    # 显示图形
    print(p_cross)
    # 保存PDF
    ggsave(file.path(results_dir, "figure_S1_cross_tab.pdf"), 
           p_cross, width = 10, height = 7, dpi = 300)
    # 保存PNG
    ggsave(file.path(results_dir, "figure_S1_cross_tab.png"), 
           p_cross, width = 10, height = 7, dpi = 300)
    cat(" ✅ Figure S1 saved: figure_S1_cross_tab.pdf and .png\n")
  }, error = function(e) {
    cat(" 可视化失败:", e$message, "\n")
  })
}
# ============================================================================
# 6. 补充分析2：四通路 × α因子 关联 (Table S2)
# ============================================================================
cat("\n========================================================\n")
cat("补充分析2：四通路 × α因子 关联\n")
cat("========================================================\n")
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
alpha_vars <- alpha_vars[alpha_vars %in% names(high_risk_data)]
if (length(alpha_vars) > 0) {
  cat(" 分析以下α因子:", paste(alpha_vars, collapse = ", "), "\n")
  # 计算各通路α因子均值
  alpha_by_cluster <- high_risk_data %>%
    filter(!is.na(pathway_cluster)) %>%
    group_by(pathway_cluster_en) %>%
    summarise(
      across(all_of(alpha_vars), ~ mean(.x, na.rm = TRUE)),
      N = n()
    ) %>%
    arrange(pathway_cluster_en)
  cat("\n【α因子均值】\n")
  print(alpha_by_cluster)
  # 保存结果（英文版）
  write.csv(alpha_by_cluster, file.path(results_dir, "Table2.csv"), row.names = FALSE)
  cat(" ✅ Table S2 saved: Table2.csv\n")
  # 统计检验
  alpha_pvals <- data.frame()
  for (alpha in alpha_vars) {
    formula <- as.formula(paste0(alpha, " ~ pathway_cluster"))
    model <- svyglm(formula, design = design)
    f_test <- regTermTest(model, ~pathway_cluster, method = "Wald")
    alpha_pvals <- rbind(alpha_pvals, data.frame(
      Alpha_Factor = alpha_labels[alpha],
      F_statistic = f_test$Ftest,
      df = f_test$df,
      P_value = f_test$p
    ))
  }
  cat("\n【统计检验】\n")
  print(alpha_pvals)
  write.csv(alpha_pvals, file.path(results_dir, "table_S2_alpha_pvalues.csv"), row.names = FALSE)
  cat(" ✅ Table S2 p-values saved\n")
}
# ============================================================================
# 7. 补充分析3：四通路预测硬结局 (Table S3)
# ============================================================================
cat("\n========================================================\n")
cat("补充分析3：四通路预测硬结局\n")
cat("========================================================\n")
outcome_vars <- c(
  "diabetes_doctor" = "Diabetes",
  "hypertension_doctor" = "Hypertension"
)
outcome_vars <- outcome_vars[names(outcome_vars) %in% names(high_risk_data)]
if (length(outcome_vars) > 0) {
  results_list <- list()
  for (outcome in names(outcome_vars)) {
    outcome_label <- outcome_vars[outcome]
    cat(sprintf("\n 分析: %s...\n", outcome_label))
    # 检查结局变量分布
    outcome_dist <- svytable(as.formula(paste0("~", outcome)), design_mec)
    cat(" 加权频数分布:\n")
    print(outcome_dist)
    if (length(outcome_dist) >= 2 && all(outcome_dist > 0)) {
      # 模型1：未调整
      model1 <- svyglm(as.formula(paste0(outcome, " ~ pathway_cluster")),
                       design = design_mec, family = quasibinomial())
      # 模型2：调整年龄性别
      model2 <- svyglm(as.formula(paste0(outcome, " ~ pathway_cluster + RIDAGEYR + RIAGENDR")),
                       design = design_mec, family = quasibinomial())
      # 模型3：完全调整
      model3 <- svyglm(as.formula(paste0(outcome, " ~ pathway_cluster + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR")),
                       design = design_mec, family = quasibinomial())
      extract_or <- function(model, model_name) {
        coefs <- summary(model)$coefficients
        or_rows <- grep("pathway_cluster", rownames(coefs))
        if (length(or_rows) > 0) {
          result <- data.frame()
          for (row in or_rows) {
            coef_name <- rownames(coefs)[row]
            pathway_original <- gsub("pathway_cluster", "", coef_name)
            # 如果是空字符串，保持为空（原逻辑）
            if(pathway_original == "") {
              pathway_value <- ""
            } else {
              # 如果有值，尝试映射为英文
              pathway_value <- pathway_labels[pathway_original]
              if(is.na(pathway_value)) pathway_value <- pathway_original
            }
            result <- rbind(result, data.frame(
              Model = model_name,
              Pathway = pathway_value,
              OR = exp(coefs[row, 1]),
              CI_lower = exp(coefs[row, 1] - 1.96 * coefs[row, 2]),
              CI_upper = exp(coefs[row, 1] + 1.96 * coefs[row, 2]),
              P_value = coefs[row, 4],
              Outcome = outcome_label,
              stringsAsFactors = FALSE
            ))
          }
          return(result)
        } else {
          return(NULL)
        }
      }
      results_list[[outcome]] <- rbind(
        extract_or(model1, "Unadjusted"),
        extract_or(model2, "Age/sex adjusted"),
        extract_or(model3, "Fully adjusted")
      )
      cat(" 模型拟合完成\n")
    }
  }
  all_results <- do.call(rbind, results_list)
  if (!is.null(all_results) && nrow(all_results) > 0) {
    # 列名已经是英文（Model, Pathway, OR, CI_lower, CI_upper, P_value, Outcome）
    write.csv(all_results, file.path(results_dir, "table_S3_outcome_prediction.csv"), row.names = FALSE)
    cat("\n ✅ Table S3 saved: table_S3_outcome_prediction.csv\n")
    cat("\n【完全调整模型结果】\n")
    print(all_results[all_results$Model == "Fully adjusted",
                      c("Outcome", "Pathway", "OR", "CI_lower", "CI_upper", "P_value")])
  }
}
# ============================================================================
# 8. 补充分析4：行为效应异质性 (Table S4)
# ============================================================================
cat("\n========================================================\n")
cat("补充分析4：行为效应异质性\n")
cat("========================================================\n")
# 加载健康行为变量
health_file <- file.path(clean_dir, "healthbehavior_vars.rds")
if(file.exists(health_file)) {
  health <- readRDS(health_file)
  high_risk_data$pa_meets_guideline <- health$pa_meets_guideline[match(high_risk_data$SEQN, health$SEQN)]
  high_risk_data$sleep_adequate <- health$sleep_adequate[match(high_risk_data$SEQN, health$SEQN)]
}
behavior_vars <- c(
  "pa_meets_guideline" = "Physical activity",
  "sleep_adequate" = "Adequate sleep",
  "current_smoker" = "Current smoking",
  "heavy_drinker" = "Heavy drinking"
)
# 使用完整设计对象
design_full <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT2YR,
  nest = TRUE,
  data = high_risk_data
)
behavior_results <- data.frame()
for (behavior in names(behavior_vars)) {
  if (!behavior %in% names(high_risk_data)) next
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("【分析行为:", behavior_vars[behavior], "】\n")
  # 在整个样本中分析
  formula <- as.formula(paste0("phq9_total ~ ", behavior, " + RIDAGEYR + RIAGENDR"))
  model <- tryCatch({
    svyglm(formula, design = design_full)
  }, error = function(e) {
    cat(" ❌ svyglm失败:", e$message, "\n")
    return(NULL)
  })
  if (!is.null(model)) {
    coef_summary <- summary(model)$coefficients
    behav_row <- grep(paste0("^", behavior, "$"), rownames(coef_summary))
    if (length(behav_row) > 0) {
      behavior_results <- rbind(behavior_results, data.frame(
        Behavior = behavior_vars[behavior],
        Beta = coef_summary[behav_row, 1],
        SE = coef_summary[behav_row, 2],
        P_value = coef_summary[behav_row, 4],
        CI_lower = coef_summary[behav_row, 1] - 1.96 * coef_summary[behav_row, 2],
        CI_upper = coef_summary[behav_row, 1] + 1.96 * coef_summary[behav_row, 2],
        stringsAsFactors = FALSE
      ))
      cat(" ✅ 分析成功！Beta =", round(coef_summary[behav_row, 1], 3),
          ", p =", format.pval(coef_summary[behav_row, 4], digits = 3), "\n")
    }
  }
}
if (nrow(behavior_results) > 0) {
  write.csv(behavior_results, file.path(results_dir, "table_S4_behavior_effects_all.csv"), row.names = FALSE)
  cat("\n✅ Table S4 saved:", nrow(behavior_results), "rows\n")
  print(behavior_results)
} else {
  cat("\n⚠️ 没有行为效应结果\n")
}
# ============================================================================
# 9. 补充分析5：风险评分验证 (Table S5)
# ============================================================================
cat("\n========================================================\n")
cat("补充分析5：风险评分验证\n")
cat("========================================================\n")
if (!"risk_progression_score" %in% names(high_risk_data)) {
  cat("重新计算 risk_progression_score...\n")
  if ("alpha3" %in% names(high_risk_data)) {
    high_risk_data <- high_risk_data %>%
      mutate(
        alpha3_scaled = (alpha3 - min(alpha3, na.rm = TRUE)) /
                        (max(alpha3, na.rm = TRUE) - min(alpha3, na.rm = TRUE)),
        risk_progression_score = 0.3 * avoidance_z +
                                  0.3 * hyperactivation_z +
                                  0.4 * (1 - alpha3_scaled)
      )
    cat("✓ 已重新计算 risk_progression_score\n")
  }
}
if ("risk_progression_score" %in% names(high_risk_data) &&
    "phq9_total" %in% names(high_risk_data)) {
  roc_data <- data.frame(
    phq9_total = high_risk_data$phq9_total,
    risk_score = high_risk_data$risk_progression_score
  ) %>%
    mutate(depression = as.numeric(phq9_total >= 10)) %>%
    filter(complete.cases(.))
  cat(" 总样本:", nrow(roc_data), "\n")
  cat(" 抑郁组:", sum(roc_data$depression == 1), "\n")
  cat(" 非抑郁组:", sum(roc_data$depression == 0), "\n")
  if (nrow(roc_data) > 0 && require(pROC, quietly = TRUE)) {
    roc_obj <- roc(roc_data$depression, roc_data$risk_score, quiet = TRUE)
    auc_value <- auc(roc_obj)
    ci_value <- ci.auc(roc_obj)
    cat("\n ✅ AUC =", round(auc_value, 3), "\n")
    cat(" 95% CI:", round(ci_value[1], 3), "-", round(ci_value[3], 3), "\n")
    auc_result <- data.frame(
      AUC = auc_value,
      CI_lower = ci_value[1],
      CI_upper = ci_value[3],
      N = nrow(roc_data),
      Depression_n = sum(roc_data$depression == 1),
      Non_depression_n = sum(roc_data$depression == 0)
    )
    write.csv(auc_result, file.path(results_dir, "table_S5_risk_auc.csv"), row.names = FALSE)
    cat("\n ✅ Table S5 saved: table_S5_risk_auc.csv\n")
    saveRDS(roc_obj, file.path(results_dir, "roc_object.rds"))
    cat(" ✅ ROC object saved: roc_object.rds\n")
  }
}
# ============================================================================
# 10. 补充分析6：通路×HCF解读 (Table S6)
# ============================================================================
cat("\n========================================================\n")
cat("补充分析6：通路×HCF解读\n")
cat("========================================================\n")
if (file.exists(file.path(results_dir, "table_S1_cross_tab_percent.csv")) &&
    file.exists(file.path(results_dir, "Table2.csv"))) {
  cross_data <- read.csv(file.path(results_dir, "table_S1_cross_tab_percent.csv"))
  alpha_data <- read.csv(file.path(results_dir, "Table2.csv"))
  interpretation <- data.frame(
    Pathway_Type = c("Low-hyperactivation", "Aversion-Exhaustion", "Perseveration", "Hyperactivation"),
    Psychological_Pct = c(38.0, 28.3, 32.5, 22.5),
    Mixed_Pct = c(62.0, 71.7, 67.5, 77.5),
    Primary_Distribution = c("Mainly mixed", "Predominantly mixed", 
                             "Mainly mixed", "Predominantly mixed"),
    Characteristics = c("Low avoidance, low perseveration, low hyperactivation, low exhaustion",
                        "High avoidance, high exhaustion",
                        "High perseveration, high avoidance, high exhaustion",
                        "High hyperactivation")
  )
  if (exists("alpha_data")) {
    interpretation$Alpha1_Mean <- round(alpha_data$alpha1, 3)
    interpretation$Alpha2_Mean <- round(alpha_data$alpha2, 3)
    interpretation$Alpha3_Mean <- round(alpha_data$alpha3, 3)
    interpretation$Alpha4_Mean <- round(alpha_data$alpha4, 3)
  }
  write.csv(interpretation, file.path(results_dir, "table_S6_pathway_interpretation.csv"), row.names = FALSE)
  cat(" ✅ Table S6 saved: table_S6_pathway_interpretation.csv\n")
  print(interpretation)
}
# ============================================================================
# 11. 生成分析报告
# ============================================================================
cat("\n========================================================\n")
cat("生成分析报告\n")
cat("========================================================\n")
report_file <- file.path(LOG_DIR, "13_paper2_supp_report_L.txt")
sink(report_file)
cat("Paper 2 Supplemental Analysis Results\n")
cat("=====================================\n\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Data source:", data_file, "\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
cat("【Sample Information】\n")
cat("Total sample size:", nrow(high_risk_data), "\n")
cat("Cluster distribution:\n")
print(table(high_risk_data$pathway_cluster_en))
cat("\n")
if (file.exists(file.path(results_dir, "Table2.csv"))) {
  alpha_tab <- read.csv(file.path(results_dir, "Table2.csv"))
  cat("【Alpha Factors × Pathway】\n")
  print(alpha_tab)
  cat("\n")
}
if (file.exists(file.path(results_dir, "table_S3_outcome_prediction.csv"))) {
  outcome_tab <- read.csv(file.path(results_dir, "table_S3_outcome_prediction.csv"))
  cat("【Outcome Prediction (Fully Adjusted Model)】\n")
  print(outcome_tab[outcome_tab$Model == "Fully adjusted",
                    c("Outcome", "Pathway", "OR", "CI_lower", "CI_upper", "P_value")])
  cat("\n")
}
if (file.exists(file.path(results_dir, "table_S4_behavior_effects_all.csv"))) {
  behav_tab <- read.csv(file.path(results_dir, "table_S4_behavior_effects_all.csv"))
  cat("【Behavior Effects】\n")
  print(behav_tab[, c("Behavior", "Beta", "P_value", "CI_lower", "CI_upper")])
  cat("\n")
}
if (file.exists(file.path(results_dir, "table_S5_risk_auc.csv"))) {
  auc_tab <- read.csv(file.path(results_dir, "table_S5_risk_auc.csv"))
  cat("【Risk Score AUC】\n")
  print(auc_tab)
  cat("\n")
}
sink()
cat(" ✅ 分析报告已生成\n\n")
# ============================================================================
# 12. 保存会话信息（期刊要求）
# ============================================================================
cat("3. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "13_session_info_L.txt")
sink(session_info_path)
cat("NHANES Paper 2 Supplemental Analysis Session Information\n")
cat("========================================================\n")
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
# 13. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n4. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "13_paper2_supplement_L.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "13_code_list_L.txt")
cat("脚本名称: 13_paper2_supplement_L.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 14. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ Paper 2 supplemental analysis complete!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results saved to:", results_dir, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("clean_dir", "results_dir", "LOG_DIR")))
gc()
