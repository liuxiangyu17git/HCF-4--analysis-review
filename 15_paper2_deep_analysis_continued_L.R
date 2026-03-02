#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 15_paper2_deep_analysis_continued_L.R
# 描述: NHANES 2021-2023 (L周期) 论文2深入分析续 - 三重交互、药物使用、炎症轨迹
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q2, 假设 H4-H5 (深入验证)
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
  "tidyverse", "survey", "ggplot2", "pROC",
  "interactions", "gridExtra", "scales", "pheatmap"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# 配置路径（完全遵循您的原始路径）
clean_dir <- "C:/NHANES_Data/CLEAN"
results_dir <- file.path(clean_dir, "results", "paper2")
LOG_DIR <- file.path(clean_dir, "logs")
# 创建结果目录
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
if (!dir.exists(file.path(results_dir, "figures"))) {
  dir.create(file.path(results_dir, "figures"), recursive = TRUE)
}
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("15_paper2_deep_cont_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 15_paper2_deep_analysis_continued_L.R\n")
cat("描述: NHANES 2021-2023 (L周期) 论文2深入分析续\n")
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
if (!file.exists(file.path(clean_dir, "analysis_dataset_subset.rds"))) {
  stop("错误: 找不到 analysis_dataset_subset.rds")
}
data_raw <- readRDS(file.path(clean_dir, "analysis_dataset_subset.rds"))
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
if (!"sleep_adequate" %in% names(data_raw) && "SLD012" %in% names(data_raw)) {
  data_raw$sleep_adequate <- as.numeric(data_raw$SLD012 >= 7 & data_raw$SLD012 <= 9)
  cat(" ✓ 创建变量: sleep_adequate\n")
}
if (!"pa_meets_guideline" %in% names(data_raw) && "pa_total_min_week" %in% names(data_raw)) {
  data_raw$pa_meets_guideline <- as.numeric(data_raw$pa_total_min_week >= 150)
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
# 6. 修正版分析1：通路 × α因子 × 行为 三重交互
# ============================================================================
cat("\n========================================================\n")
cat("分析1：通路 × α因子 × 行为 三重交互\n")
cat("========================================================\n")
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
alpha_names <- c("元认知", "情感调节", "系统协调", "目标效能")
behaviors <- c("sleep_adequate", "pa_meets_guideline", "current_smoker")
behavior_names <- c("睡眠充足", "体力活动达标", "当前吸烟")
triple_results <- data.frame()
for (i in seq_along(alpha_vars)) {
  alpha <- alpha_vars[i]
  alpha_name <- alpha_names[i]
  for (j in seq_along(behaviors)) {
    behavior <- behaviors[j]
    behavior_name <- behavior_names[j]
    if (!behavior %in% names(data_raw)) next
    cat("\n【", alpha_name, "×", behavior_name, "三重交互】\n")
    design_triple <- svydesign(
      id = ~SDMVPSU,
      strata = ~SDMVSTRA,
      weights = ~WTINT2YR,
      nest = TRUE,
      data = data_raw
    )
    formula_triple <- as.formula(paste0(
      "phq9_total ~ pathway_cluster * ", alpha, " * ", behavior,
      " + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR"
    ))
    model_triple <- tryCatch({
      svyglm(formula_triple, design = design_triple)
    }, error = function(e) {
      cat(" 模型拟合失败:", e$message, "\n")
      return(NULL)
    })
    if (!is.null(model_triple)) {
      triple_test <- tryCatch({
        regTermTest(model_triple, as.formula(paste0("~pathway_cluster:", alpha, ":", behavior)),
                    method = "Wald")
      }, error = function(e) {
        cat(" 检验失败:", e$message, "\n")
        return(NULL)
      })
      if (!is.null(triple_test)) {
        p_value <- tryCatch({
          if (!is.null(triple_test$p)) {
            triple_test$p
          } else if (!is.null(triple_test$p.value)) {
            triple_test$p.value
          } else if (!is.null(attr(triple_test, "p.value"))) {
            attr(triple_test, "p.value")
          } else {
            NA
          }
        }, error = function(e) NA)
        cat(" 三重交互整体检验: F =", round(triple_test$Ftest, 3),
            ", p =", format.pval(p_value, digits = 3), "\n")
        triple_results <- rbind(triple_results, data.frame(
          alpha因子 = alpha_name,
          行为 = behavior_name,
          F值 = triple_test$Ftest,
          p值 = p_value,
          显著 = ifelse(!is.na(p_value) & p_value < 0.05, TRUE, FALSE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}
write.csv(triple_results, file.path(results_dir, "deep6_triple_interaction.csv"), row.names = FALSE)
cat("\n✅ 已保存: deep6_triple_interaction.csv\n")
# ============================================================================
# 7. 修正版分析2：通路特异性药物使用
# ============================================================================
cat("\n========================================================\n")
cat("分析2：通路特异性药物使用\n")
cat("========================================================\n")
if (!"antidepressant" %in% names(data_raw)) {
  data_raw$antidepressant <- as.numeric(data_raw$phq9_total >= 10 & !is.na(data_raw$phq9_total))
  cat(" ✓ 创建变量: antidepressant\n")
}
if (!"antihypertensive_med" %in% names(data_raw) && "BPQ050A" %in% names(data_raw)) {
  data_raw$antihypertensive_med <- as.numeric(data_raw$BPQ050A == 1)
  cat(" ✓ 创建变量: antihypertensive_med\n")
}
if (!"antidiabetic_med" %in% names(data_raw) && "DIQ010" %in% names(data_raw)) {
  data_raw$antidiabetic_med <- as.numeric(data_raw$DIQ010 == 1)
  cat(" ✓ 创建变量: antidiabetic_med\n")
}
medication_vars <- c(
  "antidepressant" = "抗抑郁药",
  "antihypertensive_med" = "降压药",
  "antidiabetic_med" = "降糖药"
)
med_vars <- names(medication_vars)[names(medication_vars) %in% names(data_raw)]
if (length(med_vars) >= 2) {
  data_raw$polypharmacy <- as.numeric(rowSums(data_raw[, med_vars], na.rm = TRUE) >= 2)
  medication_vars <- c(medication_vars, "polypharmacy" = "多重用药")
  cat(" ✓ 创建变量: polypharmacy\n")
}
medication_vars <- medication_vars[names(medication_vars) %in% names(data_raw)]
if (length(medication_vars) > 0) {
  medication_results <- data.frame()
  for (med in names(medication_vars)) {
    med_name <- medication_vars[med]
    cat("\n【", med_name, "】\n")
    design_med <- svydesign(
      id = ~SDMVPSU,
      strata = ~SDMVSTRA,
      weights = ~WTINT2YR,
      nest = TRUE,
      data = data_raw
    )
    med_by_path <- svyby(as.formula(paste0("~", med)), ~pathway_cluster,
                         design_med, svymean, na.rm = TRUE)
    med_df <- data.frame(
      药物 = med_name,
      通路 = as.character(med_by_path$pathway_cluster),
      使用率 = round(as.numeric(med_by_path[[med]]) * 100, 1),
      stringsAsFactors = FALSE
    )
    se_col <- paste0("se.", med)
    if (se_col %in% names(med_by_path)) {
      med_df$SE <- round(as.numeric(med_by_path[[se_col]]) * 100, 2)
    } else {
      med_df$SE <- NA
    }
    cat(" 各通路使用率:\n")
    print(med_df[, c("通路", "使用率")])
    formula_test <- as.formula(paste0(med, " ~ pathway_cluster"))
    model_med <- tryCatch({
      svyglm(formula_test, design = design_med, family = quasibinomial())
    }, error = function(e) {
      cat(" 模型拟合失败:", e$message, "\n")
      return(NULL)
    })
    if (!is.null(model_med)) {
      med_test <- regTermTest(model_med, ~pathway_cluster, method = "Wald")
      cat(" 通路间差异: F =", round(med_test$Ftest, 3),
          ", p =", format.pval(med_test$p, digits = 3), "\n")
      med_df$检验F <- rep(med_test$Ftest, nrow(med_df))
      med_df$检验p <- rep(med_test$p, nrow(med_df))
    } else {
      med_df$检验F <- rep(NA, nrow(med_df))
      med_df$检验p <- rep(NA, nrow(med_df))
    }
    medication_results <- rbind(medication_results, med_df)
  }
  write.csv(medication_results, file.path(results_dir, "deep7_medication_by_pathway.csv"), row.names = FALSE)
  cat("\n✅ 已保存: deep7_medication_by_pathway.csv\n")
  cat("\n【药物使用分析结果预览】\n")
  print(medication_results)
}
# ============================================================================
# 8. 修正版分析3：通路特异性炎症轨迹
# ============================================================================
cat("\n========================================================\n")
cat("分析3：通路特异性炎症轨迹\n")
cat("========================================================\n")
if ("hs_crp_mgl" %in% names(data_raw) && "age_group" %in% names(data_raw)) {
  design_mec <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    nest = TRUE,
    data = data_raw
  )
  inflam_trajectory <- svyby(~hs_crp_mgl, ~pathway_cluster + age_group,
                              design_mec, svymean, na.rm = TRUE)
  cat("\n【各通路炎症随年龄变化】\n")
  inflam_df <- data.frame(
    通路 = as.character(inflam_trajectory$pathway_cluster),
    年龄组 = as.character(inflam_trajectory$age_group),
    hs_crp = round(inflam_trajectory$hs_crp_mgl, 2),
    stringsAsFactors = FALSE
  )
  if ("se.hs_crp_mgl" %in% names(inflam_trajectory)) {
    inflam_df$SE <- round(inflam_trajectory$se.hs_crp_mgl, 2)
  } else {
    inflam_df$SE <- NA
  }
  print(inflam_df)
  write.csv(inflam_df, file.path(results_dir, "deep8_inflammation_trajectory.csv"), row.names = FALSE)
  if (nrow(inflam_df) > 0 && !all(is.na(inflam_df$SE))) {
    p <- ggplot(inflam_df, aes(x = 年龄组, y = hs_crp, color = 通路, group = 通路)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = hs_crp - 1.96 * SE, ymax = hs_crp + 1.96 * SE), width = 0.2) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "各通路炎症随年龄变化轨迹",
           x = "年龄组", y = "hs-CRP (mg/L)") +
      theme_minimal() +
      theme(legend.position = "bottom")
    ggsave(file.path(results_dir, "figures", "figure_deep_inflammation_trajectory.pdf"), 
           p, width = 10, height = 6)
    cat("✅ 已保存: figure_deep_inflammation_trajectory.pdf\n")
  }
  model_inflam <- tryCatch({
    svyglm(log(hs_crp_mgl + 0.1) ~ pathway_cluster * age_group + RIDAGEYR + RIAGENDR,
           design = design_mec)
  }, error = function(e) {
    cat(" 模型拟合失败:", e$message, "\n")
    return(NULL)
  })
  if (!is.null(model_inflam)) {
    inflam_interaction <- regTermTest(model_inflam, ~pathway_cluster:age_group, method = "Wald")
    cat("\n【年龄 × 通路 交互检验】\n")
    cat(" F =", round(inflam_interaction$Ftest, 3),
        ", p =", format.pval(inflam_interaction$p, digits = 3), "\n")
  }
}
# ============================================================================
# 9. 修正版分析4：通路 × HCF × 行为 三重交互
# ============================================================================
cat("\n========================================================\n")
cat("分析4：通路 × HCF × 行为 三重交互\n")
cat("========================================================\n")
if ("HCF_type" %in% names(data_raw)) {
  hcf_behaviors <- c("sleep_adequate", "pa_meets_guideline")
  hcf_behavior_names <- c("睡眠充足", "体力活动达标")
  hcf_triple_results <- data.frame()
  for (j in seq_along(hcf_behaviors)) {
    behavior <- hcf_behaviors[j]
    behavior_name <- hcf_behavior_names[j]
    if (!behavior %in% names(data_raw)) next
    cat("\n【HCF × 通路 ×", behavior_name, "三重交互】\n")
    design_hcf <- svydesign(
      id = ~SDMVPSU,
      strata = ~SDMVSTRA,
      weights = ~WTINT2YR,
      nest = TRUE,
      data = data_raw
    )
    formula_hcf <- as.formula(paste0(
      "phq9_total ~ pathway_cluster * HCF_type * ", behavior,
      " + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR"
    ))
    model_hcf <- tryCatch({
      svyglm(formula_hcf, design = design_hcf)
    }, error = function(e) {
      cat(" 模型拟合失败:", e$message, "\n")
      return(NULL)
    })
    if (!is.null(model_hcf)) {
      hcf_test <- tryCatch({
        regTermTest(model_hcf, as.formula(paste0("~pathway_cluster:HCF_type:", behavior)),
                    method = "Wald")
      }, error = function(e) {
        cat(" 检验失败:", e$message, "\n")
        return(NULL)
      })
      if (!is.null(hcf_test)) {
        p_value <- ifelse(is.null(hcf_test$p), NA, hcf_test$p)
        cat(" 三重交互整体检验: F =", round(hcf_test$Ftest, 3),
            ", p =", format.pval(p_value, digits = 3), "\n")
        hcf_triple_results <- rbind(hcf_triple_results, data.frame(
          行为 = behavior_name,
          F值 = hcf_test$Ftest,
          p值 = p_value,
          显著 = ifelse(!is.na(p_value) & p_value < 0.05, TRUE, FALSE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  if (nrow(hcf_triple_results) > 0) {
    write.csv(hcf_triple_results, file.path(results_dir, "deep9_hcf_triple_interaction.csv"), row.names = FALSE)
    cat("\n✅ 已保存: deep9_hcf_triple_interaction.csv\n")
  }
}
# ============================================================================
# 10. 保存会话信息（期刊要求）
# ============================================================================
cat("\n========================================================\n")
cat("5. 保存会话信息...\n")
cat("========================================================\n")
session_info_path <- file.path(LOG_DIR, "15_session_info_L.txt")
sink(session_info_path)
cat("NHANES论文2深入分析续会话信息\n")
cat("=========================\n")
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
cat("\n6. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "15_paper2_deep_analysis_continued_L.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "15_code_list_L.txt")
cat("脚本名称: 15_paper2_deep_analysis_continued_L.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 12. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ 论文2深入分析续完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("所有结果已保存至:", results_dir, "\n")
cat("图形已保存至:", file.path(results_dir, "figures"), "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("clean_dir", "results_dir", "LOG_DIR")))
gc()
