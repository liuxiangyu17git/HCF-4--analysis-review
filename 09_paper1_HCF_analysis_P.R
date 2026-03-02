#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 09_paper1_HCF_analysis_P.R
# 描述: NHANES 2017-2020 (P周期) 论文1 - HCF宏观分型验证分析
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q1, 假设 H1-H3
# 对应研究计划: 第五部分 图1、表1、表2
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
required_packages <- c("tidyverse", "survey", "knitr", "kableExtra", "corrplot", "igraph")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# 尝试加载qgraph（如果可用）
if(!require(qgraph, quietly = TRUE)) {
  install.packages("qgraph")
  library(qgraph)
}
# 配置路径 - P周期独立目录
PROJECT_ROOT <- "C:/NHANES_Data"
CLEAN_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")
RESULTS_DIR <- file.path(CLEAN_DATA_DIR, "results")
LOG_DIR <- file.path(CLEAN_DATA_DIR, "logs")
# 创建结果目录
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("09_paper1_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 09_paper1_HCF_analysis_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 论文1 - HCF宏观分型验证\n")
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
# 2. 定义英文标签
# ============================================================================
# HCF分型英文映射
hcf_names_en <- c(
  "健康型" = "Healthy",
  "纯生理型" = "Pure Physiological",
  "纯心理型" = "Psychological",
  "身心混合型" = "Psychosomatic Mixed"
)
# ============================================================================
# 3. 加载最终数据集
# ============================================================================
cat("1. 加载最终数据集...\n")
data_file <- file.path(CLEAN_DATA_DIR, "final_analysis_dataset_P.rds")
if(!file.exists(data_file)) {
  stop("错误: final_analysis_dataset_P.rds不存在! 请先运行08_final_merge_P.R")
}
data_full <- readRDS(data_file)
cat(sprintf(" final_analysis_dataset_P.rds: %d行, %d列\n", nrow(data_full), ncol(data_full)))
# 筛选分析人群（成人非孕妇）
data <- data_full %>% filter(in_analysis == 1)
cat(sprintf(" 分析人群: %d行 (%.1f%%)\n\n", nrow(data), nrow(data)/nrow(data_full)*100))
# ============================================================================
# 4. 创建调查设计对象（P周期使用WTMECPRP）
# ============================================================================
cat("2. 创建调查设计对象...\n")
mec_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMECPRP,
  nest = TRUE,
  data = data
)
cat(" ✅ MEC调查设计对象创建完成\n")
cat(" PSU数:", length(unique(data$SDMVPSU)), "\n")
cat(" 层数:", length(unique(data$SDMVSTRA)), "\n\n")
# ============================================================================
# 5. 检查HCF分型变量
# ============================================================================
cat("3. 检查HCF分型变量...\n")
hcf_vars <- c("HCF_A", "HCF_B", "HCF_C", "HCF_D", "HCF_type")
hcf_present <- hcf_vars[hcf_vars %in% names(data)]
cat(" 存在的HCF变量:\n")
for(var in hcf_present) {
  n_missing <- sum(is.na(data[[var]]))
  pct_missing <- mean(is.na(data[[var]])) * 100
  cat(sprintf("  %s: 非缺失 %d (%.1f%%), 缺失 %d (%.1f%%)\n",
              var,
              sum(!is.na(data[[var]])),
              100 - pct_missing,
              n_missing, pct_missing))
}
# HCF_type分布（未加权）
cat("\n HCF_type分布（未加权）:\n")
hcf_table <- table(data$HCF_type, useNA = "ifany")
print(hcf_table)
cat("\n")
# ============================================================================
# 6. 表1：HCF分型的人口学特征（加权）
# ============================================================================
cat("4. 生成表1：HCF分型的人口学特征（加权）...\n")
# 定义要分析的人口学变量
demo_vars <- list(
  "Age_years" = "RIDAGEYR",
  "Gender_pct" = "RIAGENDR",
  "Race_pct" = "RIDRETH3",
  "Education_pct" = "DMDEDUC2",
  "Poverty_ratio" = "INDFMPIR",
  "Marital_status_pct" = "DMDMARTZ"
)
# 创建表1的数据框（英文列名）
table1 <- data.frame(
  Variable = character(),
  HCF_type = character(),
  N = numeric(),
  Statistic = character(),
  stringsAsFactors = FALSE
)
hcf_levels <- levels(factor(data$HCF_type))
for(hcf in hcf_levels) {
  if(is.na(hcf)) next
  subset_design <- subset(mec_design, HCF_type == hcf)
  n_hcf <- sum(data$HCF_type == hcf, na.rm = TRUE)
  # 样本量
  table1 <- rbind(table1, data.frame(
    Variable = "Sample_size",
    HCF_type = hcf,
    N = n_hcf,
    Statistic = sprintf("%d (%.1f%%)", n_hcf, n_hcf/nrow(data)*100)
  ))
  # 年龄
  age_n <- sum(!is.na(data$RIDAGEYR[data$HCF_type == hcf]))
  age_mean <- svymean(~RIDAGEYR, subset_design, na.rm = TRUE)
  age_se <- SE(age_mean)
  table1 <- rbind(table1, data.frame(
    Variable = "Age_years",
    HCF_type = hcf,
    N = age_n,
    Statistic = sprintf("%.1f (%.2f)", coef(age_mean), age_se)
  ))
  # 性别
  gender_tab <- svytable(~RIAGENDR, subset_design)
  gender_prop <- prop.table(gender_tab) * 100
  # 检查是否有男性（RIAGENDR=1）和女性（RIAGENDR=2）
  if("1" %in% names(gender_tab)) {
    table1 <- rbind(table1, data.frame(
      Variable = "Gender_Male",
      HCF_type = hcf,
      N = round(gender_tab["1"]),
      Statistic = sprintf("%.1f%%", gender_prop["1"])
    ))
  }
  if("2" %in% names(gender_tab)) {
    table1 <- rbind(table1, data.frame(
      Variable = "Gender_Female",
      HCF_type = hcf,
      N = round(gender_tab["2"]),
      Statistic = sprintf("%.1f%%", gender_prop["2"])
    ))
  }
  # 种族（简化版本，避免表格过长）
  if("RIDRETH3" %in% names(data)) {
    race_tab <- svytable(~RIDRETH3, subset_design)
    race_prop <- prop.table(race_tab) * 100
    # 只汇总主要种族
    white_idx <- which(names(race_tab) == "3")  # 非西班牙裔白人
    black_idx <- which(names(race_tab) == "4")  # 非西班牙裔黑人
    if(length(white_idx) > 0) {
      table1 <- rbind(table1, data.frame(
        Variable = "Race_NonHispanic_White",
        HCF_type = hcf,
        N = round(race_tab[white_idx]),
        Statistic = sprintf("%.1f%%", race_prop[white_idx])
      ))
    }
    if(length(black_idx) > 0) {
      table1 <- rbind(table1, data.frame(
        Variable = "Race_NonHispanic_Black",
        HCF_type = hcf,
        N = round(race_tab[black_idx]),
        Statistic = sprintf("%.1f%%", race_prop[black_idx])
      ))
    }
  }
  # 教育程度（大学及以上）
  if("DMDEDUC2" %in% names(data)) {
    college_prop <- svymean(~I(DMDEDUC2 >= 4), subset_design, na.rm = TRUE)
    college_n <- sum(data$DMDEDUC2 >= 4 & data$HCF_type == hcf, na.rm = TRUE)
    table1 <- rbind(table1, data.frame(
      Variable = "College_or_above",
      HCF_type = hcf,
      N = college_n,
      Statistic = sprintf("%.1f%%", 100 * coef(college_prop))
    ))
  }
  # 贫困（收入<1.0）
  if("INDFMPIR" %in% names(data)) {
    poverty_prop <- svymean(~I(INDFMPIR < 1), subset_design, na.rm = TRUE)
    poverty_n <- sum(data$INDFMPIR < 1 & data$HCF_type == hcf, na.rm = TRUE)
    table1 <- rbind(table1, data.frame(
      Variable = "Poverty",
      HCF_type = hcf,
      N = poverty_n,
      Statistic = sprintf("%.1f%%", 100 * coef(poverty_prop))
    ))
  }
}
# 将HCF_type列的值映射为英文
table1$HCF_type <- hcf_names_en[table1$HCF_type]
write.csv(table1, file.path(RESULTS_DIR, "paper1_table1_demographics_P.csv"), row.names = FALSE)
cat(" ✅ 表1已保存: paper1_table1_demographics_P.csv\n\n")
# ============================================================================
# 7. 表2：HCF分型与心理结局的关联
# ============================================================================
cat("5. 生成表2：HCF分型与心理结局的关联...\n")
# 定义心理结局
mental_outcomes <- list(
  "PHQ9_total" = "phq9_total",
  "Depression_diagnosis" = "depression",
  "Poor_self_rated_health" = "poor_self_rated_health"
)
# 检查变量是否存在
for(outcome in names(mental_outcomes)) {
  var_name <- mental_outcomes[[outcome]]
  if(!var_name %in% names(data)) {
    cat(sprintf(" 警告: %s (%s) 不存在\n", outcome, var_name))
  }
}
# 创建结果表（英文列名）
table2 <- data.frame(
  Outcome = character(),
  Model = character(),
  HCF_type = character(),
  beta_OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
# 以健康型为参照
data$HCF_type <- relevel(factor(data$HCF_type), ref = "健康型")
mec_design <- update(mec_design, HCF_type = factor(data$HCF_type))
# 对每个结局进行分析
for(outcome_name in names(mental_outcomes)) {
  outcome_var <- mental_outcomes[[outcome_name]]
  if(!outcome_var %in% names(data)) next
  cat(sprintf(" 分析: %s\n", outcome_name))
  # 模型：调整年龄、性别
  if(outcome_var == "phq9_total") {
    # 连续结局用线性回归
    model <- svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR")),
                    design = mec_design)
    # 提取系数
    coefs <- summary(model)$coefficients
    hcf_levels <- grep("HCF_type", rownames(coefs), value = TRUE)
    for(hcf in hcf_levels) {
      hcf_name <- gsub("HCF_type", "", hcf)
      ci <- confint(model)[hcf, ]
      table2 <- rbind(table2, data.frame(
        Outcome = outcome_name,
        Model = "Model1_age_sex",
        HCF_type = hcf_name,
        beta_OR = coefs[hcf, "Estimate"],
        CI_lower = ci[1],
        CI_upper = ci[2],
        P_value = coefs[hcf, "Pr(>|t|)"]
      ))
    }
  } else {
    # 二分类结局用逻辑回归
    model <- svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR")),
                    design = mec_design, family = quasibinomial())
    # 提取OR和CI
    coefs <- summary(model)$coefficients
    hcf_levels <- grep("HCF_type", rownames(coefs), value = TRUE)
    for(hcf in hcf_levels) {
      hcf_name <- gsub("HCF_type", "", hcf)
      or <- exp(coefs[hcf, "Estimate"])
      ci <- exp(confint(model)[hcf, ])
      table2 <- rbind(table2, data.frame(
        Outcome = outcome_name,
        Model = "Model1_age_sex",
        HCF_type = hcf_name,
        beta_OR = or,
        CI_lower = ci[1],
        CI_upper = ci[2],
        P_value = coefs[hcf, "Pr(>|t|)"]
      ))
    }
  }
}
# 将HCF_type列的值映射为英文
table2$HCF_type <- hcf_names_en[table2$HCF_type]
write.csv(table2, file.path(RESULTS_DIR, "paper1_table2_mental_outcomes_P.csv"), row.names = FALSE)
cat(" ✅ 表2已保存: paper1_table2_mental_outcomes_P.csv\n\n")
# ============================================================================
# 8. 表3：HCF分型与生理结局的关联
# ============================================================================
cat("6. 生成表3：HCF分型与生理结局的关联...\n")
# 创建生理结局变量
if("DIQ010" %in% names(data)) {
  data$diabetes <- case_when(
    data$DIQ010 == 1 ~ 1,
    data$DIQ010 == 2 ~ 0,
    TRUE ~ NA_real_
  )
  cat(" ✅ 创建 diabetes 变量\n")
}
# P周期使用BPQ020（降压药）
if("BPQ020" %in% names(data)) {
  data$hypertension <- case_when(
    data$BPQ020 == 1 ~ 1,
    data$BPQ020 == 2 ~ 0,
    TRUE ~ NA_real_
  )
  cat(" ✅ 创建 hypertension 变量 (来自BPQ020)\n")
} else if("BPQ150" %in% names(data)) {
  data$hypertension <- case_when(
    data$BPQ150 == 1 ~ 1,
    data$BPQ150 == 2 ~ 0,
    TRUE ~ NA_real_
  )
  cat(" ✅ 创建 hypertension 变量 (来自BPQ150)\n")
}
# 更新设计对象
mec_design <- update(mec_design,
                     diabetes = data$diabetes,
                     hypertension = data$hypertension)
# 定义生理结局
physical_outcomes <- list(
  "Diabetes" = "diabetes",
  "Hypertension" = "hypertension",
  "CVD" = "any_cvd",
  "CKD" = "ckd_flag"
)
# 检查变量是否存在
for(outcome in names(physical_outcomes)) {
  var_name <- physical_outcomes[[outcome]]
  if(!var_name %in% names(data)) {
    cat(sprintf(" 警告: %s (%s) 不存在\n", outcome, var_name))
  } else {
    # 检查变量分布
    cat(sprintf(" %s: 0=%d, 1=%d, NA=%d\n",
                outcome,
                sum(data[[var_name]] == 0, na.rm = TRUE),
                sum(data[[var_name]] == 1, na.rm = TRUE),
                sum(is.na(data[[var_name]]))))
  }
}
# 创建结果表（英文列名）
table3 <- data.frame(
  Outcome = character(),
  Model = character(),
  HCF_type = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
# 对每个生理结局进行分析
for(outcome_name in names(physical_outcomes)) {
  outcome_var <- physical_outcomes[[outcome_name]]
  if(!outcome_var %in% names(data)) next
  cat(sprintf("\n 分析: %s\n", outcome_name))
  # 检查结局变量是否有足够的变异
  outcome_values <- data[[outcome_var]]
  if(sum(outcome_values == 1, na.rm = TRUE) < 10) {
    cat(sprintf("  跳过: 事件数太少 (%d)\n",
                sum(outcome_values == 1, na.rm = TRUE)))
    next
  }
  # 模型：调整年龄、性别
  model <- tryCatch({
    svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR")),
           design = mec_design, family = quasibinomial())
  }, error = function(e) {
    cat(sprintf("  模型错误: %s\n", e$message))
    return(NULL)
  })
  if(!is.null(model)) {
    coefs <- summary(model)$coefficients
    hcf_levels <- grep("HCF_type", rownames(coefs), value = TRUE)
    for(hcf in hcf_levels) {
      hcf_name <- gsub("HCF_type", "", hcf)
      or <- exp(coefs[hcf, "Estimate"])
      ci <- tryCatch(exp(confint(model)[hcf, ]), error = function(e) c(NA, NA))
      table3 <- rbind(table3, data.frame(
        Outcome = outcome_name,
        Model = "Model1_age_sex",
        HCF_type = hcf_name,
        OR = or,
        CI_lower = ci[1],
        CI_upper = ci[2],
        P_value = coefs[hcf, "Pr(>|t|)"]
      ))
    }
  }
}
if(nrow(table3) > 0) {
  # 将HCF_type列的值映射为英文
  table3$HCF_type <- hcf_names_en[table3$HCF_type]
  write.csv(table3, file.path(RESULTS_DIR, "paper1_table3_physical_outcomes_P.csv"), row.names = FALSE)
  cat(" ✅ 表3已保存: paper1_table3_physical_outcomes_P.csv\n\n")
}
# ============================================================================
# 9. 表4：HCF分型与身心复合结局
# ============================================================================
cat("\n7. 生成表4：HCF分型与身心复合结局...\n")
if(all(c("depression", "diabetes", "hypertension") %in% names(data))) {
  data$depression_diabetes <- case_when(
    data$depression == 1 & data$diabetes == 1 ~ 1,
    data$depression == 0 | data$diabetes == 0 ~ 0,
    TRUE ~ NA_real_
  )
  data$depression_hypertension <- case_when(
    data$depression == 1 & data$hypertension == 1 ~ 1,
    data$depression == 0 | data$hypertension == 0 ~ 0,
    TRUE ~ NA_real_
  )
  mec_design <- update(mec_design,
                       depression_diabetes = data$depression_diabetes,
                       depression_hypertension = data$depression_hypertension)
  composite_outcomes <- list(
    "Depression_Diabetes" = "depression_diabetes",
    "Depression_Hypertension" = "depression_hypertension"
  )
  table4 <- data.frame(
    Outcome = character(),
    Model = character(),
    HCF_type = character(),
    OR = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  for(outcome_name in names(composite_outcomes)) {
    outcome_var <- composite_outcomes[[outcome_name]]
    cat(sprintf(" 分析: %s\n", outcome_name))
    event_count <- sum(data[[outcome_var]] == 1, na.rm = TRUE)
    if(event_count < 10) {
      cat(sprintf("  跳过: 事件数太少 (%d)\n", event_count))
      next
    }
    model <- tryCatch({
      svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR")),
             design = mec_design, family = quasibinomial())
    }, error = function(e) {
      cat(sprintf("  模型错误: %s\n", e$message))
      return(NULL)
    })
    if(!is.null(model)) {
      coefs <- summary(model)$coefficients
      hcf_levels <- grep("HCF_type", rownames(coefs), value = TRUE)
      for(hcf in hcf_levels) {
        hcf_name <- gsub("HCF_type", "", hcf)
        or <- exp(coefs[hcf, "Estimate"])
        ci <- tryCatch(exp(confint(model)[hcf, ]), error = function(e) c(NA, NA))
        table4 <- rbind(table4, data.frame(
          Outcome = outcome_name,
          Model = "Model1_age_sex",
          HCF_type = hcf_name,
          OR = or,
          CI_lower = ci[1],
          CI_upper = ci[2],
          P_value = coefs[hcf, "Pr(>|t|)"]
        ))
      }
    }
  }
  if(nrow(table4) > 0) {
    # 将HCF_type列的值映射为英文
    table4$HCF_type <- hcf_names_en[table4$HCF_type]
    write.csv(table4, file.path(RESULTS_DIR, "Table4_composite_P.csv"), row.names = FALSE)
    cat(" ✅ 表4已保存: Table4_composite_P.csv\n")
  }
}
# ============================================================================
# 10. 图1：四型生理-心理网络异质性可视化（P周期修正版 - 不含alt_ui_l）
# ============================================================================
cat("8. 生成图1：四型生理-心理网络异质性（P周期）...\n")
# 关闭所有图形设备
graphics.off()
# 设置输出路径
pdf_path <- file.path(RESULTS_DIR, "Figure3_P.pdf")
png_path <- file.path(RESULTS_DIR, "Figure3_P.png")
# 检查并删除旧文件
if(file.exists(pdf_path)) file.remove(pdf_path)
if(file.exists(png_path)) file.remove(png_path)
if(require(qgraph)) {
  # 定义节点 - P周期没有有效的alt_ui_l数据，所以只包含9个节点
  physio_nodes <- c("BMXBMI", "LBXHSCRP", "BPXOPLS1", "egfr_ckdepi_2021")
  psycho_nodes <- c("DPQ010", "DPQ040", "DPQ070", "DPQ090", "HUQ010")
  # 确认所有节点都存在
  physio_nodes <- physio_nodes[physio_nodes %in% names(data)]
  psycho_nodes <- psycho_nodes[psycho_nodes %in% names(data)]
  all_nodes <- c(physio_nodes, psycho_nodes)
  nodes_present <- all_nodes[all_nodes %in% names(data)]
  cat(sprintf("  可用节点: %d/9 (已移除完全缺失的alt_ui_l)\n", length(nodes_present)))
  cat("  生理节点:", paste(physio_nodes, collapse=", "), "\n")
  cat("  心理节点:", paste(psycho_nodes, collapse=", "), "\n")
  if(length(nodes_present) >= 8) {
    hcf_levels <- na.omit(unique(data$HCF_type))
    # ========================================
    # 生成PDF（高质量，用于发表）
    # ========================================
    cat("\n  生成PDF...\n")
    pdf(pdf_path, width = 15, height = 12)
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
    for(hcf in hcf_levels) {
      # 获取该分型的样本
      hcf_data <- data[data$HCF_type == hcf, nodes_present, drop = FALSE]
      if(nrow(hcf_data) < 10) {
        cat(sprintf("  跳过 %s: 样本量太少 (%d)\n", hcf, nrow(hcf_data)))
        next
      }
      # 计算相关矩阵
      cor_mat <- tryCatch({
        # 使用pairwise.complete.obs处理缺失值
        cor_mat <- cor(hcf_data, use = "pairwise.complete.obs")
        # 将对角线设为0（不显示自相关）
        diag(cor_mat) <- 0
        # 将NA转为0（表示无相关性）
        cor_mat[is.na(cor_mat)] <- 0
        # 检查是否有有效相关性
        if(all(cor_mat == 0)) {
          cat(sprintf("  警告: %s 没有有效相关性\n", hcf))
        }
        cor_mat
      }, error = function(e) {
        cat(sprintf("  相关矩阵计算错误: %s\n", e$message))
        return(NULL)
      })
      if(is.null(cor_mat)) next
      # 节点颜色：生理节点蓝色，心理节点红色
      node_colors <- ifelse(colnames(cor_mat) %in% physio_nodes, 
                           "lightblue", "lightcoral")
      # 获取英文标题
      hcf_title <- hcf_names_en[as.character(hcf)]
      if(is.na(hcf_title)) hcf_title <- as.character(hcf)
      # 绘制网络图
      qgraph(cor_mat,
             layout = "spring",
             title = hcf_title,
             color = node_colors,
             borders = TRUE,
             border.width = 1.5,
             border.color = "black",
             labels = colnames(cor_mat),
             label.cex = 1.2,
             title.cex = 1.5,
             vsize = 8,
             esize = 4,
             edge.color = "grey50",
             edge.width = 1.5,
             posCol = "darkgreen",
             negCol = "darkred",
             maximum = 1,
             cut = 0.1,
             details = TRUE,
             label.scale = FALSE,
             label.prop = 1,
             fade = FALSE)
    }
    dev.off()
    cat(sprintf("  ✅ PDF已保存: %s\n", pdf_path))
    # ========================================
    # 生成PNG（用于快速预览和网页）
    # ========================================
    cat("\n  生成PNG...\n")
    png(png_path, width = 1500, height = 1200, res = 150)
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
    for(hcf in hcf_levels) {
      hcf_data <- data[data$HCF_type == hcf, nodes_present, drop = FALSE]
      if(nrow(hcf_data) < 10) next
      # 计算相关矩阵
      cor_mat <- tryCatch({
        cor_mat <- cor(hcf_data, use = "pairwise.complete.obs")
        diag(cor_mat) <- 0
        cor_mat[is.na(cor_mat)] <- 0
        cor_mat
      }, error = function(e) {
        return(NULL)
      })
      if(is.null(cor_mat)) next
      # 节点颜色
      node_colors <- ifelse(colnames(cor_mat) %in% physio_nodes, 
                           "lightblue", "lightcoral")
      # 获取英文标题
      hcf_title <- hcf_names_en[as.character(hcf)]
      if(is.na(hcf_title)) hcf_title <- as.character(hcf)
      # 绘制网络图
      qgraph(cor_mat,
             layout = "spring",
             title = hcf_title,
             color = node_colors,
             borders = TRUE,
             border.width = 1.5,
             border.color = "black",
             labels = colnames(cor_mat),
             label.cex = 1.2,
             title.cex = 1.5,
             vsize = 8,
             esize = 4,
             edge.color = "grey50",
             edge.width = 1.5,
             posCol = "darkgreen",
             negCol = "darkred",
             maximum = 1,
             cut = 0.1,
             details = TRUE,
             label.scale = FALSE,
             label.prop = 1,
             fade = FALSE)
    }
    dev.off()
    cat(sprintf("  ✅ PNG已保存: %s\n", png_path))
    # ========================================
    # 生成节点列表说明文件
    # ========================================
    node_info <- data.frame(
      Node = nodes_present,
 Type = ifelse(nodes_present %in% physio_nodes, "Physiological", "Psychological"),
      Available_in_P_cycle = "Yes",
      Notes = ifelse(nodes_present == "alt_ui_l", "Removed due to 100% missing", "")
    )
    write.csv(node_info, file.path(RESULTS_DIR, "Figure1_nodes_P.csv"), row.names = FALSE)
    cat("\n  ✅ 节点信息已保存: Figure1_nodes_P.csv\n")
  } else {
    cat("  错误: 可用节点不足，无法生成图1\n")
  }
} else {
  cat("  错误: qgraph包未加载，无法生成图1\n")
}
cat("\n")
# ============================================================================
# 11. 网络指标计算
# ============================================================================
cat("\n9. 计算网络指标...\n")
if(require(qgraph) && require(igraph)) {
  physio_nodes <- c("BMXBMI", "LBXHSCRP", "BPXOPLS1", "alt_ui_l", "egfr_ckdepi_2021")
  psycho_nodes <- c("DPQ010", "DPQ040", "DPQ070", "DPQ090", "HUQ010")
  all_nodes <- c(physio_nodes, psycho_nodes)
  nodes_present <- all_nodes[all_nodes %in% names(data)]
  if(length(nodes_present) >= 8) {
    hcf_levels <- na.omit(unique(data$HCF_type))
    network_metrics <- data.frame()
    for(hcf in hcf_levels) {
      subset_data <- data[data$HCF_type == hcf, nodes_present, drop = FALSE]
      if(nrow(subset_data) < 10) next
      cor_mat <- tryCatch(cor(subset_data, use = "pairwise.complete.obs"),
                          error = function(e) NULL)
      if(is.null(cor_mat)) next
      g <- tryCatch({
        graph_from_adjacency_matrix(abs(cor_mat) > 0.1, mode = "undirected", diag = FALSE)
      }, error = function(e) NULL)
      if(!is.null(g) && vcount(g) > 0) {
        network_metrics <- rbind(network_metrics, data.frame(
          HCF_type = hcf_names_en[as.character(hcf)],
          density = tryCatch(edge_density(g), error = function(e) NA),
          transitivity = tryCatch(transitivity(g), error = function(e) NA),
          n_edges = ecount(g),
          n_nodes = vcount(g),
          stringsAsFactors = FALSE
        ))
      }
    }
    if(nrow(network_metrics) > 0) {
      write.csv(network_metrics, file.path(RESULTS_DIR, "Figure3_data_P.csv"), row.names = FALSE)
      cat(" ✅ 网络指标已保存: Figure3_data_P.csv\n")
    }
  }
}
# ============================================================================
# 12. 生成分析报告
# ============================================================================
cat("\n10. 生成分析报告...\n")
report_file <- file.path(LOG_DIR, "09_paper1_report_P.txt")
sink(report_file)
cat("论文1：HCF宏观分型验证分析报告 (P周期)\n")
cat("========================================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、样本信息\n")
cat(sprintf(" 总样本: %d\n", nrow(data_full)))
cat(sprintf(" 分析人群: %d (%.1f%%)\n", nrow(data), nrow(data)/nrow(data_full)*100))
cat(sprintf(" HCF分型可用: %d (%.1f%%)\n\n",
            sum(!is.na(data$HCF_type)),
            mean(!is.na(data$HCF_type))*100))
cat("二、HCF分型分布（加权）\n")
hcf_prop <- svymean(~HCF_type, mec_design, na.rm = TRUE)
print(round(hcf_prop * 100, 1))
cat("\n三、已保存文件\n")
cat(" 1. paper1_table1_demographics_P.csv\n")
cat(" 2. paper1_table2_mental_outcomes_P.csv\n")
cat(" 3. paper1_table3_physical_outcomes_P.csv\n")
cat(" 4. Figure3_P.pdf\n")
cat(" 5. Figure3_P.png\n")
cat(" 6. Figure3_data_P.csv\n")
cat(" 7. Table4_composite_P.csv\n")
sink()
cat(" ✅ 分析报告已保存\n\n")
# ============================================================================
# 13. 保存会话信息（期刊要求）
# ============================================================================
cat("11. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "09_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期论文1分析会话信息\n")
cat("==============================\n")
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
# 14. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n12. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "09_paper1_HCF_analysis_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "09_code_list_P.txt")
cat("脚本名称: 09_paper1_HCF_analysis_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 15. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期论文1分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果已保存至:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "RESULTS_DIR", "LOG_DIR")))
gc()
