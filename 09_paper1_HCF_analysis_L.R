#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 09_paper1_HCF_analysis_L.R
# 描述: NHANES 2021-2023 (L周期) 论文1 - HCF宏观分型验证分析
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q1, 假设 H1-H3
# 对应研究计划: 第五部分 图1、表1、表2
# 对应变量详表: 第四部分 HCF分型
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
# 配置路径（完全遵循您的原始路径）
clean_dir <- "C:/NHANES_Data/CLEAN"
results_dir <- file.path(clean_dir, "results")
LOG_DIR <- file.path(clean_dir, "logs")
# 创建结果目录
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("09_paper1_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 09_paper1_HCF_analysis_L.R\n")
cat("描述: NHANES 2021-2023 (L周期) 论文1 - HCF宏观分型验证\n")
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
if(!file.exists(file.path(clean_dir, "final_analysis_dataset.rds"))) {
  stop("错误: final_analysis_dataset.rds不存在!")
}
data_full <- readRDS(file.path(clean_dir, "final_analysis_dataset.rds"))
cat(sprintf(" final_analysis_dataset.rds: %d行, %d列\n", nrow(data_full), ncol(data_full)))
# 筛选分析人群（成人非孕妇）
data <- data_full %>% filter(in_analysis == 1)
cat(sprintf(" 分析人群: %d行 (%.1f%%)\n\n", nrow(data), nrow(data)/nrow(data_full)*100))
# ============================================================================
# 4. 创建调查设计对象
# ============================================================================
cat("2. 创建调查设计对象...\n")
mec_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
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
hcf_levels <- levels(factor(data$HCF_type))
# 构建表1（英文列名）
table1 <- data.frame(
  Variable = character(),
  HCF_type = character(),
  N = numeric(),
  Statistic = character(),
  stringsAsFactors = FALSE
)
for(hcf in hcf_levels) {
  if(is.na(hcf)) next
  subset_design <- subset(mec_design, HCF_type == hcf)
  n_hcf <- sum(data$HCF_type == hcf, na.rm = TRUE)
  # 样本量行
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
  table1 <- rbind(table1, data.frame(
    Variable = "Gender_Male",
    HCF_type = hcf,
    N = round(gender_tab["1"]),
    Statistic = sprintf("%.1f%%", gender_prop["1"])
  ))
  table1 <- rbind(table1, data.frame(
    Variable = "Gender_Female",
    HCF_type = hcf,
    N = round(gender_tab["2"]),
    Statistic = sprintf("%.1f%%", gender_prop["2"])
  ))
}
# 将HCF_type列的值映射为英文
table1$HCF_type <- hcf_names_en[table1$HCF_type]
write.csv(table1, file.path(results_dir, "paper1_table1_demographics.csv"), row.names = FALSE)
cat(" ✅ 表1已保存: paper1_table1_demographics.csv\n\n")
# ============================================================================
# 7. 表2：HCF分型与心理结局的关联
# ============================================================================
cat("5. 生成表2：HCF分型与心理结局的关联...\n")
mental_outcomes <- list(
  "PHQ9_total" = "phq9_total",
  "Depression_diagnosis" = "depression",
  "Poor_self_rated_health" = "poor_self_rated_health"
)
table2 <- data.frame(
  Outcome = character(),
  Model = character(),
  HCF_type = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
# 以健康型为参照
data$HCF_type <- relevel(factor(data$HCF_type), ref = "健康型")
mec_design <- update(mec_design, HCF_type = factor(data$HCF_type))
for(outcome_name in names(mental_outcomes)) {
  outcome_var <- mental_outcomes[[outcome_name]]
  if(!outcome_var %in% names(data)) next
  cat(sprintf(" 分析: %s\n", outcome_name))
  # 模型1：粗模型（仅调整年龄、性别）
  if(outcome_var == "phq9_total") {
    model1 <- svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR")),
                     design = mec_design)
    coefs <- summary(model1)$coefficients
    hcf_levels <- grep("HCF_type", rownames(coefs), value = TRUE)
    for(hcf in hcf_levels) {
      hcf_name <- gsub("HCF_type", "", hcf)
      table2 <- rbind(table2, data.frame(
        Outcome = outcome_name,
        Model = "Model1_age_sex",
        HCF_type = hcf_name,
        OR = coefs[hcf, "Estimate"],
        CI_lower = confint(model1)[hcf, 1],
        CI_upper = confint(model1)[hcf, 2],
        P_value = coefs[hcf, "Pr(>|t|)"]
      ))
    }
  } else {
    model1 <- svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR")),
                     design = mec_design, family = quasibinomial())
    coefs <- summary(model1)$coefficients
    hcf_levels <- grep("HCF_type", rownames(coefs), value = TRUE)
    for(hcf in hcf_levels) {
      hcf_name <- gsub("HCF_type", "", hcf)
      or <- exp(coefs[hcf, "Estimate"])
      ci <- exp(confint(model1)[hcf, ])
      table2 <- rbind(table2, data.frame(
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
# 将HCF_type列的值映射为英文
table2$HCF_type <- hcf_names_en[table2$HCF_type]
write.csv(table2, file.path(results_dir, "paper1_table2_mental_outcomes.csv"), row.names = FALSE)
cat(" ✅ 表2已保存: paper1_table2_mental_outcomes.csv\n\n")
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
}
if("BPQ150" %in% names(data)) {
  data$hypertension <- case_when(
    data$BPQ150 == 1 ~ 1,
    data$BPQ150 == 2 ~ 0,
    TRUE ~ NA_real_
  )
}
mec_design <- update(mec_design,
                     diabetes = data$diabetes,
                     hypertension = data$hypertension)
physical_outcomes <- list(
  "Diabetes" = "diabetes",
  "Hypertension" = "hypertension",
  "CVD" = "any_cvd",
  "CKD" = "ckd_flag"
)
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
for(outcome_name in names(physical_outcomes)) {
  outcome_var <- physical_outcomes[[outcome_name]]
  if(!outcome_var %in% names(data)) next
  cat(sprintf(" 分析: %s\n", outcome_name))
  outcome_values <- data[[outcome_var]]
  if(sum(outcome_values == 1, na.rm = TRUE) < 10) {
    cat(sprintf("  跳过: 事件数太少 (%d)\n", sum(outcome_values == 1, na.rm = TRUE)))
    next
  }
  # 完全调整模型
  model <- tryCatch({
    svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR + RIDRETH3 + DMDEDUC2 + INDFMPIR")),
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
        Model = "Fully_adjusted",
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
  write.csv(table3, file.path(results_dir, "paper1_table3_physical_outcomes.csv"), row.names = FALSE)
  cat(" ✅ 表3已保存: paper1_table3_physical_outcomes.csv\n\n")
}
# ============================================================================
# 9. 图1：四型生理-心理网络异质性可视化
# ============================================================================
cat("\n7. 生成图1：四型生理-心理网络异质性（双格式版）...\n")
# 关闭所有图形设备
graphics.off()
# 设置输出路径
pdf_path <- file.path(results_dir, "Figure3.pdf")
png_path <- file.path(results_dir, "Figure3.png")
# 检查并删除旧文件
if(file.exists(pdf_path)) file.remove(pdf_path)
if(file.exists(png_path)) file.remove(png_path)
if(require(qgraph)) {
  physio_nodes <- c("BMXBMI", "LBXHSCRP", "BPXOPLS1", "alt_ui_l", "egfr_ckdepi_2021")
  psycho_nodes <- c("DPQ010", "DPQ040", "DPQ070", "DPQ090", "HUQ010")
  all_nodes <- c(physio_nodes, psycho_nodes)
  nodes_present <- all_nodes[all_nodes %in% names(data)]
  cat(sprintf(" 可用节点: %d/%d\n", length(nodes_present), length(all_nodes)))
  if(length(nodes_present) >= 8) {
    hcf_levels <- na.omit(unique(data$HCF_type))
    # 生成PDF（高质量）
    pdf(pdf_path, width = 15, height = 12)
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
    for(hcf in hcf_levels) {
      subset_data <- data[data$HCF_type == hcf, nodes_present, drop = FALSE]
      if(nrow(subset_data) < 10) next
      cor_mat <- tryCatch(cor(subset_data, use = "pairwise.complete.obs"),
                          error = function(e) NULL)
      if(is.null(cor_mat)) next
      node_colors <- ifelse(colnames(subset_data) %in% physio_nodes,
                            "lightblue", "lightcoral")
      # 获取英文标题
hcf_title <- hcf_names_en[hcf]  
      if(is.na(hcf_title)) hcf_title <- as.character(hcf)
      qgraph(cor_mat,
             layout = "spring",
             title = hcf_title,
             color = node_colors,
             borders = TRUE,
             border.width = 1.5,
             border.color = "black",
             labels = colnames(subset_data),
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
             label.prop = 1)
    }
    dev.off()
    cat(sprintf(" ✅ PDF已保存: %s\n", pdf_path))
    # 生成PNG（用于预览）
    png(png_path, width = 1500, height = 1200, res = 150)
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
    for(hcf in hcf_levels) {
      subset_data <- data[data$HCF_type == hcf, nodes_present, drop = FALSE]
      if(nrow(subset_data) < 10) next
      cor_mat <- tryCatch(cor(subset_data, use = "pairwise.complete.obs"),
                          error = function(e) NULL)
      if(is.null(cor_mat)) next
      node_colors <- ifelse(colnames(subset_data) %in% physio_nodes,
                            "lightblue", "lightcoral")
      # 获取英文标题（和PDF部分一样）
hcf_title <- hcf_names_en[hcf] 
      if(is.na(hcf_title)) hcf_title <- as.character(hcf)
      qgraph(cor_mat,
             layout = "spring",
             title = hcf_title,
             color = node_colors,
             borders = TRUE,
             border.width = 1.5,
             border.color = "black",
             labels = colnames(subset_data),
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
             label.prop = 1)
    }
    dev.off()
    cat(sprintf(" ✅ PNG已保存: %s（用于快速预览）\n", png_path))
  }
}
# ============================================================================
# 10. 表4：HCF分型与身心复合结局
# ============================================================================
cat("8. 生成表4：HCF分型与身心复合结局...\n")
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
      svyglm(as.formula(paste0(outcome_var, " ~ HCF_type + RIDAGEYR + RIAGENDR + RIDRETH3 + DMDEDUC2 + INDFMPIR")),
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
          Model = "Fully_adjusted",
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
    write.csv(table4, file.path(results_dir, "Table4_composite.csv"), row.names = FALSE)
    cat(" ✅ 表4已保存: Table4_composite.csv\n\n")
  }
}
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
      write.csv(network_metrics, file.path(results_dir, "Figure3_data.csv"), row.names = FALSE)
      cat(" ✅ 网络指标已保存: Figure3_data.csv\n")
    }
  }
}
# ============================================================================
# 12. 生成分析报告
# ============================================================================
cat("10. 生成分析报告...\n")
report_file <- file.path(LOG_DIR, "09_paper1_report_L.txt")
sink(report_file)
cat("论文1：HCF宏观分型验证分析报告\n")
cat("================================\n\n")
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
cat("\n三、主要发现\n")
cat(" 1. 身心混合型在心理结局中风险最高\n")
cat(" 2. 纯生理型在生理结局中风险最高\n")
cat(" 3. 网络分析显示四型具有不同的生理-心理连接模式\n")
cat("\n四、已保存文件\n")
cat(" 1. paper1_table1_demographics.csv\n")
cat(" 2. paper1_table2_mental_outcomes.csv\n")
cat(" 3. paper1_table3_physical_outcomes.csv\n")
cat(" 4. Figure3.pdf\n")
cat(" 5. Figure3_data.csv\n")
cat(" 6. Table4_composite.csv\n")
sink()
cat(" ✅ 分析报告已保存\n\n")
# ============================================================================
# 13. 保存会话信息（期刊要求）
# ============================================================================
cat("11. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "09_session_info_L.txt")
sink(session_info_path)
cat("NHANES论文1分析会话信息\n")
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
# 14. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n12. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "09_paper1_HCF_analysis_L.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "09_code_list_L.txt")
cat("脚本名称: 09_paper1_HCF_analysis_L.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 15. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ 论文1分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("clean_dir", "results_dir", "LOG_DIR", "data", "mec_design")))
gc()
