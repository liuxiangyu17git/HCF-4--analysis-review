#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 08_final_merge_P.R
# 描述: NHANES 2017-2020 (P周期) 最终数据合并 - 阶段六
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: OSF https://osf.io/xxxxx (DOI: 10.17605/OSF.IO/XXXXX)
# 对应研究计划: 第四部分 4.2 变量构建总览
# 对应变量详表: 第三部分 构建阶段详表（阶段六 - 最终合并）
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
required_packages <- c("tidyverse", "survey")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# 配置路径 - P周期独立目录
PROJECT_ROOT <- "C:/NHANES_Data"
CLEAN_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")  # P周期清洗后文件目录
LOG_DIR <- file.path(CLEAN_DATA_DIR, "logs")
# 创建日志目录
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录（期刊要求）
log_file <- file.path(LOG_DIR, paste0("08_final_merge_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 08_final_merge_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 最终数据合并\n")
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
# 2. 加载master_P.rds
# ============================================================================
cat("1. 加载master_P.rds...\n")
master_file <- file.path(CLEAN_DATA_DIR, "master_P.rds")
if(!file.exists(master_file)) {
  stop("错误: master_P.rds不存在! 请先运行02_data_merge_P.R")
}
master_base <- readRDS(master_file)
cat(sprintf(" ✅ master_P.rds: %d行, %d列\n\n", nrow(master_base), ncol(master_base)))
# ============================================================================
# 3. 加载其他构建文件
# ============================================================================
cat("2. 加载其他构建文件...\n")
files_to_load <- list(
  behavior = "behavior_vars_P.rds",
  clinical = "clinical_vars_P.rds",
  social = "social_vars_P.rds",
  symptom = "symptom_vars_P.rds",
  dietbehavior = "dietbehavior_vars_P.rds",
  healthbehavior = "healthbehavior_vars_P.rds",
  biomarker_core = "02_biomarker_core_P.rds",
  biomarker_full = "02_biomarker_full_P.rds",
  alpha = "alpha_factors_P.rds",
  hcf = "HCF_typing_P.rds",
  pathway = "pathway_final_with_clusters_P.rds"
)
data_list <- list()
for(name in names(files_to_load)) {
  file_path <- files_to_load[[name]]
  full_path <- file.path(CLEAN_DATA_DIR, file_path)
  if(file.exists(full_path)) {
    data_list[[name]] <- readRDS(full_path)
    cat(sprintf(" ✅ %s: %s (%d行, %d列)\n",
                name, file_path,
                nrow(data_list[[name]]),
                ncol(data_list[[name]])))
  } else {
    cat(sprintf(" ⚠️ %s: %s 不存在\n", name, file_path))
  }
}
cat("\n")
# ============================================================================
# 4. 合并数据
# ============================================================================
cat("3. 开始合并数据...\n")
final_data <- master_base
cat(sprintf(" 基础: master_P.rds (%d行, %d列)\n", nrow(final_data), ncol(final_data)))
for(name in names(data_list)) {
  current_data <- data_list[[name]]
  if("SEQN" %in% names(current_data)) {
    new_vars <- setdiff(names(current_data), names(final_data))
    if(length(new_vars) > 0) {
      current_data_subset <- current_data %>%
        select(SEQN, all_of(new_vars))
      final_data <- final_data %>%
        left_join(current_data_subset, by = "SEQN")
      cat(sprintf(" ✅ 合并 %s: 新增 %d 个变量\n", name, length(new_vars)))
    }
  }
}
cat(sprintf("\n 最终数据集: %d行, %d列\n\n", nrow(final_data), ncol(final_data)))
# ============================================================================
# 5. 处理NHANES缺失值编码（必须保留！）
# ============================================================================
cat("4. 处理NHANES缺失值编码...\n")
nhanes_na_codes <- c(7, 9, 77, 99, 777, 999, 7777, 9999)
replace_nhanes_na <- function(x) {
  if(is.numeric(x)) {
    x[x %in% nhanes_na_codes] <- NA
  }
  return(x)
}
numeric_cols <- sapply(final_data, is.numeric)
final_data[, numeric_cols] <- lapply(final_data[, numeric_cols], replace_nhanes_na)
cat(" ✅ NHANES缺失值编码处理完成\n\n")
# ============================================================================
# 6. 添加分析变量
# ============================================================================
cat("5. 添加分析变量...\n")
final_data <- final_data %>%
  mutate(
    # 药物变量 - P周期使用BPQ020（降压药）
    antihypertensive_med = if("BPQ020" %in% names(final_data)) {
      case_when(
        BPQ020 == 1 ~ 1,
        BPQ020 == 2 ~ 0,
        TRUE ~ NA_real_
      )
    } else {
      NA_real_
    },
    antidiabetic_med = if("DIQ070" %in% names(final_data)) {
      case_when(
        DIQ070 == 1 ~ 1,
        DIQ070 == 2 ~ 0,
        TRUE ~ NA_real_
      )
    } else {
      NA_real_
    },
    # 权重选择 - P周期使用PRP权重
    analysis_weight = case_when(
      !is.na(LBXHSCRP) | !is.na(BMXBMI) ~ WTMECPRP,
      TRUE ~ WTINTPRP
    ),
    weight_source = case_when(
      !is.na(LBXHSCRP) | !is.na(BMXBMI) ~ "MEC",
      TRUE ~ "INTERVIEW"
    ),
    # 分析人群标记
    in_analysis = case_when(
      RIDAGEYR >= 18 & (RIDEXPRG != 1 | is.na(RIDEXPRG)) ~ 1,
      TRUE ~ 0
    ),
    # 年龄分组
    age_group_std = case_when(
      RIDAGEYR >= 18 & RIDAGEYR < 40 ~ "18-39",
      RIDAGEYR >= 40 & RIDAGEYR < 60 ~ "40-59",
      RIDAGEYR >= 60 ~ "60+",
      TRUE ~ NA_character_
    )
  )
cat(" ✅ 添加完成\n\n")
# ============================================================================
# 7. 验证关键变量
# ============================================================================
cat("6. 验证关键变量...\n")
core_vars <- list(
  paper1 = c("BMXBMI", "BPXOSY1", "BPXODI1", "LBXHSCRP", "LBXWBCSI",
             "antihypertensive_med", "antidiabetic_med",
             "HCF_A", "HCF_B", "HCF_C", "HCF_D", "HCF_type"),
  paper2 = c("avoidance_z", "perseveration_z", "hyperactivation_z", "exhaustion_z",
             "pathway_cluster", "pathway_cluster_name"),
  paper3 = c("alpha1", "alpha2", "alpha3", "alpha4",
             "HUQ010", "DPQ010", "DPQ040", "DPQ090"),
  common = c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3",
             "SDMVSTRA", "SDMVPSU", "WTMECPRP", "WTINTPRP", "WTDRD1")
)
for(paper in names(core_vars)) {
  vars_needed <- core_vars[[paper]]
  vars_present <- vars_needed[vars_needed %in% names(final_data)]
  cat(sprintf(" %s: %d/%d (%.1f%%)\n",
              paper, length(vars_present), length(vars_needed),
              length(vars_present)/length(vars_needed)*100))
}
cat("\n")
# ============================================================================
# 8. 保存文件
# ============================================================================
cat("7. 保存文件...\n")
# 保存完整数据集
saveRDS(final_data, file.path(CLEAN_DATA_DIR, "final_analysis_dataset_P.rds"))
cat(" ✅ final_analysis_dataset_P.rds\n")
# 保存分析子集
analysis_data <- final_data %>% filter(in_analysis == 1)
saveRDS(analysis_data, file.path(CLEAN_DATA_DIR, "analysis_dataset_subset_P.rds"))
cat(" ✅ analysis_dataset_subset_P.rds\n")
# ============================================================================
# 9. 数据字典
# ============================================================================
cat("\n8. 创建数据字典...\n")
data_dictionary <- data.frame(
  variable = names(final_data),
  class = sapply(final_data, function(x) class(x)[1]),
  missing_n = sapply(final_data, function(x) sum(is.na(x))),
  missing_pct = round(sapply(final_data, function(x) mean(is.na(x)) * 100), 1)
)
saveRDS(data_dictionary, file.path(CLEAN_DATA_DIR, "data_dictionary_P.rds"))
write.csv(data_dictionary, file.path(CLEAN_DATA_DIR, "data_dictionary_P.csv"), row.names = FALSE)
cat(" ✅ data_dictionary_P.rds 和 data_dictionary_P.csv\n")
# ============================================================================
# 10. 生成合并报告
# ============================================================================
cat("\n9. 生成合并报告...\n")
report_file <- file.path(LOG_DIR, "08_final_merge_report_P.txt")
sink(report_file)
cat("P周期最终数据合并报告\n")
cat("======================\n\n")
cat("合并时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("最终数据集: ", nrow(final_data), "行, ", ncol(final_data), "列\n")
cat("分析人群: ", sum(final_data$in_analysis == 1, na.rm = TRUE), "\n\n")
cat("变量分类:\n")
cat(" BPX变量:", length(grep("BPX", names(final_data))), "\n")
cat(" LBX变量:", length(grep("LBX", names(final_data))), "\n")
cat(" 问卷变量:", length(grep("BPQ|DIQ|HUQ|DPQ", names(final_data))), "\n")
cat(" 构建变量:", length(grep("alpha|HCF|avoid|persev|hyper|exhaust", names(final_data))), "\n")
sink()
cat(" ✅ 合并报告已保存\n\n")
# ============================================================================
# 11. 保存会话信息（期刊要求）
# ============================================================================
cat("10. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "08_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期最终合并会话信息\n")
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
# 12. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n11. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "08_final_merge_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "08_code_list_P.txt")
cat("脚本名称: 08_final_merge_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 13. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期最终数据合并完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", CLEAN_DATA_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "LOG_DIR", "final_data", "analysis_data")))
gc()
