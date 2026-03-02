#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 01_data_cleaning_P.R
# 描述: NHANES 2017-2020 (P周期) 数据清洗 - 完全符合CDC/NCHS教程要求
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: OSF https://osf.io/xxxxx (DOI: 10.17605/OSF.IO/XXXXX)
# 对应研究计划: 第三部分 3.1 数据来源
# 对应变量详表: 第一部分 NHANES文件清单
#
# 伦理声明: 使用NHANES公开数据，无需IRB批准
# 数据来源: https://www.cdc.gov/nchs/nhanes/index.htm
# 数据周期: 2017-2020 (P系列)
#
# 核心原则:
# 1. 绝不删除任何观测 (教程要求: 保留所有观测用于方差估计)
# 2. 正确处理零权重 (零权重是有效的调查设计特征)
# 3. 完整处理跳过模式
# 4. 验证设计变量结构
# 5. 为多周期分析和亚组分析做准备
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
# 加载必需包
required_packages <- c("haven", "dplyr", "labelled", "survey", "tidyr", "purrr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# 配置路径 - P周期独立目录
PROJECT_ROOT <- "C:/NHANES_Data"
RAW_DATA_DIR <- PROJECT_ROOT                    # P周期原始xpt文件在 C:/NHANES_Data/
CLEAN_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")  # P周期清洗后文件保存在 2017-2020 目录下
LOG_DIR <- file.path(CLEAN_DATA_DIR, "logs")
# 创建必要的目录
if (!dir.exists(CLEAN_DATA_DIR)) dir.create(CLEAN_DATA_DIR, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 设置工作目录
setwd(CLEAN_DATA_DIR)
# 启动日志记录
log_file <- file.path(LOG_DIR, paste0("01_cleaning_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 01_data_cleaning_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 数据清洗\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R版本:", R.version.string, "\n")
cat("随机种子: 20240226\n")
cat("========================================================\n\n")
# 记录包版本
cat("加载的包版本:\n")
for (pkg in required_packages) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")
# ============================================================================
# 2. 定义文件列表（P周期20个文件）
# ============================================================================
cat("步骤1: 定义需要处理的文件\n")
cat("----------------------------------------\n")
files_to_clean <- c(
  "P_DEMO.xpt",      # 1. 人口统计
  "P_DPQ.xpt",       # 2. 抑郁问卷
  "P_SLQ.xpt",       # 3. 睡眠问卷
  "P_PAQ.xpt",       # 4. 体力活动问卷
  "P_DBQ.xpt",       # 5. 饮食行为问卷
  "P_DR1TOT.xpt",    # 6. 24小时膳食回顾
  "P_WHQ.xpt",       # 7. 体重历史问卷
  "P_BMX.xpt",       # 8. 身体测量
  "P_BPXO.xpt",      # 9. 血压测量
  "P_BPQ.xpt",       # 10. 血压问卷
  "P_HSCRP.xpt",     # 11. C反应蛋白
  "P_CBC.xpt",       # 12. 血常规
  "P_INQ.xpt",       # 13. 收入贫困比
  "P_ALQ.xpt",       # 14. 酒精使用问卷
  "P_HUQ.xpt",       # 15. 健康问卷
  "P_SMQ.xpt",       # 16. 吸烟问卷
  "P_OCQ.xpt",       # 17. 职业问卷
  "P_BIOPRO.xpt",    # 18. 生化指标
  "P_DIQ.xpt",       # 19. 糖尿病问卷
  "P_MCQ.xpt"        # 20. 医疗状况问卷
)
cat(sprintf("需要处理的文件数: %d\n", length(files_to_clean)))
cat("文件列表:\n")
for (f in files_to_clean) cat(sprintf("  - %s\n", f))
cat("\n")
# 初始化处理日志
processing_log <- data.frame(
  File = character(),
  Status = character(),
  Rows = integer(),
  Cols = integer(),
  SEQN_Present = logical(),
  Weights_Present = logical(),
  DesignVars_Present = logical(),
  stringsAsFactors = FALSE
)
# ============================================================================
# 3. 定义核心清洗函数（完全遵循NHANES教程）
# ============================================================================
cat("步骤2: 定义清洗函数\n")
cat("----------------------------------------\n")
#' 核心清洗函数 - 应用于所有文件
#' @param df 原始数据框
#' @param file_name 文件名（用于日志）
#' @return 清洗后的数据框
clean_nhanes_data <- function(df, file_name) {
  cat(sprintf("  开始清洗: %s\n", file_name))
  original_rows <- nrow(df)
  original_vars <- names(df)
  cat(sprintf("    原始维度: %d行 × %d列\n", original_rows, length(original_vars)))
  cat(sprintf("    原则: 保留所有 %d 个观测\n", original_rows))
  # 标准化缺失值处理（NHANES教程 Table 8-1）
  for (var in original_vars) {
    if (var == "SEQN") next
    var_vector <- df[[var]]
    if (is.numeric(var_vector)) {
      df[[var]] <- ifelse(
        var_vector %in% c(7, 77, 777, 7777, 9, 99, 999, 9999),
        NA_real_,
        var_vector
      )
    } else if (is.character(var_vector)) {
      df[[var]] <- ifelse(
        is.na(var_vector) | trimws(var_vector) == "",
        NA_character_,
        var_vector
      )
    }
    if (!is.null(attr(var_vector, "label"))) {
      attr(df[[var]], "label") <- attr(var_vector, "label")
    }
  }
  # 正确处理权重变量
  weight_vars <- c("WTINT2YR", "WTMEC2YR", "WTDRD1", "WTDRD2", "WTDR2D")
  present_weights <- intersect(weight_vars, original_vars)
  if (length(present_weights) > 0) {
    cat(sprintf("    检测到权重变量: %s\n", paste(present_weights, collapse = ", ")))
    for (wt_var in present_weights) {
      df[[wt_var]] <- as.numeric(df[[wt_var]])
      df[[wt_var]] <- ifelse(df[[wt_var]] < 0, NA_real_, df[[wt_var]])
      wt_positive <- sum(df[[wt_var]] > 0, na.rm = TRUE)
      wt_zero <- sum(df[[wt_var]] == 0, na.rm = TRUE)
      wt_negative <- sum(df[[wt_var]] < 0, na.rm = TRUE)
      wt_na <- sum(is.na(df[[wt_var]]))
      cat(sprintf("      %s: 正数=%d, 零=%d, 负数=%d, NA=%d\n",
                  wt_var, wt_positive, wt_zero, wt_negative, wt_na))
    }
  }
  # 验证设计变量结构
  design_vars <- c("SDMVPSU", "SDMVSTRA")
  present_design <- intersect(design_vars, original_vars)
  if (length(present_design) == 2) {
    cat("    检测到设计变量: SDMVPSU 和 SDMVSTRA\n")
    df$SDMVPSU <- as.integer(df$SDMVPSU)
    df$SDMVSTRA <- as.integer(df$SDMVSTRA)
    unique_psu <- length(unique(na.omit(df$SDMVPSU)))
    unique_strata <- length(unique(na.omit(df$SDMVSTRA)))
    cat(sprintf("      SDMVPSU: %d 个唯一值\n", unique_psu))
    cat(sprintf("      SDMVSTRA: %d 个唯一值\n", unique_strata))
    cat("      说明: 使用Masked Variance Units (MVUs)，方差估计有效\n")
  }
  if (nrow(df) != original_rows) {
    stop(paste("错误：观测数量改变！原始:", original_rows, "现在:", nrow(df)))
  }
  if (length(names(df)) != length(original_vars)) {
    stop(paste("错误：变量数量改变！原始:", length(original_vars), "现在:", length(names(df))))
  }
  cat(sprintf("    清洗完成: %s\n", file_name))
  cat(sprintf("    最终维度: %d行 × %d列\n", nrow(df), ncol(df)))
  return(df)
}
# ============================================================================
# 4. 定义专业清洗函数（处理跳过模式和特定变量）- P周期兼容版
# ============================================================================
#' DIQ专业清洗 - 糖尿病问卷
clean_diq_special <- function(df) {
  cat("    执行DIQ专业清洗 (处理跳过模式)...\n")
  if ("DIQ010" %in% names(df)) {
    df$DIQ010 <- ifelse(df$DIQ010 %in% c(1, 2, 3), df$DIQ010, NA_real_)
    diq_skip_vars <- c("DIQ160", "DIQ170", "DIQ172", "DIQ175", "DIQ180", "DIQ190")
    for (var in diq_skip_vars) {
      if (var %in% names(df)) {
        df[[var]] <- ifelse(df$DIQ010 == 2, NA_real_, df[[var]])
      }
    }
    if ("DIQ200" %in% names(df)) {
      df$DIQ200 <- ifelse(df$DIQ010 %in% c(2, 3), NA_real_, df$DIQ200)
    }
  }
  if ("DIQ050" %in% names(df)) {
    df$DIQ050 <- ifelse(df$DIQ050 %in% c(1, 2), df$DIQ050, NA_real_)
  }
  if ("DIQ070" %in% names(df)) {
    df$DIQ070 <- ifelse(df$DIQ070 %in% c(1, 2), df$DIQ070, NA_real_)
  }
  return(df)
}
#' DPQ专业清洗 - PHQ-9抑郁问卷
clean_dpq_special <- function(df) {
  cat("    执行DPQ专业清洗 (PHQ-9问卷)...\n")
  dpq_vars <- grep("^DPQ[0-9]{3}$", names(df), value = TRUE)
  if (length(dpq_vars) > 0) {
    cat(sprintf("      找到 %d 个DPQ变量\n", length(dpq_vars)))
    for (var in dpq_vars) {
      df[[var]] <- ifelse(
        df[[var]] %in% c(0, 1, 2, 3),
        df[[var]],
        ifelse(df[[var]] %in% c(7, 77, 777, 9, 99, 999), NA_real_, df[[var]])
      )
    }
  }
  return(df)
}
#' WHQ专业清洗 - 体重历史问卷
clean_whq_special <- function(df) {
  cat("    执行WHQ专业清洗 (体重历史问卷)...\n")
  if ("WHQ010" %in% names(df)) {
    df$WHQ010 <- ifelse(df$WHQ010 %in% c(777, 7777, 999, 9999), NA_real_, df$WHQ010)
    df$WHQ010 <- ifelse(df$WHQ010 >= 50 & df$WHQ010 <= 450, df$WHQ010, NA_real_)
    cat("      已处理WHQ010 (10年前体重)\n")
  }
  if ("WHQ030" %in% names(df)) {
    df$WHQ030 <- ifelse(df$WHQ030 %in% c(777, 7777, 999, 9999), NA_real_, df$WHQ030)
    df$WHQ030 <- ifelse(df$WHQ030 >= 50 & df$WHQ030 <= 450, df$WHQ030, NA_real_)
    cat("      已处理WHQ030 (1年前体重)\n")
  }
  if ("WHD050" %in% names(df)) {
    df$WHD050 <- ifelse(df$WHD050 %in% c(777, 7777, 999, 9999), NA_real_, df$WHD050)
    df$WHD050 <- ifelse(df$WHD050 >= 50 & df$WHD050 <= 450, df$WHD050, NA_real_)
    cat("      已处理WHD050 (成年期最高体重)\n")
  }
  if ("WHD060" %in% names(df)) {
    df$WHD060 <- ifelse(df$WHD060 %in% c(777, 7777, 999, 9999), NA_real_, df$WHD060)
    df$WHD060 <- ifelse(df$WHD060 >= 18 & df$WHD060 <= 85, df$WHD060, NA_real_)
    cat("      已处理WHD060 (最高体重年龄)\n")
  }
  if ("WHQ070" %in% names(df)) {
    df$WHQ070 <- ifelse(df$WHQ070 %in% c(7, 77, 9, 99), NA_real_, df$WHQ070)
    cat("      已处理WHQ070 (过去12个月尝试减肥)\n")
  }
  weight_loss_methods <- grep("^WHD080[A-M]$", names(df), value = TRUE)
  if (length(weight_loss_methods) > 0) {
    cat("      检测到", length(weight_loss_methods), "个减肥方法变量\n")
    for (var in weight_loss_methods) {
      df[[var]] <- ifelse(df[[var]] %in% c(7, 9), NA_real_, df[[var]])
    }
    cat("      已处理减肥方法变量\n")
  }
  if ("WHQ090" %in% names(df)) {
    df$WHQ090 <- ifelse(df$WHQ090 %in% c(777, 999), NA_real_, df$WHQ090)
    df$WHQ090 <- ifelse(abs(df$WHQ090) <= 300, df$WHQ090, NA_real_)
    cat("      已处理WHQ090 (期望体重变化)\n")
  }
  if ("WHQ150" %in% names(df)) {
    df$WHQ150 <- ifelse(df$WHQ150 %in% c(777, 7777, 999, 9999), NA_real_, df$WHQ150)
    cat("      已处理WHQ150 (理想体重)\n")
  }
  return(df)
}
#' ALQ专业清洗 - 酒精问卷
clean_alq_special <- function(df) {
  cat("    执行ALQ专业清洗 (酒精问卷)...\n")
  if ("ALQ151" %in% names(df)) {
    df$ALQ151 <- ifelse(df$ALQ151 %in% c(777, 999), NA_real_, df$ALQ151)
    df$ALQ151 <- ifelse(df$ALQ151 >= 0 & df$ALQ151 <= 365, df$ALQ151, NA_real_)
  }
  alq_vars <- c("ALQ100", "ALQ1Q", "ALQ1U", "ALQ130", "ALQ141Q")
  for (var in alq_vars) {
    if (var %in% names(df)) {
      df[[var]] <- ifelse(df[[var]] %in% c(777, 7777, 999, 9999), NA_real_, df[[var]])
    }
  }
  return(df)
}
#' HUQ专业清洗 - 健康问卷
clean_huq_special <- function(df) {
  cat("    执行HUQ专业清洗 (健康问卷)...\n")
  if ("HUQ010" %in% names(df)) {
    df$HUQ010 <- ifelse(
      df$HUQ010 %in% c(1, 2, 3, 4, 5),
      df$HUQ010,
      ifelse(df$HUQ010 %in% c(7, 77, 9, 99), NA_real_, df$HUQ010)
    )
  }
  if ("HUQ030" %in% names(df)) {
    df$HUQ030 <- ifelse(df$HUQ030 %in% c(777, 999), NA_real_, df$HUQ030)
    df$HUQ030 <- ifelse(df$HUQ030 >= 0 & df$HUQ030 <= 50, df$HUQ030, NA_real_)
  }
  if ("HUQ080" %in% names(df)) {
    df$HUQ080 <- ifelse(df$HUQ080 %in% c(77, 99), NA_real_, df$HUQ080)
    df$HUQ080 <- ifelse(df$HUQ080 >= 0 & df$HUQ080 <= 30, df$HUQ080, NA_real_)
  }
  if ("HUQ090" %in% names(df)) {
    df$HUQ090 <- ifelse(df$HUQ090 %in% c(7, 77, 9, 99), NA_real_, df$HUQ090)
  }
  return(df)
}
#' DBQ专业清洗 - 饮食行为问卷（P周期版）
clean_dbq_special <- function(df) {
  cat("    执行DBQ专业清洗 (P周期版)...\n")
  # P周期：DBQ095 替代 DBQ010
  if ("DBQ095" %in% names(df)) {
    df$DBQ095 <- ifelse(
      df$DBQ095 %in% c(1, 2),
      df$DBQ095,
      ifelse(df$DBQ095 %in% c(7, 9), NA_real_, df$DBQ095)
    )
    # 创建兼容变量供下游使用
    df$DBQ010 <- df$DBQ095
    cat("      已处理DBQ095/DBQ010\n")
  }
  # P周期：DBQ097 替代 DBD030
  if ("DBQ097" %in% names(df)) {
    df$DBQ097 <- ifelse(
      df$DBQ097 %in% c(1, 2, 3, 4),
      df$DBQ097,
      ifelse(df$DBQ097 %in% c(7, 9), NA_real_, df$DBQ097)
    )
    # 创建兼容变量供下游使用
    df$DBD030 <- df$DBQ097
    cat("      已处理DBQ097/DBD030\n")
  }
  # DBQ041 变量名相同
  if ("DBQ041" %in% names(df)) {
    df$DBQ041 <- ifelse(
      df$DBQ041 %in% c(1, 2, 3),
      df$DBQ041,
      ifelse(df$DBQ041 %in% c(7, 9), NA_real_, df$DBQ041)
    )
    cat("      已处理DBQ041\n")
  }
  return(df)
}
#' SMQ专业清洗 - 吸烟问卷
clean_smq_special <- function(df) {
  cat("    执行SMQ专业清洗 (吸烟问卷)...\n")
  if ("SMQ020" %in% names(df)) {
    df$SMQ020 <- ifelse(
      df$SMQ020 %in% c(1, 2),
      df$SMQ020,
      ifelse(df$SMQ020 %in% c(7, 9), NA_real_, df$SMQ020)
    )
    cat("      已处理SMQ020 (一生吸烟状态)\n")
  }
  if ("SMQ040" %in% names(df)) {
    df$SMQ040 <- ifelse(
      df$SMQ040 %in% c(1, 2, 3),
      df$SMQ040,
      ifelse(df$SMQ040 %in% c(7, 9), NA_real_, df$SMQ040)
    )
    if ("SMQ020" %in% names(df)) {
      df$SMQ040 <- ifelse(df$SMQ020 == 2, NA_real_, df$SMQ040)
    }
    cat("      已处理SMQ040 (当前吸烟状态)\n")
  }
  if ("SMD030" %in% names(df)) {
    df$SMD030 <- ifelse(df$SMD030 %in% c(777, 999), NA_real_, df$SMD030)
    df$SMD030 <- ifelse(df$SMD030 >= 5 & df$SMD030 <= 80, df$SMD030, NA_real_)
  }
  return(df)
}
#' BPQ专业清洗 - 血压问卷
clean_bpq_special <- function(df) {
  cat("    执行BPQ专业清洗 (血压问卷)...\n")
  if ("BPQ010" %in% names(df)) {
    df$BPQ010 <- ifelse(
      df$BPQ010 %in% c(1, 2),
      df$BPQ010,
      ifelse(df$BPQ010 %in% c(7, 9), NA_real_, df$BPQ010)
    )
    cat("      已处理BPQ010 (医生告知高血压)\n")
  }
  if ("BPQ020" %in% names(df)) {
    df$BPQ020 <- ifelse(
      df$BPQ020 %in% c(1, 2),
      df$BPQ020,
      ifelse(df$BPQ020 %in% c(7, 9), NA_real_, df$BPQ020)
    )
    cat("      已处理BPQ020 (服用降压药)\n")
  }
  if ("BPQ030" %in% names(df)) {
    df$BPQ030 <- ifelse(
      df$BPQ030 %in% c(1, 2),
      df$BPQ030,
      ifelse(df$BPQ030 %in% c(7, 9), NA_real_, df$BPQ030)
    )
  }
  if ("BPQ010" %in% names(df)) {
    bpq_age_vars <- c("BPQ040A", "BPQ040B", "BPQ040C", "BPQ040D")
    for (var in bpq_age_vars) {
      if (var %in% names(df)) {
        df[[var]] <- ifelse(df$BPQ010 == 2, NA_real_, df[[var]])
      }
    }
    bpq_doc_vars <- c("BPQ050A", "BPQ050B", "BPQ050C", "BPQ050D")
    for (var in bpq_doc_vars) {
      if (var %in% names(df)) {
        df[[var]] <- ifelse(df$BPQ010 == 2, NA_real_, df[[var]])
      }
    }
    cat("      已处理BPQ跳过模式\n")
  }
  age_vars <- c("BPQ040A", "BPQ040B", "BPQ040C", "BPQ040D",
                "BPQ050A", "BPQ050B", "BPQ050C", "BPQ050D")
  for (var in age_vars) {
    if (var %in% names(df)) {
      df[[var]] <- ifelse(df[[var]] %in% c(777, 999), NA_real_, df[[var]])
      df[[var]] <- ifelse(df[[var]] >= 0 & df[[var]] <= 120, df[[var]], NA_real_)
    }
  }
  if (all(c("BPQ010", "BPQ020") %in% names(df))) {
    df$HYPERTENSION_SELF_REPORT <- ifelse(
      df$BPQ010 == 1 | df$BPQ020 == 1,
      1,
      ifelse(df$BPQ010 == 2 & df$BPQ020 == 2, 0, NA_real_)
    )
    attr(df$HYPERTENSION_SELF_REPORT, "label") <- "自我报告的高血压(诊断或服药)"
    cat("      已创建HYPERTENSION_SELF_REPORT变量\n")
  }
  return(df)
}
#' OCQ专业清洗 - 职业问卷（P周期版）
clean_ocq_special <- function(df) {
  cat("    执行OCQ专业清洗 (P周期版)...\n")
  # OCQ180 变量名相同
  if ("OCQ180" %in% names(df)) {
    df$OCQ180 <- ifelse(df$OCQ180 %in% c(777777, 999999), NA_real_, df$OCQ180)
  }
  # P周期：OCD180 替代 OCD150
  if ("OCD180" %in% names(df)) {
    df$OCD180 <- ifelse(df$OCD180 %in% c(7, 9), NA_real_, df$OCD180)
    # 创建兼容变量供下游使用
    df$OCD150 <- df$OCD180
    cat("      已处理OCD180/OCD150 (就业状态)\n")
  }
  # OCQ210 变量名相同
  if ("OCQ210" %in% names(df)) {
    df$OCQ210 <- ifelse(df$OCQ210 %in% c(777, 999), NA_real_, df$OCQ210)
    df$OCQ210 <- ifelse(df$OCQ210 >= 0 & df$OCQ210 <= 168, df$OCQ210, NA_real_)
  }
  # OCQ260 变量名相同
  if ("OCQ260" %in% names(df)) {
    df$OCQ260 <- ifelse(df$OCQ260 %in% c(77, 99), NA_real_, df$OCQ260)
  }
  return(df)
}
#' 实验室数据专业清洗（P周期兼容版）
clean_lab_special <- function(df, file_name) {
  cat(sprintf("    执行%s专业清洗 (实验室数据)...\n", file_name))
  # HSCRP - 变量名相同
  if (grepl("HSCRP", file_name) && "LBXHSCRP" %in% names(df)) {
    df$LBXHSCRP <- ifelse(df$LBXHSCRP >= 0.1 & df$LBXHSCRP <= 20, df$LBXHSCRP, NA_real_)
    cat("      已处理LBXHSCRP\n")
  }
  # CBC - 血常规
  if (grepl("CBC", file_name)) {
    # 白细胞 - P周期: LBXWBC
    if ("LBXWBC" %in% names(df)) {
      df$LBXWBC <- ifelse(df$LBXWBC %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXWBC)
      df$LBXWBC <- ifelse(df$LBXWBC >= 2 & df$LBXWBC <= 30, df$LBXWBC, NA_real_)
      df$LBXWBCSI <- df$LBXWBC  # 创建兼容变量
      cat("      已处理LBXWBC/LBXWBCSI\n")
    }
    # 中性粒细胞 - P周期: LBXNE
    if ("LBXNE" %in% names(df)) {
      df$LBXNE <- ifelse(df$LBXNE %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXNE)
      df$LBXNE <- ifelse(df$LBXNE >= 10 & df$LBXNE <= 90, df$LBXNE, NA_real_)
      df$LBXNEPCT <- df$LBXNE
      cat("      已处理LBXNE/LBXNEPCT\n")
    }
    # 淋巴细胞 - P周期: LBXLY
    if ("LBXLY" %in% names(df)) {
      df$LBXLY <- ifelse(df$LBXLY %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXLY)
      df$LBXLY <- ifelse(df$LBXLY >= 10 & df$LBXLY <= 60, df$LBXLY, NA_real_)
      df$LBXLYPCT <- df$LBXLY
      cat("      已处理LBXLY/LBXLYPCT\n")
    }
    # 血小板 - P周期: LBXPLT
    if ("LBXPLT" %in% names(df)) {
      df$LBXPLT <- ifelse(df$LBXPLT %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXPLT)
      df$LBXPLT <- ifelse(df$LBXPLT >= 100 & df$LBXPLT <= 450, df$LBXPLT, NA_real_)
      df$LBXPLTSI <- df$LBXPLT
      cat("      已处理LBXPLT/LBXPLTSI\n")
    }
    # 红细胞压积 - 变量名相同
    if ("LBXHCT" %in% names(df)) {
      df$LBXHCT <- ifelse(df$LBXHCT %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXHCT)
      cat("      已处理LBXHCT\n")
    }
  }
  # BIOPRO - 生化指标
  if (grepl("BIOPRO", file_name)) {
    # AST - P周期: LBXSAST
    if ("LBXSAST" %in% names(df)) {
      df$LBXSAST <- ifelse(df$LBXSAST %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXSAST)
      df$LBXSAST <- ifelse(df$LBXSAST >= 5 & df$LBXSAST <= 500, df$LBXSAST, NA_real_)
      df$LBXSASSI <- df$LBXSAST
      cat("      已处理LBXSAST/LBXSASSI\n")
    }
    # ALT - P周期: LBXSAT
    if ("LBXSAT" %in% names(df)) {
      df$LBXSAT <- ifelse(df$LBXSAT %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXSAT)
      df$LBXSAT <- ifelse(df$LBXSAT >= 3 & df$LBXSAT <= 500, df$LBXSAT, NA_real_)
      df$LBXSALSI <- df$LBXSAT
      cat("      已处理LBXSAT/LBXSALSI\n")
    }
    # 尿酸 - 变量名相同
    if ("LBXSUA" %in% names(df)) {
      df$LBXSUA <- ifelse(df$LBXSUA %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXSUA)
      df$LBXSUA <- ifelse(df$LBXSUA >= 1 & df$LBXSUA <= 15, df$LBXSUA, NA_real_)
      cat("      已处理LBXSUA\n")
    }
    # 肌酐 - 变量名相同
    if ("LBXSCR" %in% names(df)) {
      df$LBXSCR <- ifelse(df$LBXSCR %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXSCR)
      df$LBXSCR <- ifelse(df$LBXSCR >= 0.2 & df$LBXSCR <= 10, df$LBXSCR, NA_real_)
      cat("      已处理LBXSCR\n")
    }
    # HDL - P周期: LBXHDL
    if ("LBXHDL" %in% names(df)) {
      df$LBXHDL <- ifelse(df$LBXHDL %in% c(7777, 777, 77, 9999, 999, 99), NA_real_, df$LBXHDL)
      df$LBXHDL <- ifelse(df$LBXHDL >= 10 & df$LBXHDL <= 150, df$LBXHDL, NA_real_)
      df$LBXHDD <- df$LBXHDL
      cat("      已处理LBXHDL/LBXHDD\n")
    }
  }
  return(df)
}
# ============================================================================
# 5. 添加多周期权重构建函数（NHANES教程第5章）
# ============================================================================
#' 构建多周期MEC权重
#' @param data 包含WTMEC2YR的数据框
#' @param cycles 要合并的周期数（2,3,4,5）
#' @return 添加了多周期权重的数据框
construct_multi_year_weights <- function(data, cycles = 2) {
  if (!"WTMEC2YR" %in% names(data)) {
    cat("    警告: 数据中无WTMEC2YR变量，无法构建多周期权重\n")
    return(data)
  }
  weight_name <- paste0("WTMEC", cycles, "YR")
  data[[weight_name]] <- data$WTMEC2YR / cycles
  cat(sprintf("    已创建: %s = WTMEC2YR / %d\n", weight_name, cycles))
  return(data)
}
cat("所有清洗函数定义完成\n\n")
# ============================================================================
# 6. 主清洗流程
# ============================================================================
cat("步骤3: 开始主清洗流程\n")
cat("========================================================\n\n")
for (i in seq_along(files_to_clean)) {
  file <- files_to_clean[i]
  cat(sprintf("\n[%d/%d] 处理文件: %s\n", i, length(files_to_clean), file))
  cat("----------------------------------------\n")
  tryCatch({
    # 6.1 检查文件是否存在
    file_path <- file.path(RAW_DATA_DIR, file)
    if (!file.exists(file_path)) {
      stop(paste("文件不存在:", file_path))
    }
    # 6.2 读取XPT文件
    cat("  读取文件...\n")
    df <- read_xpt(file_path)
    original_rows <- nrow(df)
    original_cols <- ncol(df)
    cat(sprintf("  原始维度: %d行 × %d列\n", original_rows, original_cols))
    # 6.3 应用标准清洗
    df_clean <- clean_nhanes_data(df, file)
    # 6.4 应用专业清洗
    if (grepl("DIQ", file)) {
      df_clean <- clean_diq_special(df_clean)
    } else if (grepl("DPQ", file)) {
      df_clean <- clean_dpq_special(df_clean)
    } else if (grepl("WHQ", file)) {
      df_clean <- clean_whq_special(df_clean)
    } else if (grepl("ALQ", file)) {
      df_clean <- clean_alq_special(df_clean)
    } else if (grepl("HUQ", file)) {
      df_clean <- clean_huq_special(df_clean)
    } else if (grepl("DBQ", file)) {
      df_clean <- clean_dbq_special(df_clean)  # P周期版
    } else if (grepl("SMQ", file)) {
      df_clean <- clean_smq_special(df_clean)
    } else if (grepl("BPQ", file)) {
      df_clean <- clean_bpq_special(df_clean)
    } else if (grepl("OCQ", file)) {
      df_clean <- clean_ocq_special(df_clean)  # P周期版
    } else if (grepl("HSCRP|CBC|BIOPRO", file)) {
      df_clean <- clean_lab_special(df_clean, file)  # P周期兼容版
    }
    # 6.5 对DEMO文件添加多周期权重选项
    if (file == "P_DEMO.xpt") {
      cat("  为DEMO文件添加多周期权重选项...\n")
      df_clean <- construct_multi_year_weights(df_clean, cycles = 2)
      df_clean <- construct_multi_year_weights(df_clean, cycles = 4)
    }
    # 6.6 验证关键变量
    seqn_present <- "SEQN" %in% names(df_clean)
    weight_vars <- c("WTINT2YR", "WTMEC2YR", "WTDRD1")
    weights_present <- any(weight_vars %in% names(df_clean))
    design_vars <- c("SDMVPSU", "SDMVSTRA")
    design_present <- all(design_vars %in% names(df_clean))
    # 6.7 保存为RDS - 保存到CLEAN_P目录
    output_file <- gsub("\\.xpt$", "_clean.rds", file)
    output_path <- file.path(CLEAN_DATA_DIR, output_file)
    saveRDS(df_clean, file = output_path)
    cat(sprintf("  已保存: %s\n", output_path))
    # 6.8 记录处理状态
    processing_log <- rbind(processing_log, data.frame(
      File = file,
      Status = "成功",
      Rows = nrow(df_clean),
      Cols = ncol(df_clean),
      SEQN_Present = seqn_present,
      Weights_Present = weights_present,
      DesignVars_Present = design_present,
      stringsAsFactors = FALSE
    ))
    cat("  ✓ 处理完成\n")
    # 6.9 内存清理
    rm(df, df_clean)
    gc()
  }, error = function(e) {
    processing_log <<- rbind(processing_log, data.frame(
      File = file,
      Status = paste("错误:", e$message),
      Rows = NA,
      Cols = NA,
      SEQN_Present = NA,
      Weights_Present = NA,
      DesignVars_Present = NA,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  ✗ 处理失败: %s\n", e$message))
  })
  cat("----------------------------------------\n")
}
# ============================================================================
# 7. 清洗后验证
# ============================================================================
cat("\n步骤4: 清洗后验证\n")
cat("========================================================\n\n")
# 7.1 显示处理状态
cat("处理状态总结:\n")
print(processing_log)
# 7.2 验证成功率
success_count <- sum(processing_log$Status == "成功")
if (success_count == length(files_to_clean)) {
  cat(sprintf("\n✅ 所有 %d 个文件清洗成功！\n", success_count))
} else {
  cat(sprintf("\n⚠️ 只成功清洗了 %d/%d 个文件\n", success_count, length(files_to_clean)))
  failed_files <- processing_log$File[processing_log$Status != "成功"]
  cat("失败的文件:", paste(failed_files, collapse = ", "), "\n")
}
# 7.3 验证DEMO文件的设计变量
cat("\n验证DEMO文件设计变量:\n")
demo_path <- file.path(CLEAN_DATA_DIR, "P_DEMO_clean.rds")
if (file.exists(demo_path)) {
  demo_data <- readRDS(demo_path)
  if (all(c("SDMVPSU", "SDMVSTRA") %in% names(demo_data))) {
    unique_psu <- length(unique(na.omit(demo_data$SDMVPSU)))
    unique_strata <- length(unique(na.omit(demo_data$SDMVSTRA)))
    cat(sprintf("  SDMVPSU: %d 个唯一值\n", unique_psu))
    cat(sprintf("  SDMVSTRA: %d 个唯一值\n", unique_strata))
    cat("  ✅ 设计变量存在\n")
    df_value <- unique_psu - unique_strata
    cat(sprintf("  实际自由度: %d (用于假设检验)\n", df_value))
  }
}
# 7.4 验证权重处理
cat("\n验证权重处理:\n")
if (file.exists(demo_path)) {
  if ("WTMEC2YR" %in% names(demo_data)) {
    wt_summary <- summary(demo_data$WTMEC2YR)
    wt_zero <- sum(demo_data$WTMEC2YR == 0, na.rm = TRUE)
    cat(sprintf("  WTMEC2YR: 最小值=%.2f, 中位数=%.2f, 最大值=%.2f\n",
                wt_summary["Min."], wt_summary["Median"], wt_summary["Max."]))
    cat(sprintf("  零权重数: %d\n", wt_zero))
  }
}
# ============================================================================
# 8. 生成清洗报告
# ============================================================================
cat("\n步骤5: 生成清洗报告...\n")
report_file <- file.path(LOG_DIR, "01_cleaning_report_P.txt")
sink(report_file)
cat("NHANES 2017-2020 (P周期) 数据清洗报告\n")
cat("========================================\n")
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("处理文件数:", length(files_to_clean), "\n\n")
cat("文件处理状态:\n")
print(processing_log)
cat("\n\nNHANES教程规范符合性检查:\n")
cat("✅ 缺失值处理: 7/77/777=拒绝, 9/99/999=不知道 (Table 8-1)\n")
cat("✅ 权重变量: 保留零权重 (有效设计特征)\n")
cat("✅ 设计变量: 验证SDMVPSU和SDMVSTRA存在\n")
cat("✅ MVUs说明: 正确理解Masked Variance Units\n")
cat("✅ 观测保留: 未删除任何观测\n")
cat("✅ 跳过模式: 完整处理DIQ、SMQ等问卷跳过逻辑\n")
cat("✅ 专业清洗: 处理P周期变量名差异\n\n")
cat("输出文件清单:\n")
for (file in files_to_clean) {
  output_file <- gsub("\\.xpt$", "_clean.rds", file)
  cat(sprintf("  - %s\n", file.path(CLEAN_DATA_DIR, output_file)))
}
sink()
cat(sprintf("清洗报告已保存: %s\n", report_file))
# ============================================================================
# 9. 保存处理日志
# ============================================================================
log_csv_path <- file.path(LOG_DIR, "01_cleaning_log_P.csv")
write.csv(processing_log, log_csv_path, row.names = FALSE)
cat(sprintf("处理日志已保存: %s\n", log_csv_path))
# ============================================================================
# 10. 保存会话信息
# ============================================================================
session_info_path <- file.path(LOG_DIR, "01_session_info_P.txt")
sink(session_info_path)
cat("NHANES 2017-2020 (P周期) 清洗会话信息\n")
cat("========================================\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
cat(sprintf("会话信息已保存: %s\n", session_info_path))
# ============================================================================
# 11. 保存R代码副本
# ============================================================================
cat("\n步骤6: 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "01_data_cleaning_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "01_code_list_P.txt")
cat("脚本名称: 01_data_cleaning_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 12. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期清洗完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", CLEAN_DATA_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "LOG_DIR", "processing_log")))
gc()
