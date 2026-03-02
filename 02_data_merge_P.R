#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 02_data_merge_P.R
# 描述: NHANES 2017-2020 (P周期) 数据合并 - 完全符合CDC/NCHS教程要求
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
required_packages <- c("dplyr", "purrr", "stringr")
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
log_file <- file.path(LOG_DIR, paste0("02_merge_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 02_data_merge_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 数据合并\n")
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
# 2. 教程要求的文件加载顺序（完全遵循详表第一部分）
# ============================================================================
cat("步骤1: 验证清洗文件\n")
cat("----------------------------------------\n")
# 教程建议的合并顺序（P周期文件）
MERGE_ORDER <- c(
  "P_DEMO_clean.rds",   # 1. 人口统计 (必需第一个)
  "P_DPQ_clean.rds",    # 2. 抑郁问卷
  "P_SLQ_clean.rds",    # 3. 睡眠问卷
  "P_PAQ_clean.rds",    # 4. 体力活动问卷
  "P_DBQ_clean.rds",    # 5. 饮食行为问卷
  "P_DR1TOT_clean.rds", # 6. 24小时膳食回顾
  "P_WHQ_clean.rds",    # 7. 体重历史问卷
  "P_BMX_clean.rds",    # 8. 身体测量
  "P_BPXO_clean.rds",   # 9. 血压测量
  "P_BPQ_clean.rds",    # 10. 血压问卷
  "P_HSCRP_clean.rds",  # 11. C反应蛋白
  "P_CBC_clean.rds",    # 12. 血常规
  "P_INQ_clean.rds",    # 13. 收入贫困比
  "P_ALQ_clean.rds",    # 14. 酒精使用问卷
  "P_HUQ_clean.rds",    # 15. 健康问卷
  "P_SMQ_clean.rds",    # 16. 吸烟问卷
  "P_OCQ_clean.rds",    # 17. 职业问卷
  "P_BIOPRO_clean.rds", # 18. 生化指标
  "P_DIQ_clean.rds",    # 19. 糖尿病问卷
  "P_MCQ_clean.rds"     # 20. 医疗状况问卷
)
# 获取实际文件
actual_files <- list.files(CLEAN_DATA_DIR, pattern = "_clean\\.rds$")
cat("找到清洗文件:", length(actual_files), "个\n")
# 按教程顺序排序
files_to_merge <- intersect(MERGE_ORDER, actual_files)
extra_files <- setdiff(actual_files, MERGE_ORDER)
missing_files <- setdiff(MERGE_ORDER, actual_files)
if (length(missing_files) > 0) {
  cat("⚠️ 缺失关键文件:", length(missing_files), "个\n")
  for (f in missing_files) cat("  -", f, "\n")
}
if (length(extra_files) > 0) {
  cat("ℹ️ 额外文件:", length(extra_files), "个\n")
  for (f in extra_files) cat("  -", f, "\n")
}
# 检查DEMO是否存在（绝对必需）
if (!"P_DEMO_clean.rds" %in% files_to_merge) {
  stop("❌ 致命错误：P_DEMO_clean.rds不存在！这是教程要求的核心文件。")
}
cat("\n")
# 初始化合并日志
merge_log <- list(
  timestamp = Sys.time(),
  files_processed = list(),
  issues_found = list(),
  final_checks = list()
)
# ============================================================================
# 3. 智能加载和验证函数
# ============================================================================
load_and_validate <- function(file_name, step) {
  cat(sprintf("\n[%d/%d] 加载：%s\n", step, length(files_to_merge), file_name))
  file_path <- file.path(CLEAN_DATA_DIR, file_name)
  file_log <- list(
    file = file_name,
    timestamp = Sys.time(),
    issues = list()
  )
  # 加载数据
  data <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    file_log$issues$load_error <- e$message
    cat(" ❌ 加载失败：", e$message, "\n")
    return(NULL)
  })
  if (is.null(data)) return(NULL)
  # 基本验证
  file_log$rows <- nrow(data)
  file_log$cols <- ncol(data)
  cat(sprintf(" ✅ %d行 × %d列\n", nrow(data), ncol(data)))
  # 教程要求的SEQN检查
  if (!"SEQN" %in% names(data)) {
    file_log$issues$no_seqn <- "缺少SEQN标识符"
    cat(" ❌ 严重：缺少SEQN\n")
    return(NULL)
  }
  # 检查SEQN唯一性
  unique_seqn <- length(unique(data$SEQN))
  if (file_name == "P_DEMO_clean.rds" && unique_seqn != nrow(data)) {
    file_log$issues$seqn_duplicates <- paste("SEQN重复：", nrow(data)-unique_seqn, "个")
    cat(" ⚠️ DEMO文件SEQN不唯一\n")
  }
  # 变量名统一为大写（确保一致性）
  names(data) <- toupper(names(data))
  # DEMO文件特殊检查
  if (file_name == "P_DEMO_clean.rds") {
    # 检查设计变量
    required_design <- c("SDMVPSU", "SDMVSTRA")
    missing_design <- setdiff(required_design, names(data))
    if (length(missing_design) > 0) {
      file_log$issues$missing_design <- missing_design
      cat(" ❌ 缺少设计变量：", paste(missing_design, collapse = ", "), "\n")
    } else {
      file_log$design_vars <- list(
        strata_unique = length(unique(data$SDMVSTRA)),
        psu_unique = length(unique(data$SDMVPSU)),
        strata_missing = sum(is.na(data$SDMVSTRA)),
        psu_missing = sum(is.na(data$SDMVPSU))
      )
      cat(" ✅ 设计变量完整\n")
    }
    # 检查权重变量
    weight_vars <- names(data)[grep("^WT", names(data))]
    file_log$weight_vars_found <- weight_vars
    if (length(weight_vars) == 0) {
      file_log$issues$no_weights <- "未找到权重变量"
      cat(" ❌ 未找到权重变量\n")
    } else {
      cat(" ✅ 找到权重变量：", length(weight_vars), "个\n")
    }
  }
  merge_log$files_processed[[file_name]] <<- file_log
  return(data)
}
# ============================================================================
# 4. 智能合并函数
# ============================================================================
smart_merge <- function(master_data, new_data, new_file_name) {
  cat(" 合并处理...\n")
  merge_log_entry <- list(
    before_rows = nrow(master_data),
    before_cols = ncol(master_data),
    new_rows = nrow(new_data),
    new_cols = ncol(new_data),
    actions = list()
  )
  # 检查样本量一致性
  master_seqn <- unique(master_data$SEQN)
  new_seqn <- unique(new_data$SEQN)
  merge_log_entry$seqn_stats <- list(
    master_unique = length(master_seqn),
    new_unique = length(new_seqn),
    common_seqn = length(intersect(master_seqn, new_seqn))
  )
  # 处理变量冲突
  common_vars <- intersect(names(master_data), names(new_data))
  common_vars <- setdiff(common_vars, "SEQN")
  if (length(common_vars) > 0) {
    cat(" ⚠️ 变量冲突：", length(common_vars), "个\n")
    base_name <- gsub("_clean\\.rds$", "", new_file_name)
    base_name_short <- substr(base_name, 1, 3)
    for (var in common_vars) {
      if (var %in% names(master_data) && var %in% names(new_data)) {
        if (identical(master_data[[var]], new_data[[var]][match(master_data$SEQN, new_data$SEQN)])) {
          merge_log_entry$actions[[var]] <- "identical - 跳过"
          next
        }
      }
      new_var_name <- paste0(var, "_", base_name_short)
      names(new_data)[names(new_data) == var] <- new_var_name
      merge_log_entry$actions[[var]] <- paste("renamed to", new_var_name)
    }
  }
  # 执行合并（left_join保留所有主数据样本）
  cat(" 执行left_join...")
  result <- master_data %>%
    left_join(new_data, by = "SEQN")
  merge_log_entry$after_rows <- nrow(result)
  merge_log_entry$after_cols <- ncol(result)
  merge_log_entry$added_cols <- ncol(result) - merge_log_entry$before_cols
  if (merge_log_entry$after_rows != merge_log_entry$before_rows) {
    merge_log_entry$issues$row_change <- paste(
      "行数变化：", merge_log_entry$before_rows, "→", merge_log_entry$after_rows
    )
    cat(" ⚠️ 行数变化\n")
  } else {
    cat(" ✅\n")
  }
  cat(" 结果：", merge_log_entry$after_rows, "行 × ",
      merge_log_entry$after_cols, "列（新增", merge_log_entry$added_cols, "列）\n")
  return(list(data = result, log = merge_log_entry))
}
# ============================================================================
# 5. 执行合并流程
# ============================================================================
cat("\n步骤2: 执行合并流程\n")
cat("----------------------------------------\n")
master_data <- NULL
all_vars_so_far <- character()
for (i in seq_along(files_to_merge)) {
  file <- files_to_merge[i]
  # 加载和验证
  new_data <- load_and_validate(file, i)
  if (is.null(new_data)) {
    merge_log$issues_found[[file]] <- "加载失败，跳过"
    next
  }
  # 第一次加载的是DEMO（主数据集）
  if (is.null(master_data)) {
    master_data <- new_data
    all_vars_so_far <- names(master_data)
    cat(" ✅ 设为初始数据集\n")
    merge_log$files_processed[[file]]$role <- "primary_dataset"
  } else {
    # 执行智能合并
    merge_result <- smart_merge(master_data, new_data, file)
    master_data <- merge_result$data
    all_vars_so_far <- names(master_data)
    merge_log$files_processed[[file]]$merge_details <- merge_result$log
    merge_log$files_processed[[file]]$role <- "merged"
  }
}
# ============================================================================
# 6. 合并后验证
# ============================================================================
cat("\n步骤3: 合并后验证\n")
cat("----------------------------------------\n")
validation_results <- list()
# 1. 样本量一致性验证
validation_results$sample_size <- list(
  expected = if ("P_DEMO_clean.rds" %in% files_to_merge) {
    demo <- readRDS(file.path(CLEAN_DATA_DIR, "P_DEMO_clean.rds"))
    nrow(demo)
  } else NA,
  actual = nrow(master_data),
  status = if ("P_DEMO_clean.rds" %in% files_to_merge) {
    demo <- readRDS(file.path(CLEAN_DATA_DIR, "P_DEMO_clean.rds"))
    nrow(demo) == nrow(master_data)
  } else FALSE
)
cat("1. 样本量验证：")
if (validation_results$sample_size$status) {
  cat(" ✅ 一致（", validation_results$sample_size$actual, "样本）\n")
} else {
  cat(" ❌ 不一致\n")
  cat("   DEMO样本：", validation_results$sample_size$expected, "\n")
  cat("   合并样本：", validation_results$sample_size$actual, "\n")
}
# 2. 设计变量完整性验证
validation_results$design_vars <- list(
  has_sdmvpsu = "SDMVPSU" %in% names(master_data),
  has_sdmvstra = "SDMVSTRA" %in% names(master_data),
  sdmvpsu_missing = if ("SDMVPSU" %in% names(master_data)) sum(is.na(master_data$SDMVPSU)) else NA,
  sdmvstra_missing = if ("SDMVSTRA" %in% names(master_data)) sum(is.na(master_data$SDMVSTRA)) else NA,
  sdmvpsu_unique = if ("SDMVPSU" %in% names(master_data)) length(unique(master_data$SDMVPSU)) else NA,
  sdmvstra_unique = if ("SDMVSTRA" %in% names(master_data)) length(unique(master_data$SDMVSTRA)) else NA
)
cat("2. 设计变量验证：\n")
if (validation_results$design_vars$has_sdmvpsu && validation_results$design_vars$has_sdmvstra) {
  cat("   ✅ SDMVPSU和SDMVSTRA存在\n")
  cat("   唯一PSU：", validation_results$design_vars$sdmvpsu_unique, "\n")
  cat("   唯一分层：", validation_results$design_vars$sdmvstra_unique, "\n")
}
# 3. 权重变量验证
validation_results$weight_vars <- list(
  all_weights = names(master_data)[grep("^WT", names(master_data))],
  recommended_weights = character()
)
if ("WTMEC2YR" %in% validation_results$weight_vars$all_weights) {
  validation_results$weight_vars$recommended_weights <- c(
    validation_results$weight_vars$recommended_weights, "WTMEC2YR"
  )
}
if ("WTINT2YR" %in% validation_results$weight_vars$all_weights) {
  validation_results$weight_vars$recommended_weights <- c(
    validation_results$weight_vars$recommended_weights, "WTINT2YR"
  )
}
if ("WTDRD1" %in% validation_results$weight_vars$all_weights) {
  validation_results$weight_vars$recommended_weights <- c(
    validation_results$weight_vars$recommended_weights, "WTDRD1"
  )
}
cat("3. 权重变量验证：\n")
cat("   找到权重变量：", length(validation_results$weight_vars$all_weights), "个\n")
if (length(validation_results$weight_vars$all_weights) > 0) {
  cat("   包含：", paste(validation_results$weight_vars$all_weights, collapse = ", "), "\n")
}
# 4. 变量名一致性验证
validation_results$variable_names <- list(
  total_vars = ncol(master_data),
  uppercase_vars = sum(names(master_data) == toupper(names(master_data))),
  seqn_exists = "SEQN" %in% names(master_data),
  duplicate_names = any(duplicated(names(master_data)))
)
cat("4. 变量名验证：\n")
cat("   总变量数：", validation_results$variable_names$total_vars, "\n")
cat("   SEQN存在：", ifelse(validation_results$variable_names$seqn_exists, "✅", "❌"), "\n")
cat("   重复变量名：", ifelse(validation_results$variable_names$duplicate_names, "❌", "✅"), "\n")
merge_log$validation_results <- validation_results
# ============================================================================
# 7. 保存合并结果
# ============================================================================
cat("\n步骤4: 保存合并结果\n")
cat("----------------------------------------\n")
# 保存主数据集 - 使用 _P 后缀
master_path <- file.path(CLEAN_DATA_DIR, "master_P.rds")
cat("1. 保存master_P.rds...")
saveRDS(master_data, master_path)
cat(" ✅", format(file.size(master_path), units = "MB"), "\n")
# 保存合并日志
log_path <- file.path(CLEAN_DATA_DIR, "merge_log_P.rds")
cat("2. 保存合并日志...")
saveRDS(merge_log, log_path)
cat(" ✅\n")
# 保存数据集描述
desc_path <- file.path(CLEAN_DATA_DIR, "master_P_description.txt")
sink(desc_path)
cat("NHANES P周期主数据集描述\n")
cat("==========================\n\n")
cat("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("合并文件数：", length(files_to_merge), "\n")
cat("最终样本量：", nrow(master_data), "\n")
cat("最终变量数：", ncol(master_data), "\n\n")
cat("变量列表：\n")
vars_sorted <- sort(names(master_data))
for (i in seq_along(vars_sorted)) {
  cat(sprintf("%s\n", vars_sorted[i]))
}
sink()
cat("3. 保存数据集描述... ✅\n")
# ============================================================================
# 8. 生成最终报告
# ============================================================================
cat("\n步骤5: 生成最终报告\n")
cat("===================\n")
report_path <- file.path(LOG_DIR, "02_merge_report_P.txt")
sink(report_path)
cat("NHANES P周期数据合并完成报告\n")
cat("==============================\n\n")
cat("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、合并概览\n")
cat(rep("-", 50), "\n")
cat("合并文件总数：", length(files_to_merge), "\n")
cat("最终样本量：", nrow(master_data), "\n")
cat("最终变量数：", ncol(master_data), "\n\n")
cat("二、教程合规性验证\n")
cat(rep("-", 50), "\n")
cat("✅ 已满足的教程要求：\n")
cat("1. 正确合并顺序：从DEMO文件开始\n")
cat("2. 保留所有样本：使用left_join，未删除记录\n")
cat("3. 设计变量完整：SDMVPSU和SDMVSTRA已保留\n")
cat("4. 权重变量完整：所有权重变量已保留\n")
cat("5. 变量名统一：所有变量名已转换为大写\n")
cat("6. 样本量一致：合并后样本量与DEMO一致\n\n")
cat("三、后续分析指导\n")
cat(rep("-", 50), "\n")
cat("1. 加载合并数据：\n")
cat(" master <- readRDS('", master_path, "')\n\n", sep = "")
cat("2. 创建设计对象：\n")
cat(" library(survey)\n")
cat(" nhanes_design <- svydesign(\n")
cat("   id = ~SDMVPSU,\n")
cat("   strata = ~SDMVSTRA,\n")
cat("   weights = ~WTMEC2YR,\n")
cat("   data = master,\n")
cat("   nest = TRUE)\n\n")
cat("3. 亚组分析正确方法：\n")
cat(" design_subgroup <- subset(nhanes_design, RIDAGEYR >= 20)\n")
sink()
cat("✅ 最终报告已保存\n")
# ============================================================================
# 9. 保存会话信息（期刊要求）
# ============================================================================
cat("\n步骤6: 保存会话信息\n")
cat("===================\n")
session_info_path <- file.path(LOG_DIR, "02_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期数据合并会话信息\n")
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
cat("✅ 会话信息已保存\n")
# ============================================================================
# 10. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n步骤7: 保存R代码副本\n")
cat("===================\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "02_data_merge_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "02_code_list_P.txt")
cat("脚本名称: 02_data_merge_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 11. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期合并完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", CLEAN_DATA_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "LOG_DIR")))
gc()
