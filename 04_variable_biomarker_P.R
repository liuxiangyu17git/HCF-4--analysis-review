#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 04_variable_biomarker_P.R
# 描述: NHANES 2017-2020 (P周期) 生物标志物构建 - 阶段二
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: OSF https://osf.io/xxxxx (DOI: 10.17605/OSF.IO/XXXXX)
# 对应研究计划: 第四部分 4.2 变量构建总览
# 对应变量详表: 第三部分 构建阶段详表（阶段二 - 生物标志物）
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
required_packages <- c("dplyr", "tidyr", "survey", "haven", "srvyr")
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
log_file <- file.path(LOG_DIR, paste0("04_biomarker_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 04_variable_biomarker_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 生物标志物构建\n")
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
# 2. 数据加载和验证
# ============================================================================
cat("1. 加载主数据集并验证完整性...\n")
master_file <- file.path(CLEAN_DATA_DIR, "master_P.rds")
if(!file.exists(master_file)) {
  stop("❌ 错误：master_P.rds不存在！请先运行03_variable_basic_P.R")
}
file_info <- file.info(master_file)
cat("数据文件：", master_file, "\n")
cat("文件大小：", round(file_info$size/1024/1024, 2), "MB\n")
cat("修改时间：", file_info$mtime, "\n\n")
master <- readRDS(master_file)
cat("数据完整性验证：\n")
cat("样本量：", nrow(master), "\n")
cat("变量数：", ncol(master), "\n")
# 检查关键必需变量
essential_vars <- c("SEQN", "SDMVPSU", "SDMVSTRA", "WTMECPRP",
                    "RIDAGEYR", "RIAGENDR", "RIDEXPRG")
missing_essential <- setdiff(essential_vars, names(master))
if(length(missing_essential) > 0) {
  cat("❌ 严重错误：缺失必需变量：", paste(missing_essential, collapse = ", "), "\n")
  stop("无法继续分析。")
} else {
  cat("✅ 所有必需变量存在\n")
}
# 检查调查周期
if("SDDSRVYR" %in% names(master)) {
  survey_years <- unique(master$SDDSRVYR)
  year_mapping <- c("1" = "1999-2000", "2" = "2001-2002", "3" = "2003-2004",
                    "4" = "2005-2006", "5" = "2007-2008", "6" = "2009-2010",
                    "7" = "2011-2012", "8" = "2013-2014", "9" = "2015-2016",
                    "10" = "2017-2018", "11" = "2019-2020")
  actual_years <- year_mapping[as.character(survey_years)]
  cat("调查周期：", paste(actual_years, collapse = ", "), "\n")
}
cat("\n")
# ============================================================================
# 3. 创建亚组指示变量（NHANES教程完全合规）
# ============================================================================
cat("2. 创建亚组指示变量（教程完全合规）...\n")
data <- master %>%
  mutate(
    adult_indicator = ifelse(RIDAGEYR >= 18, 1, 0),
    is_pregnant_excluded = case_when(
      RIDEXPRG == 1 ~ 1,
      RIDEXPRG == 2 ~ 0,
      RIDEXPRG == 3 ~ 0,
      RIDEXPRG %in% c(7, 9) ~ NA_real_,
      is.na(RIDEXPRG) ~ 0,
      TRUE ~ NA_real_
    ),
    target_population = ifelse(adult_indicator == 1 & is_pregnant_excluded == 0, 1, 0),
    in_analysis = target_population
  )
cat("详细亚组分布：\n")
cat("总样本：", nrow(data), "\n")
cat("年龄分布：\n")
cat(" • <18岁：", sum(data$RIDAGEYR < 18, na.rm = TRUE),
    " (", round(mean(data$RIDAGEYR < 18, na.rm = TRUE)*100, 1), "%)\n")
cat(" • ≥18岁：", sum(data$adult_indicator == 1, na.rm = TRUE),
    " (", round(mean(data$adult_indicator == 1, na.rm = TRUE)*100, 1), "%)\n\n")
cat("目标人群（成人非孕妇）：", sum(data$target_population == 1, na.rm = TRUE),
    " (", round(mean(data$target_population == 1, na.rm = TRUE)*100, 1), "%)\n\n")
# ============================================================================
# 4. 构建生物标志物变量（P周期兼容版）
# ============================================================================
cat("3. 构建生物标志物变量...\n")
# 4.1 肝功能指标 - P周期使用LBXSAST/LBXSAT
cat(" 3.1 肝功能指标...\n")
# AST - P周期: LBXSAST
if("LBXSAST" %in% names(data)) {
  data <- data %>%
    mutate(
      ast_ui_l = case_when(
        LBXSAST %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXSAST < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXSAST)
      )
    )
  # 创建兼容变量
  data$LBXSASSI <- data$ast_ui_l
} else if("LBXSASSI" %in% names(data)) {
  data <- data %>%
    mutate(
      ast_ui_l = case_when(
        LBXSASSI %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXSASSI < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXSASSI)
      )
    )
} else {
  data$ast_ui_l <- NA_real_
}
# ALT - P周期: LBXSAT
if("LBXSAT" %in% names(data)) {
  data <- data %>%
    mutate(
      alt_ui_l = case_when(
        LBXSAT %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXSAT < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXSAT)
      )
    )
  # 创建兼容变量
  data$LBXSALSI <- data$alt_ui_l
} else if("LBXSALSI" %in% names(data)) {
  data <- data %>%
    mutate(
      alt_ui_l = case_when(
        LBXSALSI %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXSALSI < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXSALSI)
      )
    )
} else {
  data$alt_ui_l <- NA_real_
}
data <- data %>%
  mutate(
    ast_alt_ratio = case_when(
      !is.na(ast_ui_l) & !is.na(alt_ui_l) & alt_ui_l > 0 ~ ast_ui_l / alt_ui_l,
      TRUE ~ NA_real_
    ),
    elevated_alt_aga2022 = case_when(
      RIAGENDR == 1 & !is.na(alt_ui_l) & alt_ui_l > 30 ~ 1,
      RIAGENDR == 2 & !is.na(alt_ui_l) & alt_ui_l > 19 ~ 1,
      !is.na(alt_ui_l) ~ 0,
      TRUE ~ NA_real_
    ),
    elevated_alt_traditional = case_when(
      !is.na(alt_ui_l) & alt_ui_l > 40 ~ 1,
      !is.na(alt_ui_l) ~ 0,
      TRUE ~ NA_real_
    ),
    liver_injury = case_when(
      elevated_alt_traditional == 1 ~ 1,
      !is.na(elevated_alt_traditional) ~ 0,
      TRUE ~ NA_real_
    )
  )
# 4.2 尿酸指标 - LBXSUA变量名相同
cat(" 3.2 尿酸指标...\n")
if("LBXSUA" %in% names(data)) {
  data <- data %>%
    mutate(
      uric_acid_mgdl = case_when(
        LBXSUA %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXSUA < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXSUA)
      ),
      hyperuricemia = case_when(
        !is.na(uric_acid_mgdl) & !is.na(RIAGENDR) &
          ((RIAGENDR == 1 & uric_acid_mgdl > 7.0) |
           (RIAGENDR == 2 & uric_acid_mgdl > 6.0)) ~ 1,
        !is.na(uric_acid_mgdl) & !is.na(RIAGENDR) ~ 0,
        TRUE ~ NA_real_
      )
    )
} else {
  data$uric_acid_mgdl <- NA_real_
  data$hyperuricemia <- NA_real_
}
# 4.3 非空腹血脂指标 - P周期使用LBXHDL
cat(" 3.3 非空腹血脂指标...\n")
if("LBXHDL" %in% names(data)) {
  data <- data %>%
    mutate(
      hdl_cholesterol_mgdl = case_when(
        LBXHDL %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXHDL < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXHDL)
      ),
      low_hdl = case_when(
        !is.na(hdl_cholesterol_mgdl) & !is.na(RIAGENDR) &
          ((RIAGENDR == 1 & hdl_cholesterol_mgdl < 40) |
           (RIAGENDR == 2 & hdl_cholesterol_mgdl < 50)) ~ 1,
        !is.na(hdl_cholesterol_mgdl) & !is.na(RIAGENDR) ~ 0,
        TRUE ~ NA_real_
      )
    )
  # 创建兼容变量
  data$LBXHDD <- data$hdl_cholesterol_mgdl
  cat(" ✅ HDL胆固醇构建完成（来自LBXHDL）\n")
} else if("LBXHDD" %in% names(data)) {
  data <- data %>%
    mutate(
      hdl_cholesterol_mgdl = case_when(
        LBXHDD %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXHDD < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXHDD)
      ),
      low_hdl = case_when(
        !is.na(hdl_cholesterol_mgdl) & !is.na(RIAGENDR) &
          ((RIAGENDR == 1 & hdl_cholesterol_mgdl < 40) |
           (RIAGENDR == 2 & hdl_cholesterol_mgdl < 50)) ~ 1,
        !is.na(hdl_cholesterol_mgdl) & !is.na(RIAGENDR) ~ 0,
        TRUE ~ NA_real_
      )
    )
  cat(" ✅ HDL胆固醇构建完成（来自LBXHDD）\n")
} else {
  data$hdl_cholesterol_mgdl <- NA_real_
  data$low_hdl <- NA_real_
  cat(" ⚠️ HDL变量不存在，跳过HDL胆固醇构建\n")
}
# 4.4 炎症标志物
cat(" 3.4 炎症标志物...\n")
# C反应蛋白 - LBXHSCRP变量名相同
if("LBXHSCRP" %in% names(data)) {
  data <- data %>%
    mutate(
      hs_crp_mgl = case_when(
        LBXHSCRP %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXHSCRP < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXHSCRP)
      ),
      high_hs_crp = case_when(
        !is.na(hs_crp_mgl) & hs_crp_mgl > 3.0 ~ 1,
        !is.na(hs_crp_mgl) ~ 0,
        TRUE ~ NA_real_
      )
    )
} else {
  data$hs_crp_mgl <- NA_real_
}
# NLR - P周期使用LBXNE/LBXLY
if(all(c("LBXNE", "LBXLY") %in% names(data))) {
  data <- data %>%
    mutate(
      neutrophil_percent = case_when(
        LBXNE %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXNE < 0 | LBXNE > 100 ~ NA_real_,
        TRUE ~ as.numeric(LBXNE)
      ),
      lymphocyte_percent = case_when(
        LBXLY %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXLY < 0 | LBXLY > 100 ~ NA_real_,
        TRUE ~ as.numeric(LBXLY)
      ),
      nlr_ratio = case_when(
        !is.na(neutrophil_percent) & !is.na(lymphocyte_percent) &
          lymphocyte_percent > 0 ~ neutrophil_percent / lymphocyte_percent,
        TRUE ~ NA_real_
      )
    )
  # 创建兼容变量
  data$LBXNEPCT <- data$neutrophil_percent
  data$LBXLYPCT <- data$lymphocyte_percent
} else if(all(c("LBXNEPCT", "LBXLYPCT") %in% names(data))) {
  data <- data %>%
    mutate(
      neutrophil_percent = case_when(
        LBXNEPCT %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXNEPCT < 0 | LBXNEPCT > 100 ~ NA_real_,
        TRUE ~ as.numeric(LBXNEPCT)
      ),
      lymphocyte_percent = case_when(
        LBXLYPCT %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXLYPCT < 0 | LBXLYPCT > 100 ~ NA_real_,
        TRUE ~ as.numeric(LBXLYPCT)
      ),
      nlr_ratio = case_when(
        !is.na(neutrophil_percent) & !is.na(lymphocyte_percent) &
          lymphocyte_percent > 0 ~ neutrophil_percent / lymphocyte_percent,
        TRUE ~ NA_real_
      )
    )
} else {
  data$neutrophil_percent <- NA_real_
  data$lymphocyte_percent <- NA_real_
  data$nlr_ratio <- NA_real_
}
# 白细胞计数 - P周期使用LBXWBC
if("LBXWBC" %in% names(data)) {
  data <- data %>%
    mutate(
      white_blood_cells_x10e9_l = case_when(
        LBXWBC %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXWBC < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXWBC)
      ),
      leukocytosis = case_when(
        !is.na(white_blood_cells_x10e9_l) & white_blood_cells_x10e9_l > 11.0 ~ 1,
        !is.na(white_blood_cells_x10e9_l) ~ 0,
        TRUE ~ NA_real_
      )
    )
  # 创建兼容变量
  data$LBXWBCSI <- data$white_blood_cells_x10e9_l
} else if("LBXWBCSI" %in% names(data)) {
  data <- data %>%
    mutate(
      white_blood_cells_x10e9_l = case_when(
        LBXWBCSI %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXWBCSI < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXWBCSI)
      ),
      leukocytosis = case_when(
        !is.na(white_blood_cells_x10e9_l) & white_blood_cells_x10e9_l > 11.0 ~ 1,
        !is.na(white_blood_cells_x10e9_l) ~ 0,
        TRUE ~ NA_real_
      )
    )
} else {
  data$white_blood_cells_x10e9_l <- NA_real_
}
data <- data %>%
  mutate(
    elevated_ast = case_when(
      !is.na(ast_ui_l) & ast_ui_l > 40 ~ 1,
      !is.na(ast_ui_l) ~ 0,
      TRUE ~ NA_real_
    )
  )
# 4.5 肾功能指标 - LBXSCR变量名相同
cat(" 3.5 肾功能指标...\n")
if("LBXSCR" %in% names(data)) {
  data <- data %>%
    mutate(
      serum_creatinine_mgdl = case_when(
        LBXSCR %in% c(7777, 777, 77, 9999, 999, 99, 77777, 99999) ~ NA_real_,
        LBXSCR < 0 ~ NA_real_,
        TRUE ~ as.numeric(LBXSCR)
      )
    )
  if(all(c("serum_creatinine_mgdl", "RIDAGEYR", "RIAGENDR") %in% names(data))) {
    data <- data %>%
      mutate(
        scr = serum_creatinine_mgdl,
        age = RIDAGEYR,
        sex = RIAGENDR,
        kappa = case_when(
          sex == 2 & scr <= 0.7 ~ 0.7,
          sex == 2 & scr > 0.7 ~ 0.7,
          sex == 1 & scr <= 0.9 ~ 0.9,
          sex == 1 & scr > 0.9 ~ 0.9,
          TRUE ~ NA_real_
        ),
        alpha = case_when(
          sex == 2 & scr <= 0.7 ~ -0.241,
          sex == 2 & scr > 0.7 ~ -1.200,
          sex == 1 & scr <= 0.9 ~ -0.302,
          sex == 1 & scr > 0.9 ~ -1.200,
          TRUE ~ NA_real_
        ),
        egfr_ckdepi_2021 = case_when(
          !is.na(scr) & !is.na(age) & !is.na(sex) & !is.na(kappa) & !is.na(alpha) ~
            142 * (scr/kappa)^alpha * 0.9938^age,
          TRUE ~ NA_real_
        ),
        egfr_ckdepi_2021 = ifelse(!is.na(egfr_ckdepi_2021),
                                   pmin(pmax(egfr_ckdepi_2021, 5), 200),
                                   NA_real_),
        ckd_stage_kdigo = case_when(
          !is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 >= 90 ~ "G1",
          !is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 >= 60 ~ "G2",
          !is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 >= 45 ~ "G3a",
          !is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 >= 30 ~ "G3b",
          !is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 >= 15 ~ "G4",
          !is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 < 15 ~ "G5",
          TRUE ~ NA_character_
        ),
        ckd_flag = ifelse(!is.na(egfr_ckdepi_2021) & egfr_ckdepi_2021 < 60, 1, 0)
      ) %>%
      select(-scr, -age, -sex, -kappa, -alpha)
    cat(" ✅ eGFR计算完成\n")
  }
} else {
  cat(" ⚠️ LBXSCR不存在，跳过肾功能指标构建\n")
}
# 4.6 临床诊断变量
cat(" 3.6 临床诊断变量...\n")
if("DIQ010" %in% names(data)) {
  data <- data %>%
    mutate(
      diabetes_doctor = case_when(
        DIQ010 == 1 ~ 1,
        DIQ010 == 2 ~ 0,
        DIQ010 %in% c(7, 9) ~ NA_real_,
        TRUE ~ NA_real_
      )
    )
} else {
  data$diabetes_doctor <- NA_real_
}
if("BPQ050A" %in% names(data)) {
  data <- data %>%
    mutate(
      hypertension_doctor = case_when(
        BPQ050A == 1 ~ 1,
        BPQ050A == 2 ~ 0,
        BPQ050A %in% c(7, 9) ~ NA_real_,
        TRUE ~ NA_real_
      )
    )
} else {
  data$hypertension_doctor <- NA_real_
}
if("BPQ080" %in% names(data)) {
  data <- data %>%
    mutate(
      high_cholesterol_doctor = case_when(
        BPQ080 == 1 ~ 1,
        BPQ080 == 2 ~ 0,
        BPQ080 %in% c(7, 9) ~ NA_real_,
        TRUE ~ NA_real_
      )
    )
} else {
  data$high_cholesterol_doctor <- NA_real_
}
# 4.7 补充变量
cat(" 3.7 补充变量...\n")
# 白细胞计数（已有，这里只是确保变量名一致）
if("white_blood_cells_x10e9_l" %in% names(data)) {
  data$white_blood_cells <- data$white_blood_cells_x10e9_l
  cat(" ✅ 白细胞计数变量已就绪\n")
}
data <- data %>%
  mutate(
    elevated_ast = case_when(
      !is.na(ast_ui_l) & ast_ui_l > 40 ~ 1,
      !is.na(ast_ui_l) ~ 0,
      TRUE ~ NA_real_
    ),
    elevated_alt_traditional = case_when(
      !is.na(alt_ui_l) & alt_ui_l > 40 ~ 1,
      !is.na(alt_ui_l) ~ 0,
      TRUE ~ NA_real_
    ),
    high_nlr = case_when(
      !is.na(nlr_ratio) & nlr_ratio > 3.0 ~ 1,
      !is.na(nlr_ratio) ~ 0,
      TRUE ~ NA_real_
    ),
    liver_injury = case_when(
      elevated_ast == 1 | elevated_alt_traditional == 1 ~ 1,
      !is.na(elevated_ast) & !is.na(elevated_alt_traditional) ~ 0,
      TRUE ~ NA_real_
    )
  )
if("low_hdl" %in% names(data)) {
  data <- data %>%
    mutate(
      dyslipidemia = case_when(
        low_hdl == 1 ~ 1,
        low_hdl == 0 ~ 0,
        TRUE ~ NA_real_
      )
    )
  cat(" ✅ dyslipidemia构建完成\n")
} else {
  data$dyslipidemia <- NA_real_
}
# 4.8 复合代谢评分
cat(" 3.8 复合代谢评分...\n")
has_egfr <- "egfr_ckdepi_2021" %in% names(data)
data <- data %>%
  mutate(
    hcf_a_score = 
      (BMXBMI >= 30) +
      (hypertension_doctor == 1) +
      (diabetes_doctor == 1) +
      (high_cholesterol_doctor == 1),
    hcf_a_layer = ifelse(hcf_a_score >= 2, 1, 0),
    has_valid_mec = case_when(
      !is.na(WTMECPRP) & WTMECPRP > 0 ~ 1,
      TRUE ~ 0
    ),
    has_biomarker_data = case_when(
      !is.na(ast_ui_l) | !is.na(alt_ui_l) |
        {if(has_egfr) !is.na(egfr_ckdepi_2021) else FALSE} |
        !is.na(hs_crp_mgl) ~ 1,
      TRUE ~ 0
    )
  )
cat(" ✅ 复合代谢评分完成\n\n")
# ============================================================================
# 5. 创建分析数据集 - 完整版和核心版
# ============================================================================
cat("4. 创建最终分析数据集...\n")
# 完整版变量列表（48个）
full_vars <- c(
  "SEQN", "SDMVPSU", "SDMVSTRA", "WTMECPRP",
  "adult_indicator", "is_pregnant_excluded", "target_population", "in_analysis",
  "has_valid_mec", "has_biomarker_data",
  "RIDAGEYR", "RIAGENDR", "RIDEXPRG",
  "ast_ui_l", "alt_ui_l", "ast_alt_ratio", "elevated_alt_aga2022",
  "elevated_alt_traditional", "elevated_ast", "liver_injury",
  "uric_acid_mgdl", "hyperuricemia",
  "hdl_cholesterol_mgdl", "low_hdl", "dyslipidemia",
  "hs_crp_mgl", "high_hs_crp", "white_blood_cells_x10e9_l", "leukocytosis",
  "neutrophil_percent", "lymphocyte_percent", "nlr_ratio", "high_nlr",
  "serum_creatinine_mgdl", "egfr_ckdepi_2021", "ckd_stage_kdigo", "ckd_flag",
  "diabetes_doctor", "hypertension_doctor", "high_cholesterol_doctor",
  "hcf_a_score", "hcf_a_layer"
)
existing_full <- full_vars[full_vars %in% names(data)]
analysis_data_full <- data[, existing_full]
# 核心版变量列表（16个）- 完全遵循详表
core_vars <- c(
  "SEQN", "SDMVPSU", "SDMVSTRA", "WTMECPRP",
  "target_population", "in_analysis",
  "RIDAGEYR", "RIAGENDR",
  "diabetes_doctor", "hypertension_doctor", "high_cholesterol_doctor",
  "hs_crp_mgl", "high_hs_crp", "white_blood_cells_x10e9_l", "nlr_ratio", "high_nlr",
  "egfr_ckdepi_2021", "ckd_flag",
  "alt_ui_l", "elevated_alt_traditional"
)
existing_core <- core_vars[core_vars %in% names(data)]
analysis_data_core <- data[, existing_core]
# 保存文件 - 使用 _P 后缀
saveRDS(analysis_data_full, file.path(CLEAN_DATA_DIR, "02_biomarker_full_P.rds"))
saveRDS(analysis_data_core, file.path(CLEAN_DATA_DIR, "02_biomarker_core_P.rds"))
cat(sprintf(" ✅ 完整版已保存: 02_biomarker_full_P.rds (%d个变量)\n", ncol(analysis_data_full)))
cat(sprintf(" ✅ 核心版已保存: 02_biomarker_core_P.rds (%d个变量)\n", ncol(analysis_data_core)))
# 生成版本说明
version_note <- data.frame(
  输出文件 = c("02_biomarker_full_P.rds", "02_biomarker_core_P.rds"),
  变量数量 = c(ncol(analysis_data_full), ncol(analysis_data_core)),
  构建日期 = format(Sys.Date(), "%Y-%m-%d"),
  适用场景 = c("敏感性分析/扩展研究/后备论文", "三篇论文主分析"),
  强制规范 = c("❌ 非主分析使用", "✅ 主分析必须使用")
)
write.csv(version_note, file.path(CLEAN_DATA_DIR, "02_biomarker_version_control_P.csv"), row.names = FALSE)
cat(" ✅ 版本说明已保存: 02_biomarker_version_control_P.csv\n\n")
# ============================================================================
# 6. 数据质量评估
# ============================================================================
cat("5. 数据质量评估...\n")
biomarker_vars <- c("ast_ui_l", "alt_ui_l", "egfr_ckdepi_2021",
                    "hs_crp_mgl", "nlr_ratio", "diabetes_doctor")
missing_analysis <- data.frame(
  Variable = biomarker_vars,
  N_Total = nrow(data),
  N_Missing = sapply(biomarker_vars, function(v) {
    if(v %in% names(data)) sum(is.na(data[[v]])) else NA
  })
)
missing_analysis$Pct_Missing <- round(missing_analysis$N_Missing / missing_analysis$N_Total * 100, 1)
cat("缺失值分析：\n")
print(missing_analysis)
high_missing <- missing_analysis[missing_analysis$Pct_Missing > 10 & !is.na(missing_analysis$Pct_Missing), ]
if(nrow(high_missing) > 0) {
  cat(" ⚠️ 高缺失率变量（>10%）：", paste(high_missing$Variable, collapse = ", "), "\n")
} else {
  cat(" ✅ 所有生物标志物变量缺失率≤10%\n")
}
cat("\n")
# ============================================================================
# 7. 生成方法学描述文件
# ============================================================================
cat("6. 生成投稿材料...\n")
methods_file <- file.path(CLEAN_DATA_DIR, "methods_description_nhanes_P.txt")
sink(methods_file)
cat("NHANES数据分析方法（P周期 - 2017-2020）\n")
cat("=============================================\n\n")
cat("数据来源\n")
cat("--------\n")
cat("本研究使用美国国家健康与营养调查（NHANES）2017-2020周期数据。\n\n")
cat("样本选择\n")
cat("--------\n")
cat(sprintf("分析限于成人（年龄≥18岁）非孕妇（n=%d）。\n", 
            sum(data$target_population == 1, na.rm = TRUE)))
cat("按照NHANES教程要求，未删除任何观测，而是创建亚组指示变量用于分析。\n\n")
cat("变量定义\n")
cat("--------\n")
cat("1. 肝功能：丙氨酸氨基转移酶（ALT）、天冬氨酸氨基转移酶（AST）\n")
cat("2. 肾功能：使用CKD-EPI 2021公式计算估算肾小球滤过率（eGFR）\n")
cat("3. 炎症标志物：高敏感C反应蛋白（hsCRP）、中性粒细胞-淋巴细胞比值（NLR）\n")
cat("4. 代谢状态：通过医生诊断变量（糖尿病、高血压、高胆固醇）\n\n")
cat("统计分析\n")
cat("--------\n")
cat("所有分析均使用R和survey包进行。复杂调查设计通过Taylor级数线性化方法处理。\n")
sink()
cat(" ✅ 方法学描述已保存\n\n")
# ============================================================================
# 8. 保存会话信息（期刊要求）
# ============================================================================
cat("7. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "04_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期生物标志物构建会话信息\n")
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
# 9. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n8. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "04_variable_biomarker_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "04_code_list_P.txt")
cat("脚本名称: 04_variable_biomarker_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 10. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期生物标志物构建完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", CLEAN_DATA_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "LOG_DIR")))
gc()
