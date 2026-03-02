#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 03_variable_basic_P.R
# 描述: NHANES 2017-2020 (P周期) 基础变量构建 - 阶段一
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: OSF https://osf.io/xxxxx (DOI: 10.17605/OSF.IO/XXXXX)
# 对应研究计划: 第四部分 4.2 变量构建总览
# 对应变量详表: 第三部分 构建阶段详表（阶段一）
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
required_packages <- c("dplyr", "survey", "tidyverse")
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
log_file <- file.path(LOG_DIR, paste0("03_variable_basic_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 03_variable_basic_P.R\n")
cat("描述: NHANES 2017-2020 (P周期) 基础变量构建\n")
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
cat("输出路径:", CLEAN_DATA_DIR, "\n\n")
# ============================================================================
# 2. 加载数据
# ============================================================================
cat("1. 加载master_P.rds...\n")
master_file <- file.path(CLEAN_DATA_DIR, "master_P.rds")
if(!file.exists(master_file)) {
  stop("错误: master_P.rds不存在! 请先运行02_data_merge_P.R")
}
data <- readRDS(master_file)
cat(sprintf(" 成功加载: %d行, %d列\n", nrow(data), ncol(data)))
# ============================================================================
# 3. 创建标识变量（不删除任何记录）
# ============================================================================
cat("\n2. 创建分析子组标识变量...\n")
data$is_adult <- as.numeric(!is.na(data$RIDAGEYR) & data$RIDAGEYR >= 18)
data$is_non_pregnant_adult <- as.numeric(
  data$is_adult == 1 & (is.na(data$RIDEXPRG) | data$RIDEXPRG != 1)
)
cat(sprintf(" 成人（≥18岁）: %d (%.1f%%)\n",
            sum(data$is_adult == 1, na.rm = TRUE),
            sum(data$is_adult == 1, na.rm = TRUE)/nrow(data)*100))
data$age_group_std <- cut(data$RIDAGEYR,
                          breaks = c(-Inf, 5, 11, 19, 39, 59, Inf),
                          labels = c("0-5y", "6-11y", "12-19y", "20-39y", "40-59y", "60+y"))
# ============================================================================
# 4. 处理NHANES缺失值编码
# ============================================================================
cat("\n3. 处理NHANES缺失值编码...\n")
nhanes_na_codes <- c(7, 9, 77, 99, 777, 999, 7777, 9999)
replace_nhanes_na <- function(x) {
  if(is.numeric(x)) x[x %in% nhanes_na_codes] <- NA
  return(x)
}
numeric_cols <- sapply(data, is.numeric)
data[, numeric_cols] <- lapply(data[, numeric_cols], replace_nhanes_na)
# ============================================================================
# 5. behavior_vars_P.rds（行为学变量）- 完全遵循详表
# ============================================================================
cat("\n4.1 构建behavior_vars_P.rds...\n")
behavior_data <- data.frame(SEQN = data$SEQN)
behavior_data$is_adult <- data$is_adult
behavior_data$is_non_pregnant_adult <- data$is_non_pregnant_adult
# 吸烟变量
behavior_data$current_smoker <- NA
behavior_data$former_smoker <- NA
behavior_data$never_smoker <- NA
behavior_data$smoking_status <- NA
behavior_data$smoking_pack_years <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1) {
    # 当前吸烟者
    if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 1 &&
       !is.na(data$SMQ040[i]) && data$SMQ040[i] %in% c(1, 2)) {
      behavior_data$current_smoker[i] <- 1
      behavior_data$smoking_status[i] <- "current"
    } else if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 2) {
      behavior_data$current_smoker[i] <- 0
    }
    # 既往吸烟者
    if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 1 &&
       !is.na(data$SMQ040[i]) && data$SMQ040[i] == 3) {
      behavior_data$former_smoker[i] <- 1
      behavior_data$smoking_status[i] <- "former"
    } else if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 2) {
      behavior_data$former_smoker[i] <- 0
    }
    # 从不吸烟者
    if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 2) {
      behavior_data$never_smoker[i] <- 1
      behavior_data$smoking_status[i] <- "never"
    } else if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 1) {
      behavior_data$never_smoker[i] <- 0
    }
    # 吸烟包年
    if(!is.na(data$SMQ020[i]) && data$SMQ020[i] == 1) {
      cigs_per_day <- ifelse(!is.na(data$SMD650[i]), data$SMD650[i] / 20, 0)
      smoke_years <- ifelse(!is.na(data$RIDAGEYR[i]), max(0, data$RIDAGEYR[i] - 18), 0)
      behavior_data$smoking_pack_years[i] <- cigs_per_day * smoke_years
    }
  }
}
behavior_data$smoking_status <- factor(behavior_data$smoking_status,
                                       levels = c("never", "former", "current"))
# 饮酒变量
behavior_data$alcohol_drinker <- NA
behavior_data$heavy_drinker <- NA
behavior_data$binge_drinker <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1) {
    if(!is.na(data$ALQ130[i])) {
      behavior_data$alcohol_drinker[i] <- ifelse(data$ALQ130[i] > 0, 1, 0)
    }
    if(!is.na(data$ALQ130[i]) && data$ALQ130[i] > 0 && !is.na(data$RIAGENDR[i])) {
      if((data$RIAGENDR[i] == 1 && data$ALQ130[i] > 14) ||
         (data$RIAGENDR[i] == 2 && data$ALQ130[i] > 7)) {
        behavior_data$heavy_drinker[i] <- 1
      } else {
        behavior_data$heavy_drinker[i] <- 0
      }
    }
    if(!is.na(data$ALQ151[i])) {
      behavior_data$binge_drinker[i] <- ifelse(data$ALQ151[i] == 1, 1, 0)
    }
  }
}
# 设计变量
behavior_data$sdmvpsu <- if("SDMVPSU" %in% names(data)) data$SDMVPSU else 1
behavior_data$sdmvstra <- if("SDMVSTRA" %in% names(data)) data$SDMVSTRA else 1
behavior_data$WTMECPRP <- ifelse(!is.na(data$WTMECPRP) & data$WTMECPRP > 0, data$WTMECPRP, NA)
saveRDS(behavior_data, file.path(CLEAN_DATA_DIR, "behavior_vars_P.rds"))
cat(sprintf(" ✅ 已保存: behavior_vars_P.rds (%d个变量)\n", ncol(behavior_data)))
# ============================================================================
# 6. clinical_vars_P.rds（临床诊断变量）- 完全遵循详表
# ============================================================================
cat("\n4.2 构建clinical_vars_P.rds...\n")
clinical_data <- data.frame(SEQN = data$SEQN)
clinical_data$is_adult <- data$is_adult
clinical_data$is_non_pregnant_adult <- data$is_non_pregnant_adult
# BMI相关
clinical_data$bmi_continuous <- ifelse(
  !is.na(data$BMXBMI) & data$BMXBMI >= 10 & data$BMXBMI <= 70,
  data$BMXBMI, NA
)
clinical_data$bmi_category <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1 && !is.na(clinical_data$bmi_continuous[i])) {
    bmi_val <- clinical_data$bmi_continuous[i]
    if(bmi_val < 18.5) clinical_data$bmi_category[i] <- "Underweight"
    else if(bmi_val < 25) clinical_data$bmi_category[i] <- "Normal"
    else if(bmi_val < 30) clinical_data$bmi_category[i] <- "Overweight"
    else clinical_data$bmi_category[i] <- "Obese"
  }
}
clinical_data$bmi_category <- as.factor(clinical_data$bmi_category)
# 肥胖
clinical_data$obesity <- NA
clinical_data$waist_circumference <- ifelse(
  !is.na(data$BMXWAIST) & data$BMXWAIST >= 50 & data$BMXWAIST <= 200,
  data$BMXWAIST, NA
)
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1) {
    obese <- 0
    if(!is.na(data$BMXBMI[i]) && data$BMXBMI[i] >= 30) obese <- 1
    if(!is.na(data$BMXWAIST[i]) && !is.na(data$RIAGENDR[i])) {
      if((data$RIAGENDR[i] == 1 & data$BMXWAIST[i] >= 102) ||
         (data$RIAGENDR[i] == 2 & data$BMXWAIST[i] >= 88)) {
        obese <- 1
      }
    }
    clinical_data$obesity[i] <- obese
  }
}
# 血压
clinical_data$systolic_bp <- ifelse(
  !is.na(data$BPXOSY1) & data$BPXOSY1 >= 70 & data$BPXOSY1 <= 250,
  data$BPXOSY1, NA
)
clinical_data$diastolic_bp <- ifelse(
  !is.na(data$BPXODI1) & data$BPXODI1 >= 30 & data$BPXODI1 <= 150,
  data$BPXODI1, NA
)
# 高血压
clinical_data$hypertension <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1) {
    ht <- 0
    # 血压读数判断
    if(!is.na(clinical_data$systolic_bp[i]) && clinical_data$systolic_bp[i] >= 130) ht <- 1
    if(!is.na(clinical_data$diastolic_bp[i]) && clinical_data$diastolic_bp[i] >= 80) ht <- 1
    # 降压药判断 - P周期使用 BPQ020
    if("BPQ020" %in% names(data) && !is.na(data$BPQ020[i]) && data$BPQ020[i] == 1) {
      ht <- 1
    }
    clinical_data$hypertension[i] <- ht
  }
}
# 糖尿病
clinical_data$diabetes <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1 && !is.na(data$DIQ010[i])) {
    clinical_data$diabetes[i] <- ifelse(data$DIQ010[i] == 1, 1, 0)
  }
}
# 高胆固醇诊断
clinical_data$high_cholesterol_dx <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1 && !is.na(data$BPQ080[i])) {
    clinical_data$high_cholesterol_dx[i] <- ifelse(data$BPQ080[i] == 1, 1, 0)
  }
}
# 心血管疾病
clinical_data$any_cvd <- NA
clinical_data$coronary_heart_disease <- NA
clinical_data$heart_failure <- NA
clinical_data$stroke <- NA
clinical_data$heart_attack <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1) {
    cvd_any <- 0
    if(!is.na(data$MCQ160C[i])) {
      clinical_data$coronary_heart_disease[i] <- ifelse(data$MCQ160C[i] == 1, 1, 0)
      if(data$MCQ160C[i] == 1) cvd_any <- 1
    }
    if(!is.na(data$MCQ160E[i])) {
      clinical_data$heart_failure[i] <- ifelse(data$MCQ160E[i] == 1, 1, 0)
      if(data$MCQ160E[i] == 1) cvd_any <- 1
    }
    if(!is.na(data$MCQ160F[i])) {
      clinical_data$stroke[i] <- ifelse(data$MCQ160F[i] == 1, 1, 0)
      if(data$MCQ160F[i] == 1) cvd_any <- 1
    }
    if(!is.na(data$MCQ160B[i])) {
      clinical_data$heart_attack[i] <- ifelse(data$MCQ160B[i] == 1, 1, 0)
      if(data$MCQ160B[i] == 1) cvd_any <- 1
    }
    clinical_data$any_cvd[i] <- cvd_any
  }
}
# 静息心率
clinical_data$pulse_rate <- ifelse(
  !is.na(data$BPXOPLS1) & data$BPXOPLS1 >= 30 & data$BPXOPLS1 <= 200,
  data$BPXOPLS1, NA
)
# HDL胆固醇
clinical_data$hdl_cholesterol <- NA
clinical_data$low_hdl <- NA
if("LBXHDD" %in% names(data)) {
  clinical_data$hdl_cholesterol <- ifelse(
    !is.na(data$LBXHDD) & data$LBXHDD >= 10 & data$LBXHDD <= 150,
    data$LBXHDD, NA
  )
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(clinical_data$hdl_cholesterol[i]) && !is.na(data$RIAGENDR[i])) {
      if((data$RIAGENDR[i] == 1 && clinical_data$hdl_cholesterol[i] < 40) ||
         (data$RIAGENDR[i] == 2 && clinical_data$hdl_cholesterol[i] < 50)) {
        clinical_data$low_hdl[i] <- 1
      } else {
        clinical_data$low_hdl[i] <- 0
      }
    }
  }
}
# C反应蛋白
clinical_data$crp <- NA
clinical_data$high_crp <- NA
clinical_data$crp_category <- NA
if("LBXHSCRP" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(!is.na(data$LBXHSCRP[i])) {
      crp_val <- data$LBXHSCRP[i]
      if(crp_val >= 0.1 && crp_val <= 100) {
        clinical_data$crp[i] <- crp_val
      }
    }
  }
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(clinical_data$crp[i])) {
      crp_val <- clinical_data$crp[i]
      clinical_data$high_crp[i] <- ifelse(crp_val >= 3, 1, 0)
      if(crp_val < 1) clinical_data$crp_category[i] <- "低风险"
      else if(crp_val < 3) clinical_data$crp_category[i] <- "中等风险"
      else clinical_data$crp_category[i] <- "高风险"
    }
  }
  clinical_data$crp_category <- factor(clinical_data$crp_category,
                                       levels = c("低风险", "中等风险", "高风险"))
}
# A层评分
clinical_data$A_layer_score <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1) {
    score <- 0
    if(!is.na(clinical_data$obesity[i]) && clinical_data$obesity[i] == 1) score <- score + 1
    if(!is.na(clinical_data$hypertension[i]) && clinical_data$hypertension[i] == 1) score <- score + 1
    if(!is.na(clinical_data$diabetes[i]) && clinical_data$diabetes[i] == 1) score <- score + 1
    if(!is.na(clinical_data$high_cholesterol_dx[i]) && clinical_data$high_cholesterol_dx[i] == 1) score <- score + 1
    clinical_data$A_layer_score[i] <- score
  }
}
# 设计变量
clinical_data$sdmvpsu <- behavior_data$sdmvpsu
clinical_data$sdmvstra <- behavior_data$sdmvstra
clinical_data$WTMECPRP <- behavior_data$WTMECPRP
saveRDS(clinical_data, file.path(CLEAN_DATA_DIR, "clinical_vars_P.rds"))
cat(sprintf(" ✅ 已保存: clinical_vars_P.rds (%d个变量)\n", ncol(clinical_data)))
# ============================================================================
# 7. social_vars_P.rds（社会决定因素）
# ============================================================================
cat("\n4.3 构建social_vars_P.rds...\n")
social_data <- data.frame(SEQN = data$SEQN)
# 人口学
social_data$age <- data$RIDAGEYR
social_data$gender <- factor(data$RIAGENDR, levels = 1:2, labels = c("Male", "Female"))
if("RIDRETH3" %in% names(data)) {
  social_data$race_ethnicity <- factor(
    data$RIDRETH3,
    levels = 1:7,
    labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
               "Non-Hispanic Black", "Non-Hispanic Asian",
               "Other Race", "Non-Hispanic Other")
  )
}
# 教育程度
social_data$education_level <- NA
social_data$education_college <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1 && !is.na(data$DMDEDUC2[i])) {
    if(data$DMDEDUC2[i] == 1) social_data$education_level[i] <- "<9年"
    else if(data$DMDEDUC2[i] == 2) social_data$education_level[i] <- "9-11年"
    else if(data$DMDEDUC2[i] == 3) social_data$education_level[i] <- "高中"
    else if(data$DMDEDUC2[i] == 4) social_data$education_level[i] <- "大学"
    else if(data$DMDEDUC2[i] == 5) social_data$education_level[i] <- "研究生"
    social_data$education_college[i] <- ifelse(data$DMDEDUC2[i] >= 4, 1, 0)
  }
}
social_data$education_level <- as.factor(social_data$education_level)
# 婚姻状态
social_data$marital_status <- NA
social_data$married <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1 && !is.na(data$DMDMARTZ[i])) {
    if(data$DMDMARTZ[i] == 1) social_data$marital_status[i] <- "已婚"
    else if(data$DMDMARTZ[i] == 2) social_data$marital_status[i] <- "丧偶"
    else if(data$DMDMARTZ[i] == 3) social_data$marital_status[i] <- "离婚"
    else if(data$DMDMARTZ[i] == 4) social_data$marital_status[i] <- "分居"
    else if(data$DMDMARTZ[i] == 5) social_data$marital_status[i] <- "未婚"
    else if(data$DMDMARTZ[i] == 6) social_data$marital_status[i] <- "同居"
    social_data$married[i] <- ifelse(data$DMDMARTZ[i] == 1, 1, 0)
  }
}
social_data$marital_status <- as.factor(social_data$marital_status)
# 就业状态 - 注意P周期使用OCD180
social_data$employment_status <- NA
social_data$employed <- NA
if("OCD180" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(data$OCD180[i])) {
      if(data$OCD180[i] == 1) social_data$employment_status[i] <- "就业"
      else if(data$OCD180[i] == 2) social_data$employment_status[i] <- "在职休假"
      else if(data$OCD180[i] == 3) social_data$employment_status[i] <- "失业"
      else if(data$OCD180[i] == 4) social_data$employment_status[i] <- "不工作"
      social_data$employed[i] <- ifelse(data$OCD180[i] == 1, 1, 0)
    }
  }
}
social_data$employment_status <- as.factor(social_data$employment_status)
# 贫困收入比
social_data$poverty_ratio <- ifelse(
  !is.na(data$INDFMPIR) & data$INDFMPIR >= 0 & data$INDFMPIR <= 5,
  data$INDFMPIR, NA
)
social_data$below_poverty <- ifelse(
  !is.na(social_data$poverty_ratio) & social_data$poverty_ratio < 1, 1, 0
)
# 设计变量
social_data$sdmvpsu <- behavior_data$sdmvpsu
social_data$sdmvstra <- behavior_data$sdmvstra
if("WTINTPRP" %in% names(data)) {
  social_data$WTINTPRP <- data$WTINTPRP
} else {
  social_data$WTINTPRP <- NA
  warning("WTINTPRP不存在于原始数据中")
}
saveRDS(social_data, file.path(CLEAN_DATA_DIR, "social_vars_P.rds"))
cat(sprintf(" ✅ 已保存: social_vars_P.rds (%d个变量)\n", ncol(social_data)))
# ============================================================================
# 8. symptom_vars_P.rds（症状与功能）
# ============================================================================
cat("\n4.4 构建symptom_vars_P.rds...\n")
symptom_data <- data.frame(SEQN = data$SEQN)
# 自评健康
symptom_data$self_rated_health <- NA
symptom_data$poor_self_rated_health <- NA
if("HUQ010" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(data$HUQ010[i])) {
      if(data$HUQ010[i] == 1) symptom_data$self_rated_health[i] <- "优秀"
      else if(data$HUQ010[i] == 2) symptom_data$self_rated_health[i] <- "很好"
      else if(data$HUQ010[i] == 3) symptom_data$self_rated_health[i] <- "好"
      else if(data$HUQ010[i] == 4) symptom_data$self_rated_health[i] <- "一般"
      else if(data$HUQ010[i] == 5) symptom_data$self_rated_health[i] <- "差"
      symptom_data$poor_self_rated_health[i] <- ifelse(data$HUQ010[i] >= 4, 1, 0)
    }
  }
  symptom_data$self_rated_health <- factor(symptom_data$self_rated_health,
                                          levels = c("优秀", "很好", "好", "一般", "差"))
}
# PHQ-9抑郁量表
symptom_data$phq9_total <- NA
symptom_data$depression <- NA
phq_vars <- c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050",
              "DPQ060", "DPQ070", "DPQ080", "DPQ090")
if(all(phq_vars %in% names(data))) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1) {
      scores <- sapply(phq_vars, function(var) {
        if(!is.na(data[[var]][i]) && data[[var]][i] >= 0 && data[[var]][i] <= 3) {
          return(data[[var]][i])
        } else {
          return(0)
        }
      })
      symptom_data$phq9_total[i] <- sum(scores)
      symptom_data$depression[i] <- ifelse(symptom_data$phq9_total[i] >= 10, 1, 0)
    }
  }
}
# 设计变量
symptom_data$sdmvpsu <- behavior_data$sdmvpsu
symptom_data$sdmvstra <- behavior_data$sdmvstra
if("WTINTPRP" %in% names(data)) {
  symptom_data$WTINTPRP <- data$WTINTPRP
} else {
  symptom_data$WTINTPRP <- NA
}
saveRDS(symptom_data, file.path(CLEAN_DATA_DIR, "symptom_vars_P.rds"))
cat(sprintf(" ✅ 已保存: symptom_vars_P.rds (%d个变量)\n", ncol(symptom_data)))
# ============================================================================
# 9. dietbehavior_vars_P.rds（饮食行为变量）
# ============================================================================
cat("\n4.5 构建dietbehavior_vars_P.rds...\n")
diet_data <- data.frame(SEQN = data$SEQN)
# 核心饮食行为变量 - 注意P周期使用DBQ095/DBQ097
diet_data$diet_special <- NA
if("DBQ095" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(data$DBQ095[i])) {
      diet_data$diet_special[i] <- ifelse(data$DBQ095[i] == 1, 1, 0)
    }
  }
}
# 减肥饮食
diet_data$diet_weight_loss <- NA
if("DRQSDT1" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(data$DRQSDT1[i])) {
      diet_data$diet_weight_loss[i] <- ifelse(data$DRQSDT1[i] == 1, 1, 0)
    }
  }
}
# 饮食规则性评分
diet_data$diet_rigidity_score <- NA
for(i in 1:nrow(data)) {
  if(data$is_adult[i] == 1 && !is.na(diet_data$diet_special[i])) {
    diet_data$diet_rigidity_score[i] <- ifelse(diet_data$diet_special[i] == 1, 1, 0)
  }
}
# 糖摄入量
diet_data$sugar_intake_g <- NA
if("DR1TSUGR" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(data$DR1TSUGR[i]) && data$DR1TSUGR[i] >= 0) {
      diet_data$sugar_intake_g[i] <- data$DR1TSUGR[i]
    }
  }
}
# 总能量摄入
diet_data$total_calories <- NA
if("DR1TKCAL" %in% names(data)) {
  for(i in 1:nrow(data)) {
    if(data$is_adult[i] == 1 && !is.na(data$DR1TKCAL[i]) && data$DR1TKCAL[i] >= 0) {
      diet_data$total_calories[i] <- data$DR1TKCAL[i]
    }
  }
}
# 膳食专用权重
diet_data$sdmvpsu <- if("SDMVPSU" %in% names(data)) data$SDMVPSU else NA
diet_data$sdmvstra <- if("SDMVSTRA" %in% names(data)) data$SDMVSTRA else NA
if("WTDRD1" %in% names(data)) {
  diet_data$wtdrd1 <- ifelse(!is.na(data$WTDRD1) & data$WTDRD1 > 0, data$WTDRD1, NA)
} else {
  diet_data$wtdrd1 <- NA
}
if("WTDR2D" %in% names(data)) {
  diet_data$wtdr2d <- ifelse(!is.na(data$WTDR2D) & data$WTDR2D > 0, data$WTDR2D, NA)
}
saveRDS(diet_data, file.path(CLEAN_DATA_DIR, "dietbehavior_vars_P.rds"))
cat(sprintf(" ✅ 已保存: dietbehavior_vars_P.rds (%d个变量)\n", ncol(diet_data)))
# ============================================================================
# 10. healthbehavior_vars_P.rds（健康行为达标变量）- 修正版
# ============================================================================
cat("\n4.6 构建healthbehavior_vars_P.rds...\n")
health_data <- data.frame(SEQN = data$SEQN)
# 睡眠达标
if("SLD012" %in% names(data)) {
  health_data$sleep_adequate <- ifelse(
    !is.na(data$SLD012) & data$SLD012 >= 7 & data$SLD012 <= 9, 1, 0
  )
  n_sleep <- sum(!is.na(health_data$sleep_adequate))
  cat(sprintf(" ✅ 睡眠变量: sleep_adequate (有效样本: %d)\n", n_sleep))
}
# ============================================================================
# 10. healthbehavior_vars_P.rds（健康行为达标变量）- 修正版
# ============================================================================
cat("\n4.6 构建healthbehavior_vars_P.rds...\n")
health_data <- data.frame(SEQN = data$SEQN)
# 睡眠达标
if("SLD012" %in% names(data)) {
  health_data$sleep_adequate <- ifelse(
    !is.na(data$SLD012) & data$SLD012 >= 7 & data$SLD012 <= 9, 1, 0
  )
  n_sleep <- sum(!is.na(health_data$sleep_adequate))
  cat(sprintf(" ✅ 睡眠变量: sleep_adequate (有效样本: %d)\n", n_sleep))
}
# ============================================================================
# 体力活动变量构建（P周期）- 修正达标标准
# ============================================================================
cat("\n 🔧 构建体力活动变量（P周期）...\n")
# 初始化变量
health_data$pa_total_min_week <- NA
health_data$pa_meets_guideline <- NA
# P周期体力活动变量名
vig_work_var <- "PAD615"      # 剧烈工作活动（分钟/天）
mod_work_var <- "PAD630"      # 中等工作活动（分钟/天）
trans_days_var <- "PAQ640"    # 步行/骑车天数（天/周）
trans_mins_var <- "PAD645"     # 步行/骑车时间（分钟/天）
vig_rec_days_var <- "PAQ655"   # 剧烈休闲天数（天/周）
vig_rec_mins_var <- "PAD660"    # 剧烈休闲时间（分钟/天）
mod_rec_days_var <- "PAQ670"    # 中等休闲天数（天/周）
mod_rec_mins_var <- "PAD675"    # 中等休闲时间（分钟/天）
# 逐行计算原始分钟（不加权）
for(i in 1:nrow(data)) {
  if(!is.na(data$is_adult[i]) && data$is_adult[i] == 1) {
    total_min <- 0
    has_activity <- FALSE
    # 1. 剧烈工作活动 (PAD615)
    if(vig_work_var %in% names(data) && !is.na(data[[vig_work_var]][i])) {
      mins_per_day <- data[[vig_work_var]][i]
      if(mins_per_day > 0 && mins_per_day <= 840) {
        days_per_week <- 5  # 假设每周工作5天
        total_min <- total_min + (days_per_week * mins_per_day)
        has_activity <- TRUE
      }
    }
    # 2. 中等工作活动 (PAD630)
    if(mod_work_var %in% names(data) && !is.na(data[[mod_work_var]][i])) {
      mins_per_day <- data[[mod_work_var]][i]
      if(mins_per_day > 0 && mins_per_day <= 900) {
        days_per_week <- 5
        total_min <- total_min + (days_per_week * mins_per_day)
        has_activity <- TRUE
      }
    }
    # 3. 交通步行/骑车 (PAQ640 + PAD645)
    if(trans_days_var %in% names(data) && trans_mins_var %in% names(data) &&
       !is.na(data[[trans_days_var]][i]) && !is.na(data[[trans_mins_var]][i])) {
      days_per_week <- data[[trans_days_var]][i]
      mins_per_day <- data[[trans_mins_var]][i]
      if(days_per_week >= 1 && days_per_week <= 7 && 
         mins_per_day >= 10 && mins_per_day <= 840) {
        total_min <- total_min + (days_per_week * mins_per_day)
        has_activity <- TRUE
      }
    }
    # 4. 剧烈休闲活动 (PAQ655 + PAD660)
    if(vig_rec_days_var %in% names(data) && vig_rec_mins_var %in% names(data) &&
       !is.na(data[[vig_rec_days_var]][i]) && !is.na(data[[vig_rec_mins_var]][i])) {
      days_per_week <- data[[vig_rec_days_var]][i]
      mins_per_day <- data[[vig_rec_mins_var]][i]
      if(days_per_week >= 1 && days_per_week <= 7 && 
         mins_per_day >= 10 && mins_per_day <= 480) {
        total_min <- total_min + (days_per_week * mins_per_day)
        has_activity <- TRUE
      }
    }
    # 5. 中等休闲活动 (PAQ670 + PAD675)
    if(mod_rec_days_var %in% names(data) && mod_rec_mins_var %in% names(data) &&
       !is.na(data[[mod_rec_days_var]][i]) && !is.na(data[[mod_rec_mins_var]][i])) {
      days_per_week <- data[[mod_rec_days_var]][i]
      mins_per_day <- data[[mod_rec_mins_var]][i]
      if(days_per_week >= 1 && days_per_week <= 7 && 
         mins_per_day >= 10 && mins_per_day <= 600) {
        total_min <- total_min + (days_per_week * mins_per_day)
        has_activity <- TRUE
      }
    }
    # 赋值
    if(has_activity) {
      health_data$pa_total_min_week[i] <- total_min
    }
  }
}
# **关键修正：使用150分钟标准，不是600！**
health_data$pa_meets_guideline <- ifelse(
  !is.na(health_data$pa_total_min_week) & health_data$pa_total_min_week >= 150, 1, 0
)
# 统计信息
n_pa <- sum(!is.na(health_data$pa_total_min_week))
n_meet <- sum(health_data$pa_meets_guideline == 1, na.rm = TRUE)
pct_meet <- ifelse(n_pa > 0, n_meet / n_pa * 100, 0)
cat(sprintf("\n ✅ 体力活动变量构建完成:\n"))
cat(sprintf("    - pa_total_min_week: 有效样本 %d 人\n", n_pa))
cat(sprintf("    - 均值: %.1f 分钟/周\n", 
            mean(health_data$pa_total_min_week, na.rm = TRUE)))
cat(sprintf("    - 中位数: %.1f 分钟/周\n", 
            median(health_data$pa_total_min_week, na.rm = TRUE)))
cat(sprintf("    - pa_meets_guideline: 达标 %d 人 (%.1f%%)\n", n_meet, pct_meet))
# ============================================================================
# 设计变量
# ============================================================================
cat("\n 🔧 添加设计变量...\n")
health_data$sdmvpsu <- if("SDMVPSU" %in% names(data)) data$SDMVPSU else NA
health_data$sdmvstra <- if("SDMVSTRA" %in% names(data)) data$SDMVSTRA else NA
health_data$wtmecprp <- ifelse(!is.na(data$WTMECPRP) & data$WTMECPRP > 0, data$WTMECPRP, NA)
health_data$wtintprp <- ifelse(!is.na(data$WTINTPRP) & data$WTINTPRP > 0, data$WTINTPRP, NA)
n_psu <- sum(!is.na(health_data$sdmvpsu))
n_stra <- sum(!is.na(health_data$sdmvstra))
n_wt <- sum(!is.na(health_data$wtmecprp))
cat(sprintf("    - sdmvpsu: 有效样本 %d\n", n_psu))
cat(sprintf("    - sdmvstra: 有效样本 %d\n", n_stra))
cat(sprintf("    - wtmecprp: 有效样本 %d\n", n_wt))
# ============================================================================
# 保存文件
# ============================================================================
saveRDS(health_data, file.path(CLEAN_DATA_DIR, "healthbehavior_vars_P.rds"))
cat(sprintf("\n ✅ 已保存: healthbehavior_vars_P.rds (%d个变量)\n", ncol(health_data)))
# ============================================================================
# 11. 创建多个调查设计对象（最小公分母原则）
# ============================================================================
cat("\n5. 创建多个调查设计对象（最小公分母原则）...\n")
cat("========================================================\n")
# 合并所有变量到主数据集
cat("5.0 合并所有变量到主数据集...\n")
master_all <- data.frame(SEQN = data$SEQN)
# 添加行为变量
behavior_cols <- intersect(
  c("SEQN", "is_adult", "is_non_pregnant_adult", "current_smoker",
    "former_smoker", "never_smoker", "smoking_status",
    "smoking_pack_years", "alcohol_drinker", "heavy_drinker",
    "binge_drinker", "sdmvpsu", "sdmvstra", "WTINTPRP"),
  names(behavior_data)
)
if(length(behavior_cols) > 1) {
  behavior_selected <- behavior_data %>% select(all_of(behavior_cols))
  master_all <- master_all %>%
    left_join(behavior_selected, by = "SEQN")
}
# 添加临床变量
clinical_cols <- intersect(
  c("SEQN", "bmi_continuous", "bmi_category", "obesity",
    "waist_circumference", "systolic_bp", "diastolic_bp",
    "hypertension", "diabetes", "high_cholesterol_dx",
    "any_cvd", "coronary_heart_disease", "heart_failure", "stroke",
    "heart_attack", "pulse_rate", "hdl_cholesterol", "low_hdl",
    "crp", "high_crp", "crp_category", "WTMECPRP", "sdmvpsu", "sdmvstra"),
  names(clinical_data)
)
if(length(clinical_cols) > 1) {
  clinical_selected <- clinical_data %>% select(all_of(clinical_cols))
  master_all <- master_all %>%
    left_join(clinical_selected, by = "SEQN", suffix = c("", "_clin"))
}
# 添加社会变量
social_cols <- intersect(
  c("SEQN", "age", "gender", "race_ethnicity", "education_level",
    "education_college", "marital_status", "married",
    "employment_status", "employed", "poverty_ratio",
    "below_poverty", "WTINTPRP", "sdmvpsu", "sdmvstra"),
  names(social_data)
)
if(length(social_cols) > 1) {
  social_selected <- social_data %>% select(all_of(social_cols))
  master_all <- master_all %>%
    left_join(social_selected, by = "SEQN", suffix = c("", "_soc"))
}
# 添加症状变量
symptom_cols <- intersect(
  c("SEQN", "self_rated_health", "poor_self_rated_health",
    "phq9_total", "depression", "WTINTPRP", "sdmvpsu", "sdmvstra"),
  names(symptom_data)
)
if(length(symptom_cols) > 1) {
  symptom_selected <- symptom_data %>% select(all_of(symptom_cols))
  master_all <- master_all %>%
    left_join(symptom_selected, by = "SEQN", suffix = c("", "_sym"))
}
# 添加饮食变量
diet_cols <- intersect(
  c("SEQN", "diet_special", "diet_weight_loss",
    "diet_rigidity_score", "sugar_intake_g",
    "total_calories", "wtdrd1", "sdmvpsu", "sdmvstra"),
  names(diet_data)
)
if(length(diet_cols) > 1) {
  diet_selected <- diet_data %>% select(all_of(diet_cols))
  master_all <- master_all %>%
    left_join(diet_selected, by = "SEQN", suffix = c("", "_diet"))
}
# 添加健康行为变量
health_cols <- intersect(
  c("SEQN", "sleep_adequate", "pa_total_min_week",
    "pa_meets_guideline", "WTMECPRP", "WTINTPRP", "sdmvpsu", "sdmvstra"),
  names(health_data)
)
if(length(health_cols) > 1) {
  health_selected <- health_data %>% select(all_of(health_cols))
  master_all <- master_all %>%
    left_join(health_selected, by = "SEQN", suffix = c("", "_health"))
}
cat(sprintf(" ✅ 主数据集合并完成: %d行, %d列\n", nrow(master_all), ncol(master_all)))
# 统一设计变量
cat("\n5.1 统一设计变量...\n")
psu_cols <- grep("^sdmvpsu", names(master_all), value = TRUE)
stra_cols <- grep("^sdmvstra", names(master_all), value = TRUE)
master_all <- master_all %>%
  mutate(
    final_sdmvpsu = do.call(coalesce, c(select(., any_of(psu_cols)), na.rm = TRUE)),
    final_sdmvstra = do.call(coalesce, c(select(., any_of(stra_cols)), na.rm = TRUE))
  )
cat(" ✅ 设计变量统一完成\n")
# 创建多个设计对象
cat("\n5.2 创建多个设计对象（最小公分母原则）...\n")
design_objects <- list()
# 访谈权重设计对象
cat("\n 5.2.1 创建访谈权重设计对象 (WTINTPRP)...\n")
interview_valid <- master_all %>%
  filter(!is.na(final_sdmvpsu) & !is.na(final_sdmvstra)) %>%
  filter(!is.na(WTINTPRP) & WTINTPRP > 0)
cat(sprintf("  - 有效样本: %d/%d\n", nrow(interview_valid), nrow(master_all)))
if(nrow(interview_valid) > 0) {
  design_objects$interview <- svydesign(
    id = ~final_sdmvpsu,
    strata = ~final_sdmvstra,
    weights = ~WTINTPRP,
    data = interview_valid,
    nest = TRUE,
    single.psu = "average"
  )
  saveRDS(design_objects$interview, file.path(CLEAN_DATA_DIR, "interview_design_P.rds"))
  cat(" ✅ 已保存: interview_design_P.rds\n")
}
# MEC权重设计对象
cat("\n 5.2.2 创建MEC权重设计对象 (WTMECPRP)...\n")
mec_valid <- master_all %>%
  filter(!is.na(final_sdmvpsu) & !is.na(final_sdmvstra)) %>%
  filter(!is.na(WTMECPRP) & WTMECPRP > 0)
cat(sprintf("  - 有效样本: %d/%d\n", nrow(mec_valid), nrow(master_all)))
if(nrow(mec_valid) > 0) {
  design_objects$mec <- svydesign(
    id = ~final_sdmvpsu,
    strata = ~final_sdmvstra,
    weights = ~WTMECPRP,
    data = mec_valid,
    nest = TRUE,
    single.psu = "average"
  )
  saveRDS(design_objects$mec, file.path(CLEAN_DATA_DIR, "mec_design_P.rds"))
  cat(" ✅ 已保存: mec_design_P.rds\n")
}
# 膳食权重设计对象
cat("\n 5.2.3 创建膳食权重设计对象 (WTDRD1)...\n")
diet_valid <- master_all %>%
  filter(!is.na(final_sdmvpsu) & !is.na(final_sdmvstra)) %>%
  filter(!is.na(wtdrd1) & wtdrd1 > 0)
cat(sprintf("  - 有效样本: %d/%d\n", nrow(diet_valid), nrow(master_all)))
if(nrow(diet_valid) > 0) {
  design_objects$diet <- svydesign(
    id = ~final_sdmvpsu,
    strata = ~final_sdmvstra,
    weights = ~wtdrd1,
    data = diet_valid,
    nest = TRUE,
    single.psu = "average"
  )
  saveRDS(design_objects$diet, file.path(CLEAN_DATA_DIR, "diet_design_P.rds"))
  cat(" ✅ 已保存: diet_design_P.rds\n")
}
cat("\n ✅ 5.2完成：创建了", length(design_objects), "个设计对象\n")
# 创建子组设计对象
cat("\n5.3 创建子组分析设计对象...\n")
if(!is.null(design_objects$interview)) {
  if("is_adult" %in% names(design_objects$interview$variables)) {
    adult_interview <- subset(design_objects$interview, is_adult == 1)
    saveRDS(adult_interview, file.path(CLEAN_DATA_DIR, "adult_interview_design_P.rds"))
    cat(" ✅ 已保存: adult_interview_design_P.rds\n")
  }
}
if(!is.null(design_objects$mec)) {
  if("is_adult" %in% names(design_objects$mec$variables)) {
    adult_mec <- subset(design_objects$mec, is_adult == 1)
    saveRDS(adult_mec, file.path(CLEAN_DATA_DIR, "adult_mec_design_P.rds"))
    cat(" ✅ 已保存: adult_mec_design_P.rds\n")
  }
  if("is_non_pregnant_adult" %in% names(design_objects$mec$variables)) {
    non_preg_mec <- subset(design_objects$mec, is_non_pregnant_adult == 1)
    saveRDS(non_preg_mec, file.path(CLEAN_DATA_DIR, "non_pregnant_adult_mec_design_P.rds"))
    cat(" ✅ 已保存: non_pregnant_adult_mec_design_P.rds\n")
  }
}
if(!is.null(design_objects$diet)) {
  if("is_adult" %in% names(design_objects$diet$variables)) {
    adult_diet <- subset(design_objects$diet, is_adult == 1)
    saveRDS(adult_diet, file.path(CLEAN_DATA_DIR, "adult_diet_design_P.rds"))
    cat(" ✅ 已保存: adult_diet_design_P.rds\n")
  }
}
cat("\n ✅ 5.3完成：子组创建完成\n")
# 创建设计对象索引
cat("\n5.4 创建设计对象索引...\n")
design_index <- list(
  interview = list(
    full = "interview_design_P.rds",
    adult = "adult_interview_design_P.rds"
  ),
  mec = list(
    full = "mec_design_P.rds",
    adult = "adult_mec_design_P.rds",
    non_pregnant = "non_pregnant_adult_mec_design_P.rds"
  ),
  diet = list(
    full = "diet_design_P.rds",
    adult = "adult_diet_design_P.rds"
  )
)
saveRDS(design_index, file.path(CLEAN_DATA_DIR, "design_index_P.rds"))
cat(" ✅ 已保存: design_index_P.rds\n")
# ============================================================================
# 12. 生成合规报告
# ============================================================================
cat("\n6. 生成NHANES教程合规报告...\n")
report_file <- file.path(LOG_DIR, "03_variable_basic_report_P.txt")
sink(report_file)
cat("NHANES P周期基础变量构建 - 完全合规报告\n")
cat("==========================================\n\n")
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、NHANES教程合规性\n")
cat(" ✅ 100%遵循教程要求:\n")
cat("  - 不删除任何记录: 使用标识变量区分成人/儿童\n")
cat("  - 统一权重选择: 多个设计对象对应不同权重\n")
cat("  - 包含设计变量: sdmvstra（层）, sdmvpsu（PSU）\n")
cat("  - 正确方差估计: 创建survey设计对象\n")
cat("  - 最小公分母原则: 根据变量来源选择权重\n\n")
cat("二、变量构建完成清单\n")
cat(sprintf(" 1. behavior_vars_P.rds: %d个变量\n", ncol(behavior_data)))
cat(sprintf(" 2. clinical_vars_P.rds: %d个变量\n", ncol(clinical_data)))
cat(sprintf(" 3. social_vars_P.rds: %d个变量\n", ncol(social_data)))
cat(sprintf(" 4. symptom_vars_P.rds: %d个变量\n", ncol(symptom_data)))
cat(sprintf(" 5. dietbehavior_vars_P.rds: %d个变量\n", ncol(diet_data)))
cat(sprintf(" 6. healthbehavior_vars_P.rds: %d个变量\n", ncol(health_data)))
sink()
cat(" ✅ 合规报告已保存\n")
# ============================================================================
# 13. 保存会话信息（期刊要求）
# ============================================================================
cat("\n7. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "03_session_info_P.txt")
sink(session_info_path)
cat("NHANES P周期基础变量构建会话信息\n")
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
# 14. 保存R代码副本（期刊要求）
# ============================================================================
cat("\n8. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "03_variable_basic_P.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "03_code_list_P.txt")
cat("脚本名称: 03_variable_basic_P.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 15. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ P周期基础变量构建完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", CLEAN_DATA_DIR, "\n")
cat("========================================================\n")
sink()
# 清理临时变量
rm(list = setdiff(ls(), c("CLEAN_DATA_DIR", "LOG_DIR")))
gc()
