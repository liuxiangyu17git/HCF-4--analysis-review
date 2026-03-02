#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 16_paper3_alpha_analysis_L.R
# 描述: NHANES 2021-2023 (L周期) 论文3 - α因子正式分析
# 整合版本：包含α1-α4与多系统疾病关联、中介网络、系统特异性分析
# 完全避免循环论证，α3不参与代谢系统预测
# 
# 研究项目: HCF-4α框架：心理生理健康的三层次整合模型
# 对应预注册: 研究问题 Q3, 假设 H6-H10
# 对应研究计划: 第五部分 图3、表3、图4、表4
# 对应变量详表: 第三部分 构建阶段详表（阶段三 - α因子）
#
# 伦理声明: 使用NHANES公开数据，无需IRB批准
# 数据来源: https://www.cdc.gov/nchs/nhanes/index.htm
# 数据周期: 2021-2023 (L系列)
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
required_packages <- c(
  "tidyverse", "survey", "lavaan", "semTools",
  "psych", "ggplot2", "pheatmap", "gridExtra",
  "interactions", "emmeans", "jtools", "corrplot"
)
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
# ============================================================================
# 英文标签定义（只加这一部分，不改任何逻辑）
# ============================================================================
# α因子英文标签
alpha_labels <- c(
  "alpha1" = "α₁ Metacognitive",
  "alpha2" = "α₂ Emotional Regulation",
  "alpha3" = "α₃ Systemic Coordination",
  "alpha4" = "α₄ Goal Efficacy"
)
# 通路英文标签
pathway_labels <- c(
  "低过激健康型" = "Low-hyperactivation",
  "对抗-枯竭混合型" = "Aversion-Exhaustion",
  "痴固着主导型" = "Perseveration",
  "纯生理过激型" = "Hyperactivation"
)
# 系统英文标签
system_labels <- c(
  "代谢" = "Metabolic",
  "心血管" = "Cardiovascular",
  "肾脏" = "Renal", 
  "肝脏" = "Hepatic",
  "神经精神" = "Neuropsychiatric"
)
# 疾病英文标签
disease_labels <- c(
  "糖尿病" = "Diabetes",
  "肥胖" = "Obesity",
  "高胆固醇" = "High cholesterol",
  "高血压" = "Hypertension",
  "心血管疾病" = "CVD",
  "心力衰竭" = "Heart failure",
  "中风" = "Stroke",
  "心肌梗死" = "Myocardial infarction",
  "慢性肾病" = "CKD",
  "ALT升高" = "Elevated ALT",
  "抑郁" = "Depression",
  "自评健康差" = "Poor self-rated health"
)
# ============================================================================
# 2. 配置路径
# ============================================================================
clean_dir <- "C:/NHANES_Data/CLEAN"
results_dir <- file.path(clean_dir, "results", "paper3")
LOG_DIR <- file.path(clean_dir, "logs")
# 创建结果目录
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
# 启动日志记录
log_file <- file.path(LOG_DIR, paste0("16_paper3_alpha_L_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 16_paper3_alpha_analysis_L.R\n")
cat("描述: NHANES 2021-2023 (L周期) 论文3 - α因子正式分析（整合版）\n")
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
# 3. 加载数据
# ============================================================================
cat("1. 加载数据...\n")
data_file <- file.path(clean_dir, "analysis_dataset_subset.rds")
if(!file.exists(data_file)) {
  stop("错误: analysis_dataset_subset.rds 不存在! 请先运行08_final_merge_L.R")
}
data_full <- readRDS(data_file)
cat(sprintf("已加载: analysis_dataset_subset.rds\n"))
cat(sprintf("样本量: %d\n", nrow(data_full)))
cat(sprintf("变量数: %d\n\n", ncol(data_full)))
data <- data_full
cat(sprintf("分析样本: %d\n\n", nrow(data)))
# ============================================================================
# 4. 创建设计对象
# ============================================================================
cat("2. 创建设计对象...\n")
design_int <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTINT2YR,
  nest = TRUE,
  data = data
)
design_mec <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = data
)
cat(" ✅ 访谈权重设计对象创建成功\n")
cat(" ✅ MEC权重设计对象创建成功\n\n")
# ============================================================================
# 5. 检查α因子变量
# ============================================================================
cat("3. 检查α因子变量...\n")
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
alpha_vars <- alpha_vars[alpha_vars %in% names(data)]
if(length(alpha_vars) == 0) {
  stop("错误: 没有找到α因子变量，请先运行05_alpha_factors_L.R")
}
cat(sprintf(" ✅ 找到α因子变量: %s\n\n", paste(alpha_vars, collapse = ", ")))
alpha_desc <- data.frame(
  Variable = alpha_labels[alpha_vars],
  Mean = sapply(alpha_vars, function(v) mean(data[[v]], na.rm = TRUE)),
  SD = sapply(alpha_vars, function(v) sd(data[[v]], na.rm = TRUE)),
  Missing_Pct = sapply(alpha_vars, function(v) mean(is.na(data[[v]])) * 100)
)
print(alpha_desc)
write.csv(alpha_desc, file.path(results_dir, "alpha_descriptives.csv"), row.names = FALSE)
cat(" ✅ 描述统计已保存\n\n")
# ============================================================================
# 6. 定义多系统疾病结局（8大系统）
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("4. 定义多系统疾病结局\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
# 定义疾病列表（所有系统）
disease_definitions <- list(
  # 代谢系统
  "diabetes" = list(var = "diabetes_doctor", label = "Diabetes", system = "代谢", type = "binary"),
  "obesity" = list(var = "obesity", label = "Obesity", system = "代谢", type = "binary"),
  "high_cholesterol" = list(var = "high_cholesterol_doctor", label = "High cholesterol", system = "代谢", type = "binary"),
  # 心血管系统
  "hypertension" = list(var = "hypertension", label = "Hypertension", system = "心血管", type = "binary"),
  "cvd" = list(var = "any_cvd", label = "CVD", system = "心血管", type = "binary"),
  "heart_failure" = list(var = "heart_failure", label = "Heart failure", system = "心血管", type = "binary"),
  "stroke" = list(var = "stroke", label = "Stroke", system = "心血管", type = "binary"),
  "heart_attack" = list(var = "heart_attack", label = "Myocardial infarction", system = "心血管", type = "binary"),
  # 肾脏系统
  "ckd" = list(var = "ckd_flag", label = "CKD", system = "肾脏", type = "binary"),
  # 肝脏系统
  "elevated_alt" = list(var = "elevated_alt_traditional", label = "Elevated ALT", system = "肝脏", type = "binary"),
  # 神经精神系统
  "depression" = list(var = "depression", label = "Depression", system = "神经精神", type = "binary"),
  "poor_health" = list(var = "poor_self_rated_health", label = "Poor self-rated health", system = "神经精神", type = "binary")
)
# 检查可用变量
available_diseases <- list()
cat("\n可用疾病结局:\n")
for(disease in names(disease_definitions)) {
  def <- disease_definitions[[disease]]
  if(def$var %in% names(data)) {
    available_diseases[[disease]] <- def
    event_count <- sum(data[[def$var]] == 1, na.rm = TRUE)
    cat(sprintf("  ✅ %s [%s]: %s (n=%d)\n", 
                disease, system_labels[def$system], def$label, event_count))
  }
}
# ============================================================================
# 7. 分析1：α因子与多系统疾病的关联（Table 3）
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("5. 分析1：α因子与多系统疾病的关联（Table 3）\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
alpha_vars_all <- alpha_vars
alpha_disease_results <- data.frame()
for(disease in names(available_diseases)) {
  def <- available_diseases[[disease]]
  disease_var <- def$var
  system <- def$system
  disease_label_en <- def$label
  # 根据系统选择安全α因子
  if(system %in% c("代谢")) {
    safe_alphas <- c("alpha1", "alpha2", "alpha4")  # 代谢系统排除α3
    cat(sprintf("\n%s: 使用α1/α2/α4 (排除α3避免循环论证)\n", disease_label_en))
  } else {
    safe_alphas <- alpha_vars_all  # 其他系统保留所有α因子
    cat(sprintf("\n%s: 使用所有α因子\n", disease_label_en))
  }
  event_count <- sum(data[[disease_var]] == 1, na.rm = TRUE)
  if(event_count < 20) {
    cat("  跳过: 事件数不足\n")
    next
  }
  formula <- as.formula(paste0(disease_var, " ~ ", 
                                paste(safe_alphas, collapse = " + "),
                                " + RIDAGEYR + RIAGENDR"))
  model <- tryCatch({
    svyglm(formula, design = design_mec, family = quasibinomial())
  }, error = function(e) {
    cat("  模型失败:", e$message, "\n")
    return(NULL)
  })
  if(!is.null(model)) {
    coefs <- summary(model)$coefficients
    for(alpha in safe_alphas) {
      if(alpha %in% rownames(coefs)) {
        alpha_disease_results <- rbind(alpha_disease_results, data.frame(
          System = system_labels[system],
          Disease = disease_label_en,
          Alpha_Factor = alpha_labels[alpha],
          OR = exp(coefs[alpha, "Estimate"]),
          CI_lower = exp(coefs[alpha, "Estimate"] - 1.96 * coefs[alpha, "Std. Error"]),
          CI_upper = exp(coefs[alpha, "Estimate"] + 1.96 * coefs[alpha, "Std. Error"]),
          P_value = coefs[alpha, "Pr(>|t|)"],
          Event_n = event_count,
          stringsAsFactors = FALSE
        ))
      }
    }
    cat("  ✅ 完成\n")
  }
}
# 获取痴固着型α因子均值
cat("\n获取痴固着型α因子均值...\n")
perseveration_alphas <- data.frame(
  System = "Perseveration type mean",
  Disease = "Alpha factor mean",
  `α₁ Metacognitive` = NA, 
  `α₂ Emotional Regulation` = NA, 
  `α₃ Systemic Coordination` = NA, 
  `α₄ Goal Efficacy` = NA,
  check.names = FALSE
)
if("pathway_cluster" %in% names(data)) {
  # 痴固着主导型是聚类3
  perseveration_data <- data %>% filter(pathway_cluster == 3)
  if(nrow(perseveration_data) > 0) {
    perseveration_alphas <- data.frame(
      System = "Perseveration type mean",
      Disease = "Alpha factor mean",
      `α₁ Metacognitive` = round(mean(perseveration_data$alpha1, na.rm = TRUE), 3),
      `α₂ Emotional Regulation` = round(mean(perseveration_data$alpha2, na.rm = TRUE), 3),
      `α₃ Systemic Coordination` = round(mean(perseveration_data$alpha3, na.rm = TRUE), 3),
      `α₄ Goal Efficacy` = round(mean(perseveration_data$alpha4, na.rm = TRUE), 3),
      check.names = FALSE
    )
    cat("\n痴固着主导型α因子均值:\n")
    print(perseveration_alphas)
  } else {
    cat("\n⚠️ 没有找到痴固着主导型数据\n")
  }
}
# 生成完整的Table 3（宽格式）
cat("\n生成完整的Table 3（包含痴固着型均值）...\n")
if(nrow(alpha_disease_results) > 0) {
  # 转换为宽格式
  table3 <- alpha_disease_results %>%
    select(System, Disease, Alpha_Factor, OR) %>%
    pivot_wider(names_from = Alpha_Factor, values_from = OR, id_cols = c(System, Disease))
  # 确保列名顺序一致
  expected_cols <- c("α₁ Metacognitive", "α₂ Emotional Regulation", 
                     "α₃ Systemic Coordination", "α₄ Goal Efficacy")
  for(col in expected_cols) {
    if(!col %in% names(table3)) {
      table3[[col]] <- NA
    }
  }
  # 按正确顺序排列
  table3 <- table3 %>%
    select(System, Disease, all_of(expected_cols))
  # 添加痴固着型均值行（使用相同的列名）
  table3_with_persev <- bind_rows(
    table3,
    perseveration_alphas
  )
  # 保存完整的Table 3
  write.csv(table3_with_persev, 
            file.path(results_dir, "Table3.csv"), 
            row.names = FALSE)
  cat("✅ 完整的Table 3已保存: Table3.csv\n")
} else {
  cat("⚠️ alpha_disease_results 为空，无法生成 Table 3\n")
}
# ============================================================================
# 8. 分析2：α因子与健康结局的关联（原分析1保留）
# ============================================================================
cat("\n========================================================\n")
cat("6. 分析2：α因子与健康结局的关联\n")
cat("========================================================\n")
outcomes <- list(
  "phq9_total" = "PHQ-9 total (continuous)",
  "depression" = "Depression (binary)",
  "poor_self_rated_health" = "Poor self-rated health",
  "hs_crp_mgl" = "Inflammation (continuous)",
  "BPXOPLS1" = "Resting heart rate"
)
alpha_outcome_continuous <- data.frame()
alpha_outcome_binary <- data.frame()
for(i in 1:length(outcomes)) {
  outcome_name <- names(outcomes)[i]
  outcome_label <- outcomes[[i]]
  if(!outcome_name %in% names(data)) {
    cat(sprintf("\n⚠️ 跳过: %s - 变量不存在\n", outcome_label))
    next
  }
  cat(sprintf("\n分析: %s\n", outcome_label))
  if(outcome_name %in% c("hs_crp_mgl", "BPXOPLS1", "BMXBMI")) {
    current_design <- design_mec
  } else {
    current_design <- design_int
  }
  formula <- as.formula(paste0(outcome_name, " ~ ", 
                                paste(alpha_vars, collapse = " + "),
                                " + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR"))
  if(outcome_name %in% c("depression", "poor_self_rated_health")) {
    model <- tryCatch({
      svyglm(formula, design = current_design, family = quasibinomial())
    }, error = function(e) {
      cat(" 模型失败:", e$message, "\n")
      return(NULL)
    })
    if(!is.null(model)) {
      coefs <- summary(model)$coefficients
      for(alpha in alpha_vars) {
        if(alpha %in% rownames(coefs)) {
          tmp <- data.frame(
            Outcome = outcome_label,
            Alpha_Factor = alpha_labels[alpha],
            OR = exp(coefs[alpha, "Estimate"]),
            CI_lower = exp(coefs[alpha, "Estimate"] - 1.96 * coefs[alpha, "Std. Error"]),
            CI_upper = exp(coefs[alpha, "Estimate"] + 1.96 * coefs[alpha, "Std. Error"]),
            P_value = coefs[alpha, "Pr(>|t|)"],
            Model_Type = "Logistic regression",
            stringsAsFactors = FALSE
          )
          alpha_outcome_binary <- rbind(alpha_outcome_binary, tmp)
        }
      }
      cat(" ✅ 二分类模型成功\n")
    }
  } else {
    model <- tryCatch({
      svyglm(formula, design = current_design)
    }, error = function(e) {
      cat(" 模型失败:", e$message, "\n")
      return(NULL)
    })
    if(!is.null(model)) {
      coefs <- summary(model)$coefficients
      for(alpha in alpha_vars) {
        if(alpha %in% rownames(coefs)) {
          tmp <- data.frame(
            Outcome = outcome_label,
            Alpha_Factor = alpha_labels[alpha],
            Beta = coefs[alpha, "Estimate"],
            SE = coefs[alpha, "Std. Error"],
            P_value = coefs[alpha, "Pr(>|t|)"],
            Model_Type = "Linear regression",
            stringsAsFactors = FALSE
          )
          alpha_outcome_continuous <- rbind(alpha_outcome_continuous, tmp)
        }
      }
      cat(" ✅ 连续模型成功\n")
    }
  }
}
if(nrow(alpha_outcome_continuous) > 0) {
  write.csv(alpha_outcome_continuous,
            file.path(results_dir, "alpha_outcome_continuous.csv"),
            row.names = FALSE)
  cat("\n✅ 已保存连续结果: alpha_outcome_continuous.csv\n")
}
if(nrow(alpha_outcome_binary) > 0) {
  write.csv(alpha_outcome_binary,
            file.path(results_dir, "alpha_outcome_binary.csv"),
            row.names = FALSE)
  cat("✅ 已保存二分类结果: alpha_outcome_binary.csv\n")
}
alpha_outcome_all <- bind_rows(alpha_outcome_continuous, alpha_outcome_binary)
write.csv(alpha_outcome_all,
          file.path(results_dir, "eTable4.csv"),
          row.names = FALSE)
cat("✅ 已保存全部结果: eTable4.csv\n\n")
# ============================================================================
# 9. 分析3：α3的独特价值（预测身心失联）
# ============================================================================
cat("\n========================================================\n")
cat("7. 分析3：α3的独特价值（预测身心失联）\n")
cat("========================================================\n")
data <- data %>%
  mutate(
    metabolic_inflammation_decoupling = as.numeric(
      BMXBMI >= 30 & hs_crp_mgl < 3 & !is.na(BMXBMI) & !is.na(hs_crp_mgl)
    ),
    symptom_function_separation = as.numeric(
      phq9_total >= 10 & poor_self_rated_health == 0 & !is.na(phq9_total) & !is.na(poor_self_rated_health)
    ),
    autonomic_perception_separation = as.numeric(
      BPXOPLS1 > 90 & DPQ040 < 2 & !is.na(BPXOPLS1) & !is.na(DPQ040)
    )
  )
design_mec <- update(design_mec,
                     metabolic_inflammation_decoupling = data$metabolic_inflammation_decoupling,
                     symptom_function_separation = data$symptom_function_separation,
                     autonomic_perception_separation = data$autonomic_perception_separation)
decoupling_outcomes <- c(
  "metabolic_inflammation_decoupling" = "Metabolic-inflammation decoupling",
  "symptom_function_separation" = "Symptom-function separation",
  "autonomic_perception_separation" = "Autonomic-perception separation"
)
decoupling_results <- data.frame()
for(outcome_name in names(decoupling_outcomes)) {
  cat(sprintf("\n分析: %s\n", decoupling_outcomes[outcome_name]))
  formula2 <- as.formula(paste0(outcome_name, " ~ alpha1 + alpha2 + alpha3 + alpha4 + RIDAGEYR + RIAGENDR"))
  model2 <- svyglm(formula2, design = design_mec, family = quasibinomial())
  test_alpha3 <- regTermTest(model2, "alpha3", method = "Wald")
  coefs <- summary(model2)$coefficients
  if("alpha3" %in% rownames(coefs)) {
    decoupling_results <- rbind(decoupling_results, data.frame(
      Outcome = decoupling_outcomes[outcome_name],
      Alpha3_OR = exp(coefs["alpha3", "Estimate"]),
      CI_lower = exp(coefs["alpha3", "Estimate"] - 1.96 * coefs["alpha3", "Std. Error"]),
      CI_upper = exp(coefs["alpha3", "Estimate"] + 1.96 * coefs["alpha3", "Std. Error"]),
      P_value = coefs["alpha3", "Pr(>|t|)"],
      Wald_F = test_alpha3$Ftest,
      Wald_p = test_alpha3$p
    ))
  }
}
if(nrow(decoupling_results) > 0) {
  write.csv(decoupling_results, file.path(results_dir, "alpha3_decoupling.csv"), row.names = FALSE)
  cat("\n✅ 已保存: alpha3_decoupling.csv\n")
  print(decoupling_results)
}
# ============================================================================
# 10. 分析4：α因子的增量预测价值
# ============================================================================
cat("\n========================================================\n")
cat("8. 分析4：α因子的增量预测价值\n")
cat("========================================================\n")
cat("\n预测抑郁 (PHQ-9总分)\n")
model_step1 <- svyglm(phq9_total ~ RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                      design = design_int)
model_step2 <- svyglm(phq9_total ~ RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR +
                        BMXBMI + hs_crp_mgl + BPXOPLS1,
                      design = design_mec)
model_step3 <- svyglm(phq9_total ~ RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR +
                        BMXBMI + hs_crp_mgl + BPXOPLS1 +
                        alpha1 + alpha2 + alpha3 + alpha4,
                      design = design_mec)
cat("\n【AIC比较】\n")
aic_values <- c(AIC(model_step1), AIC(model_step2), AIC(model_step3))
aic_results <- data.frame(
  Model = c("Model 1: Demographics", "Model 2: +Physiological", "Model 3: +Alpha factors"),
  AIC = aic_values,
  Delta_AIC = c(0, aic_values[2] - aic_values[1], aic_values[3] - aic_values[2])
)
print(aic_results)
write.csv(aic_results, file.path(results_dir, "aic_comparisons.csv"), row.names = FALSE)
cat("\n【BIC比较】\n")
bic_values <- c(
  BIC(model_step1, maximal = model_step3),
  BIC(model_step2, maximal = model_step3),
  BIC(model_step3, maximal = model_step3)
)
bic_results <- data.frame(
  Model = c("Model 1: Demographics", "Model 2: +Physiological", "Model 3: +Alpha factors"),
  BIC = bic_values,
  Delta_BIC = c(0, bic_values[2] - bic_values[1], bic_values[3] - bic_values[2])
)
print(bic_results)
write.csv(bic_results, file.path(results_dir, "bic_comparisons.csv"), row.names = FALSE)
cat("\n【嵌套模型比较】\n")
cat("\n检验生理指标的增量贡献 (模型1 vs 模型2):\n")
phys_vars <- c("BMXBMI", "hs_crp_mgl", "BPXOPLS1")
test_phys <- regTermTest(model_step2, phys_vars, method = "Wald")
print(test_phys)
cat("\n检验α因子的增量贡献 (模型2 vs 模型3):\n")
test_alpha <- regTermTest(model_step3, alpha_vars, method = "Wald")
print(test_alpha)
incremental_tests <- data.frame(
  Test = c("Physiological indicators", "Alpha factors"),
  F_value = c(test_phys$Ftest, test_alpha$Ftest),
  df = c(test_phys$df, test_alpha$df),
  P_value = c(test_phys$p, test_alpha$p)
)
write.csv(incremental_tests, file.path(results_dir, "incremental_tests.csv"), row.names = FALSE)
cat("\n【效应量比较】\n")
n <- nrow(model_step3$data)
df1_phys <- 3
df1_alpha <- 4
f2_phys <- (test_phys$Ftest * df1_phys) / (n - df1_phys - 1)
f2_alpha <- (test_alpha$Ftest * df1_alpha) / (n - df1_alpha - 1)
effect_sizes <- data.frame(
  Effect = c("Physiological indicators", "Alpha factors"),
  F_value = c(test_phys$Ftest, test_alpha$Ftest),
  df = c(df1_phys, df1_alpha),
  Cohen_f2 = round(c(f2_phys, f2_alpha), 4)
)
print(effect_sizes)
write.csv(effect_sizes, file.path(results_dir, "effect_sizes.csv"), row.names = FALSE)
cat("\n【分析4结论】\n")
if(test_alpha$p < 0.05) {
  cat(sprintf("✅ α因子在生理指标之上具有显著的增量预测价值\n F = %.3f, p = %.4f\n", 
              test_alpha$Ftest, test_alpha$p))
} else {
  cat("⚠️ α因子的增量预测价值未达到统计学显著性\n")
}
# ============================================================================
# 8a. 新增：增量预测价值的优化分析（JAMA推荐格式）
# ============================================================================
cat("\n========================================================\n")
cat("8a. 增量预测价值优化分析（JAMA推荐格式）\n")
cat("========================================================\n")
# 计算R²（如果原代码中没有计算）
r2_model1 <- cor(model_step1$fitted.values, model_step1$y, use = "complete.obs")^2
r2_model2 <- cor(model_step2$fitted.values, model_step2$y, use = "complete.obs")^2
r2_model3 <- cor(model_step3$fitted.values, model_step3$y, use = "complete.obs")^2
# 计算增量结果
incremental_predictive_value <- data.frame(
  Comparison = c("Adding physiological indicators", "Adding alpha factors"),
  Delta_AIC = c(aic_values[2] - aic_values[1], aic_values[3] - aic_values[2]),  # 改这里
  Delta_BIC = c(bic_values[2] - bic_values[1], bic_values[3] - bic_values[2]),  # 改这里
  F_value = c(test_phys$Ftest, test_alpha$Ftest),
  df = c(test_phys$df, test_alpha$df),
  P_value = c(test_phys$p, test_alpha$p),
  Cohen_f2 = round(c(f2_phys, f2_alpha), 4)
)
# 保存增量结果
write.csv(incremental_predictive_value, 
          file.path(results_dir, "eTable25.csv"), 
          row.names = FALSE)
# 生成详细比较表格
supp_table_model_comparison <- data.frame(
  Model = c(
    "Model 1: Demographics",
    "Model 2: Model 1 + Physiological indicators",
    "Model 3: Model 2 + Alpha factors"
  ),
  Covariates = c(
    "Age, sex, education, poverty",
    "Age, sex, education, poverty, BMI, hs-CRP, heart rate",
    "Age, sex, education, poverty, BMI, hs-CRP, heart rate, α₁-α₄"
  ),
  AIC = round(aic_values, 2),
  BIC = round(bic_values, 2),
  R2 = round(c(r2_model1, r2_model2, r2_model3), 3)
)
# 保存详细表格
write.csv(supp_table_model_comparison, 
          file.path(results_dir, "eTable25_details.csv"), 
          row.names = FALSE)
cat("\n✅ 增量预测价值优化分析完成\n")
cat("   新增文件: eTable25.csv\n")
cat("   新增文件: eTable25_details.csv\n")
# ============================================================================
# 11. 分析5：系统特异性疾病负担 - 修正循环论证
# ============================================================================
cat("\n========================================================\n")
cat("9. 分析5：系统特异性疾病负担\n")
cat("========================================================\n")
# 创建系统疾病负担变量
disease_systems <- list(
  metabolic = c("diabetes_doctor", "obesity", "high_cholesterol_doctor"),
  cardiovascular = c("hypertension", "any_cvd", "heart_failure", "stroke", "heart_attack"),
  renal = c("ckd_flag"),
  hepatic = c("elevated_alt_traditional"),
  mental = c("depression", "poor_self_rated_health")
)
for(sys in names(disease_systems)) {
  vars <- disease_systems[[sys]][disease_systems[[sys]] %in% names(data)]
  if(length(vars) > 0) {
    data[[paste0("disease_count_", sys)]] <- rowSums(data[, vars], na.rm = TRUE)
    design_mec$variables[[paste0("disease_count_", sys)]] <- data[[paste0("disease_count_", sys)]]
  }
}
system_burden <- data.frame()
for(alpha in alpha_vars) {
  for(sys in names(disease_systems)) {
    # 代谢系统排除α3（避免循环论证）
    if(sys == "metabolic" && alpha == "alpha3") {
      cat(sprintf("跳过 α3 对 %s 系统的分析（避免循环论证）\n", sys))
      next
    }
    count_var <- paste0("disease_count_", sys)
    if(!count_var %in% names(design_mec$variables)) next
    formula <- as.formula(paste0(count_var, " ~ ", alpha, " + RIDAGEYR + RIAGENDR"))
    model <- tryCatch({
      svyglm(formula, design = design_mec)
    }, error = function(e) NULL)
    if(!is.null(model)) {
      coefs <- summary(model)$coefficients
      if(alpha %in% rownames(coefs)) {
        system_burden <- rbind(system_burden, data.frame(
          Alpha_Factor = alpha_labels[alpha],
          System = system_labels[sys],
          Beta = coefs[alpha, "Estimate"],
          SE = coefs[alpha, "Std. Error"],
          P_value = coefs[alpha, "Pr(>|t|)"],
          CI_lower = coefs[alpha, "Estimate"] - 1.96 * coefs[alpha, "Std. Error"],
          CI_upper = coefs[alpha, "Estimate"] + 1.96 * coefs[alpha, "Std. Error"]
        ))
      }
    }
  }
}
if(nrow(system_burden) > 0) {
  write.csv(system_burden, file.path(results_dir, "eTable5.csv"), row.names = FALSE)
  cat("\n✅ 已保存: eTable5.csv\n")
  print(system_burden)
}
# ============================================================================
# 12. 分析6：二联中介分析（Figure 4）
# ============================================================================
cat("\n========================================================\n")
cat("10. 分析6：二联中介分析（Figure 4）\n")
cat("========================================================\n")
# 准备中介数据 - 使用心血管疾病负担（与α3无关）
med_data <- data %>%
  select(disease_count_cardiovascular,
         alpha2, alpha4,
         BMXBMI, hs_crp_mgl,
         RIDAGEYR, RIAGENDR) %>%
  rename(age = RIDAGEYR, sex = RIAGENDR, bmi = BMXBMI,
         crp = hs_crp_mgl) %>%
  mutate(log_crp = log(crp + 0.1)) %>%
  filter(complete.cases(.))
cat("完整案例样本量:", nrow(med_data), "\n")
# 模型1: BMI → α2 → 心血管疾病负担
cat("\n--- 模型1: BMI → α2 → Cardiovascular disease burden ---\n")
model1 <- '
  disease_count_cardiovascular ~ c*bmi + age + sex
  alpha2 ~ a*bmi + age + sex
  disease_count_cardiovascular ~ b*alpha2 + age + sex
  indirect := a*b
  total := c + (a*b)
  prop_mediated := indirect / total * 100
'
fit1 <- sem(model1, data = med_data, estimator = "MLR")
summary(fit1, fit.measures = TRUE, standardized = TRUE)
# 模型2: BMI → α4 → 心血管疾病负担
cat("\n--- 模型2: BMI → α4 → Cardiovascular disease burden ---\n")
model2 <- '
  disease_count_cardiovascular ~ c*bmi + age + sex
  alpha4 ~ a*bmi + age + sex
  disease_count_cardiovascular ~ b*alpha4 + age + sex
  indirect := a*b
  total := c + (a*b)
  prop_mediated := indirect / total * 100
'
fit2 <- sem(model2, data = med_data, estimator = "MLR")
summary(fit2, fit.measures = TRUE, standardized = TRUE)
# 模型3: CRP → α2 → 心血管疾病负担
cat("\n--- 模型3: CRP → α2 → Cardiovascular disease burden ---\n")
model3 <- '
  disease_count_cardiovascular ~ c*log_crp + age + sex
  alpha2 ~ a*log_crp + age + sex
  disease_count_cardiovascular ~ b*alpha2 + age + sex
  indirect := a*b
  total := c + (a*b)
  prop_mediated := indirect / total * 100
'
fit3 <- sem(model3, data = med_data, estimator = "MLR")
summary(fit3, fit.measures = TRUE, standardized = TRUE)
# 收集结果
mediation_results <- bind_rows(
  parameterEstimates(fit1) %>% filter(op == ":=") %>% mutate(Model = "BMI → α₂ → CVD"),
  parameterEstimates(fit2) %>% filter(op == ":=") %>% mutate(Model = "BMI → α₄ → CVD"),
  parameterEstimates(fit3) %>% filter(op == ":=") %>% mutate(Model = "CRP → α₂ → CVD")
)
write.csv(mediation_results, file.path(results_dir, "eTable6.csv"), row.names = FALSE)
cat("\n✅ 已保存: eTable6.csv\n")
# ============================================================================
# 13. 分析7：α因子与四通路的交互
# ============================================================================
cat("\n========================================================\n")
cat("11. 分析7：α因子与四通路的交互\n")
cat("========================================================\n")
if("pathway_cluster" %in% names(data)) {
  data$pathway_cluster <- factor(data$pathway_cluster)
  design_int <- update(design_int, pathway_cluster = data$pathway_cluster)
  interaction_results <- data.frame()
  for(alpha in alpha_vars) {
    formula <- as.formula(paste0("phq9_total ~ ", alpha, " * pathway_cluster + RIDAGEYR + RIAGENDR"))
    model <- svyglm(formula, design = design_int)
    f_test <- regTermTest(model, ~pathway_cluster, method = "Wald")
    interaction_results <- rbind(interaction_results, data.frame(
      Alpha_Factor = alpha_labels[alpha],
      F_value = f_test$Ftest,
      P_value = f_test$p,
      df = f_test$df
    ))
    cat("\n", alpha_labels[alpha], "的简单斜率:\n")
    emm <- emmeans(model, ~ pathway_cluster, at = setNames(list(mean(data[[alpha]], na.rm = TRUE)), alpha))
    print(pairs(emm))
  }
  write.csv(interaction_results, file.path(results_dir, "Table4.csv"), row.names = FALSE)
  cat("\n✅ 已保存: Table4.csv\n")
  print(interaction_results)
}
# ============================================================================
# 14. 分析8：α因子与HCF分型的交互
# ============================================================================
cat("\n========================================================\n")
cat("12. 分析8：α因子与HCF分型的交互\n")
cat("========================================================\n")
if("HCF_type" %in% names(data)) {
  data$HCF_type <- factor(data$HCF_type)
  design_int <- update(design_int, HCF_type = data$HCF_type)
  hcf_interaction <- data.frame()
  for(alpha in alpha_vars) {
    formula <- as.formula(paste0("phq9_total ~ ", alpha, " * HCF_type + RIDAGEYR + RIAGENDR"))
    model <- svyglm(formula, design = design_int)
    f_test <- regTermTest(model, ~HCF_type, method = "Wald")
    hcf_interaction <- rbind(hcf_interaction, data.frame(
      Alpha_Factor = alpha_labels[alpha],
      F_value = f_test$Ftest,
      P_value = f_test$p,
      df = f_test$df
    ))
  }
  # 保存为 eTable 18
  write.csv(hcf_interaction, file.path(results_dir, "eTable18.csv"), row.names = FALSE)
  cat("\n✅ 已保存: eTable18.csv (α因子×HCF分型交互)\n")
  print(hcf_interaction)
}
# ============================================================================
# 15. 可视化：Figure 2 - 四通路α因子分布图
# ============================================================================
cat("\n========================================================\n")
cat("13. 生成Figure 2 - 四通路α因子分布图\n")
cat("========================================================\n")
paper2_results_dir <- file.path(clean_dir, "results", "paper2")
alpha_pathway_file <- file.path(paper2_results_dir, "table_S2_alpha_by_cluster.csv")
if(file.exists(alpha_pathway_file)) {
  alpha_pathway_data <- read.csv(alpha_pathway_file)
  # 直接使用现有的列 - 保持原来正确的逻辑
  alpha_pathway_long <- alpha_pathway_data %>%
    pivot_longer(cols = starts_with("alpha"),
                 names_to = "alpha_factor",
                 values_to = "value") %>%
    mutate(
      alpha_label = case_when(
        alpha_factor == "alpha1" ~ "α₁ Metacognitive",
        alpha_factor == "alpha2" ~ "α₂ Emotional Regulation",
        alpha_factor == "alpha3" ~ "α₃ Systemic Coordination",
        alpha_factor == "alpha4" ~ "α₄ Goal Efficacy"
      ),
      # 直接使用 pathway_cluster_en，它已经是英文
      pathway_cluster = factor(pathway_cluster_en,
                               levels = c("Low-hyperactivation", 
                                          "Aversion-Exhaustion", 
                                          "Perseveration", 
                                          "Hyperactivation"))
    )
  # 检查数据
  cat("数据行数:", nrow(alpha_pathway_long), "\n")
  p_alpha_pathway <- ggplot(alpha_pathway_long,
                            aes(x = pathway_cluster, y = value, fill = alpha_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Figure 2. Alpha Factor Distribution Across Pathways",
         x = "Pathway Type", y = "Alpha Factor Mean (Standardized)",
         fill = "Alpha Factor") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  # 显示图形
  print(p_alpha_pathway)
# 保存图表（PDF和PNG）
ggsave(file.path(results_dir, "Figure2.pdf"), p_alpha_pathway, width = 10, height = 6)
ggsave(file.path(results_dir, "Figure2.png"), p_alpha_pathway, width = 10, height = 6, dpi = 300)
# 保存数据（CSV）
write.csv(alpha_pathway_long, file.path(results_dir, "alpha_pathway_long.csv"), row.names = FALSE)
cat(" ✅ Figure 2 saved: Figure2.pdf/.png\n")
cat(" ✅ alpha_pathway_long.csv saved\n")
}
# ============================================================================
# 16. 可视化：Figure 4 - 中介网络图（最终优化版）
# ============================================================================
cat("\n========================================================\n")
cat("14. 生成Figure 4 - 中介网络图（最终优化版）\n")
cat("========================================================\n")
# 加载必要的包
library(qgraph)
# 定义节点
nodes <- c("BMI", "CRP", "α₂", "α₄", "CVD")
# 定义边 - 使用深绿色
edges <- data.frame(
  from = c("BMI", "BMI", "CRP", "α₂", "α₄"),
  to = c("α₂", "α₄", "α₂", "CVD", "CVD"),
  weight = c(0.0102, 0.00066, 0.0684, 1.0, 0.8),
  color = c("#006400", "#006400", "#006400", "#006400", "#006400")
)
# 创建邻接矩阵
adj_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes))
rownames(adj_matrix) <- nodes
colnames(adj_matrix) <- nodes
for(i in 1:nrow(edges)) {
  adj_matrix[edges$from[i], edges$to[i]] <- edges$weight[i]
}
# 绘制网络图 - PDF
pdf(file.path(results_dir, "Figure4.pdf"), 
    width = 10, height = 8,
    family = "Helvetica")
qgraph(adj_matrix,
       layout = "spring",
       labels = nodes,
       label.cex = 1.4,        # 减小节点标签（原1.6）
       label.font = 2,
       color = c("lightblue", "lightblue", "lightgreen", "lightgreen", "lightcoral"),
       borders = TRUE,
       border.width = 2,
       vsize = 9,               # 减小节点大小（原11）
       esize = 8,
       edge.color = edges$color,
       edge.width = 2.5,
       edge.labels = c("12.8%", "3.0%", "19.5%", "", ""),
       edge.label.cex = 1.3,    # 稍微减小边标签
       edge.label.font = 2,
       edge.label.position = c(0.3, 0.7, 0.5, 0.5, 0.5),
       title = "Figure 4. α₂ Mediation Network",
       title.cex = 1.4,         # 减小标题（原1.6）
       title.font = 2,
       cut = 0.001,
       minimum = 0.0005,
       maximum = 1,
       details = TRUE,
       posCol = "#006400",
       fade = FALSE,
       # 增大底部文字的显示
       label.scale = FALSE,
       label.prop = 1,
       # 调整边距，给底部更多空间
       mar = c(6, 3, 3, 3))
dev.off()
cat(" ✅ Figure 4 PDF saved: Figure4.pdf\n")
# 绘制PNG（高分辨率）
png(file.path(results_dir, "Figure4.png"), 
    width = 1500, height = 1200, res = 250)
qgraph(adj_matrix,
       layout = "spring",
       labels = nodes,
       label.cex = 1.4,
       label.font = 2,
       color = c("lightblue", "lightblue", "lightgreen", "lightgreen", "lightcoral"),
       borders = TRUE,
       border.width = 2,
       vsize = 9,
       esize = 8,
       edge.color = edges$color,
       edge.width = 2.5,
       edge.labels = c("12.8%", "3.0%", "19.5%", "", ""),
       edge.label.cex = 1.3,
       edge.label.font = 2,
       edge.label.position = c(0.3, 0.7, 0.5, 0.5, 0.5),
       title = "Figure 4. α₂ Mediation Network",
       title.cex = 1.4,
       title.font = 2,
       cut = 0.001,
       minimum = 0.0005,
       maximum = 1,
       details = TRUE,
       posCol = "#006400",
       fade = FALSE,
       label.scale = FALSE,
       label.prop = 1,
       mar = c(6, 3, 3, 3))
dev.off()
cat(" ✅ Figure 4 PNG saved: Figure4.png\n")
# ============================================================================
# 17. 生成分析报告
# ============================================================================
cat("\n========================================================\n")
cat("15. 生成分析报告\n")
cat("========================================================\n")
report_file <- file.path(LOG_DIR, "16_paper3_report_L.txt")
sink(report_file)
cat("Alpha Factor Analysis Report\n")
cat("============================\n\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、Alpha factor descriptives\n")
print(alpha_desc)
cat("\n二、Alpha factors and multi-system diseases (Table 3)\n")
if(exists("table3") && nrow(table3) > 0) print(table3)
cat("\n三、Alpha factors predicting health outcomes\n")
if(exists("alpha_outcome_all") && nrow(alpha_outcome_all) > 0) print(head(alpha_outcome_all))
cat("\n四、α3 unique value\n")
if(exists("decoupling_results") && nrow(decoupling_results) > 0) print(decoupling_results)
cat("\n五、Incremental predictive value\n")
print(incremental_tests)
print(effect_sizes)
cat("\n六、System-specific burden\n")
if(exists("system_burden") && nrow(system_burden) > 0) print(system_burden)
cat("\n七、Mediation analysis\n")
if(exists("mediation_results") && nrow(mediation_results) > 0) print(mediation_results)
cat("\n八、Interaction with pathways\n")
if(exists("interaction_results") && nrow(interaction_results) > 0) print(interaction_results)
cat("\n九、Interaction with HCF types\n")
if(exists("hcf_interaction") && nrow(hcf_interaction) > 0) print(hcf_interaction)
sink()
cat(" ✅ 分析报告已保存\n\n")
# ============================================================================
# 18. 保存会话信息
# ============================================================================
cat("16. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "16_session_info_L.txt")
sink(session_info_path)
cat("NHANES Paper 3 Alpha Factor Analysis Session Information\n")
cat("======================================================\n")
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
# 19. 保存R代码副本
# ============================================================================
cat("\n17. 保存R代码副本...\n")
scripts_dir <- file.path("C:/NHANES_Data", "scripts")
if (!dir.exists(scripts_dir)) {
  dir.create(scripts_dir, recursive = TRUE)
}
code_save_path <- file.path(scripts_dir, "16_paper3_alpha_analysis_L.R")
cat("\n⚠️  请手动将当前脚本保存到以下位置：\n")
cat(sprintf("   %s\n\n", code_save_path))
cat("   这是JAMA Psychiatry的明确要求：所有分析代码必须保存并公开。\n")
code_list_path <- file.path(LOG_DIR, "16_code_list_L.txt")
cat("脚本名称: 16_paper3_alpha_analysis_L.R\n", file = code_list_path)
cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = code_list_path, append = TRUE)
cat("建议保存位置:", code_save_path, "\n", file = code_list_path, append = TRUE)
cat(" ✅ 代码清单已保存\n")
# ============================================================================
# 20. 清理临时变量
# ============================================================================
cat("\n18. 清理临时变量...\n")
rm(list = setdiff(ls(), c("clean_dir", "results_dir", "LOG_DIR")))
gc()
# ============================================================================
# 21. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ Paper 3 alpha factor analysis complete!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results saved to:", results_dir, "\n")
cat("========================================================\n")
sink()
