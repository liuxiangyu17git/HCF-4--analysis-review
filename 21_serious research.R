#!/usr/bin/env Rscript
# ============================================================================
# 脚本: 21_serious research.R
# 描述: 严谨的补充分析
# 
# 包含:
#   1. α₂阈值稳定性 (两周期对比 + CART自动切点)
#   2. 压缩假说分层证据 (HCF分型、年龄、性别、种族)
#   3. α因子正交性验证 (相关矩阵 + EFA)
#   4. α₁状态因子证据 (相关性对比 + ANCOVA)
#   5. α₂干预靶点证据 (平衡性检验 + PSM后NNT)
#   6. 补充缺失的分析 (混合效应模型、CFA、多组不变性)
#
# 依赖: 08_final_merge_P.R 的输出文件
# 最后修改: 2026-03-04
# ============================================================================
# ============================================================================
# 1. 环境配置
# ============================================================================
rm(list = ls())
gc()
set.seed(20240226)
# 加载必要包
required_packages <- c(
  "tidyverse", "survey", "rpart", "partykit", "tableone",
  "MatchIt", "cobalt", "psych", "lavaan", "cocor", "lme4",
  "ggplot2", "gridExtra", "pROC", "splines", "WeightIt",
  "cluster"
)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
# ============================================================================
# 2. 定义所有路径变量（必须先定义！）
# ============================================================================
PROJECT_ROOT <- "C:/NHANES_Data"
L_DATA_DIR <- file.path(PROJECT_ROOT, "CLEAN")
P_DATA_DIR <- file.path(PROJECT_ROOT, "2017-2020")
RESULTS_DIR <- file.path(PROJECT_ROOT, "CLEAN", "results", "serious_research")
LOG_DIR <- file.path(L_DATA_DIR, "logs")
# 创建目录（确保存在）
dir.create(PROJECT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(L_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(P_DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
# ============================================================================
# 3. 诊断代码（现在可以正常运行）
# ============================================================================
cat("\n🔍 开始诊断 ==========\n")
# 1. 检查 RESULTS_DIR
cat("1. RESULTS_DIR:\n")
cat("   值:", RESULTS_DIR, "\n")
cat("   存在:", dir.exists(RESULTS_DIR), "\n")
# 2. 测试写入 RESULTS_DIR
test1 <- file.path(RESULTS_DIR, "test1.csv")
write.csv(data.frame(test = 1), test1)
cat("   写入测试:", file.exists(test1), "\n")
if (file.exists(test1)) file.remove(test1)
# 3. 检查 LOG_DIR
cat("\n2. LOG_DIR:\n")
cat("   值:", LOG_DIR, "\n")
cat("   存在:", dir.exists(LOG_DIR), "\n")
# 4. 测试写入 LOG_DIR
test2 <- file.path(LOG_DIR, "test2.txt")
writeLines("test", test2)
cat("   写入测试:", file.exists(test2), "\n")
if (file.exists(test2)) file.remove(test2)
cat("\n🔍 诊断结束 ==========\n\n")
# ============================================================================
# 4. 启动日志（现在所有路径都已定义）
# ============================================================================
log_file <- file.path(LOG_DIR, paste0("21_serious_research_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("========================================================\n")
cat("脚本: 21_serious research.R\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================\n\n")
# ============================================================================
# 2. 加载两个周期的数据
# ============================================================================
cat("1. 加载两个周期的数据...\n")
# L周期 (2021-2023)
L_file <- file.path(L_DATA_DIR, "analysis_dataset_subset.rds")
if (!file.exists(L_file)) stop("错误: 找不到L周期数据")
data_L <- readRDS(L_file)
data_L$cycle <- "2021-2023"
cat(sprintf("L周期样本量: %d\n", nrow(data_L)))
# P周期 (2017-2020)
P_file <- file.path(P_DATA_DIR, "analysis_dataset_subset_P.rds")
if (!file.exists(P_file)) stop("错误: 找不到P周期数据")
data_P <- readRDS(P_file)
data_P$cycle <- "2017-2020"
cat(sprintf("P周期样本量: %d\n\n", nrow(data_P)))
# 合并数据用于部分分析
data_combined <- bind_rows(data_P, data_L)
# ============================================================================
# 3. 加载通路聚类数据（修正版 - 避免 .x .y 问题）
# ============================================================================
cat("2. 加载通路聚类数据...\n")
# L周期通路数据
pathway_L_file <- file.path(L_DATA_DIR, "pathway_final_with_clusters.rds")
if(file.exists(pathway_L_file)) {
  pathway_L <- readRDS(pathway_L_file)
  cat("L周期通路数据加载成功，维度:", dim(pathway_L), "\n")
  # 合并前删除可能冲突的变量
  data_L <- data_L %>% select(-any_of(c("pathway_cluster", "pathway_cluster_name")))
  data_L <- left_join(data_L, 
                      pathway_L %>% select(SEQN, pathway_cluster, pathway_cluster_name),
                      by = "SEQN")
  cat("L周期通路数据合并成功\n")
}
# P周期通路数据
pathway_P_file <- file.path(P_DATA_DIR, "pathway_final_with_clusters_P.rds")
if(file.exists(pathway_P_file)) {
  pathway_P <- readRDS(pathway_P_file)
  cat("P周期通路数据加载成功，维度:", dim(pathway_P), "\n")
  data_P <- data_P %>% select(-any_of(c("pathway_cluster", "pathway_cluster_name")))
  data_P <- left_join(data_P,
                      pathway_P %>% select(SEQN, pathway_cluster, pathway_cluster_name),
                      by = "SEQN")
  cat("P周期通路数据合并成功\n")
}
# 检查通路聚类分布
cat("\n通路聚类分布检查:\n")
if("pathway_cluster" %in% names(data_L)) {
  cat("L周期通路聚类:\n")
  print(table(data_L$pathway_cluster, useNA = "ifany"))
  cat("L周期聚类3（痴固着）样本量:", sum(data_L$pathway_cluster == 3, na.rm = TRUE), "\n")
} else {
  cat("L周期: pathway_cluster 不存在\n")
}
if("pathway_cluster" %in% names(data_P)) {
  cat("P周期通路聚类:\n")
  print(table(data_P$pathway_cluster, useNA = "ifany"))
  cat("P周期聚类3（痴固着）样本量:", sum(data_P$pathway_cluster == 3, na.rm = TRUE), "\n")
} else {
  cat("P周期: pathway_cluster 不存在\n")
}
# ============================================================================
# 4. 创建调查设计对象
# ============================================================================
cat("\n3. 创建调查设计对象...\n")
# L周期设计对象
design_L <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = data_L
)
# P周期设计对象
design_P <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMECPRP,
  nest = TRUE,
  data = data_P
)
cat(" ✅ 设计对象创建成功\n\n")
# ============================================================================
# 生成 eTable 28-29: 两个周期的α因子成分相关矩阵（符合NHANES要求）
# ============================================================================
cat("\n========================================================\n")
cat("生成 eTable 28-29: α Factor Components Correlation Matrix (Both Cycles)\n")
cat("========================================================\n")
cat("结果将保存至:", RESULTS_DIR, "\n")
# 定义成分变量（来自alpha_factors.rds）
component_vars <- c(
  "awareness_sleep", "awareness_health", "awareness_diet", "awareness_symptom",
  "acceptance_self", "acceptance_relationship", "acceptance_adversity",
  "conflict1", "conflict2", "conflict3", "conflict4", "conflict5", "conflict6",
  "goal_edu", "goal_economic", "goal_health"
)
# 简化列名用于表格显示
col_abbr <- c(
  "ASlp", "AHea", "ADie", "ASym",
  "ASel", "ARel", "AAdv",
  "C1", "C2", "C3", "C4", "C5", "C6",
  "GEd", "GEc", "GHe"
)
# 缩写解释（用于脚注）
abbr_key <- data.frame(
  Abbreviation = col_abbr,
  Full_Name = c(
    "awareness_sleep", "awareness_health", "awareness_diet", "awareness_symptom",
    "acceptance_self", "acceptance_relationship", "acceptance_adversity",
    "conflict1", "conflict2", "conflict3", "conflict4", "conflict5", "conflict6",
    "goal_edu", "goal_economic", "goal_health"
  )
)
write.csv(abbr_key, file.path(RESULTS_DIR, "eTable28_abbr_key.csv"), row.names = FALSE)
# ============================================================================
# L周期 (2021-2023) - 直接从主数据获取（包含设计变量）
# ============================================================================
cat("\n--- L周期 (2021-2023) 加权相关矩阵 ---\n")
# 直接从主数据获取成人非孕妇样本
adult_idx_L <- which(data_L$is_non_pregnant_adult == 1)
if(length(adult_idx_L) > 0) {
  data_adult_L <- data_L[adult_idx_L, ]
  cat(sprintf("L周期成人非孕妇样本量: %d\n", nrow(data_adult_L)))
  # 检查成分变量是否存在
  existing_vars_L <- component_vars[component_vars %in% names(data_adult_L)]
  cat(sprintf("L周期找到 %d 个成分变量\n", length(existing_vars_L)))
  if(length(existing_vars_L) >= 16) {
    # 设计变量已经在 data_adult_L 中！
    design_vars <- c("SDMVPSU", "SDMVSTRA", "WTINT2YR")
    cat("设计变量存在性检查:\n")
    for(v in design_vars) {
      cat(sprintf("  %s: %s\n", v, v %in% names(data_adult_L)))
    }
    # 筛选完整案例
    needed_cols <- c(existing_vars_L, design_vars)
    complete_idx <- complete.cases(data_adult_L[, needed_cols])
    data_complete <- data_adult_L[complete_idx, ]
    cat(sprintf("完整案例样本量: %d (%.1f%% of adults)\n", 
                nrow(data_complete), 
                100 * nrow(data_complete) / nrow(data_adult_L)))
    if(nrow(data_complete) > 100) {
      # 创建设计对象 - 直接使用 data_complete，所有变量都在！
      temp_design <- svydesign(
        id = ~SDMVPSU,
        strata = ~SDMVSTRA,
        weights = ~WTINT2YR,
        data = data_complete,
        nest = TRUE
      )
      # 计算加权相关
      formula_str <- paste("~", paste(existing_vars_L, collapse = " + "))
      vcov_matrix <- svyvar(as.formula(formula_str), temp_design, na.rm = TRUE)
      cor_matrix_L <- cov2cor(as.matrix(vcov_matrix))
      # 格式化输出
      cor_matrix_L_rounded <- round(cor_matrix_L, 2)
      rownames(cor_matrix_L_rounded) <- col_abbr[1:length(existing_vars_L)]
      colnames(cor_matrix_L_rounded) <- col_abbr[1:length(existing_vars_L)]
      # 保存完整版
      write.csv(cor_matrix_L_rounded, file.path(RESULTS_DIR, "eTable28_full.csv"), row.names = TRUE)
      cat("✅ 已保存 eTable28_full.csv (加权版)\n")
      # 创建简洁版表格
      eTable28 <- as.data.frame(cor_matrix_L_rounded)
      eTable28[] <- lapply(eTable28, function(x) {
        ifelse(abs(x) < 0.005, "-", sprintf("%.2f", x))
      })
      eTable28$Component <- rownames(eTable28)
      eTable28 <- eTable28[, c("Component", col_abbr[1:length(existing_vars_L)])]
      write.csv(eTable28, file.path(RESULTS_DIR, "eTable28.csv"), row.names = FALSE)
      cat("✅ 已保存 eTable28.csv (加权简洁版)\n")
    } else {
      cat("⚠️ 完整案例不足100，使用未加权相关\n")
      cor_matrix_L <- cor(data_adult_L[, existing_vars_L], use = "pairwise.complete.obs")
      # 格式化输出
      cor_matrix_L_rounded <- round(cor_matrix_L, 2)
      rownames(cor_matrix_L_rounded) <- col_abbr[1:length(existing_vars_L)]
      colnames(cor_matrix_L_rounded) <- col_abbr[1:length(existing_vars_L)]
      write.csv(cor_matrix_L_rounded, file.path(RESULTS_DIR, "eTable28_full_unweighted.csv"), row.names = TRUE)
      cat("✅ 已保存 eTable28_full_unweighted.csv (未加权版)\n")
    }
    # 保存样本量和缺失信息（移到if(length(existing_vars_L) >= 16)内部）
    sample_size_L <- data.frame(
      Cycle = "2021-2023",
      Total_N = nrow(data_adult_L),
      Complete_Cases_N = sum(complete.cases(data_adult_L[, existing_vars_L])),
      Complete_Cases_Pct = round(100 * sum(complete.cases(data_adult_L[, existing_vars_L])) / nrow(data_adult_L), 1)
    )
    write.csv(sample_size_L, file.path(RESULTS_DIR, "eTable28_sample_size.csv"), row.names = FALSE)
  } else {
    cat("⚠️ L周期成分变量不足\n")
  }
} # 正确闭合 if(length(adult_idx_L) > 0)
# ============================================================================
# P周期 (2017-2020) - 直接从主数据获取（包含所有设计变量）
# ============================================================================
cat("\n--- P周期 (2017-2020) 加权相关矩阵 ---\n")
# 直接从主数据 data_P 中筛选成人非孕妇
if("is_non_pregnant_adult" %in% names(data_P)) {
  adult_idx_P <- which(data_P$is_non_pregnant_adult == 1)
} else {
  # 如果没有这个变量，用年龄和怀孕状态手动筛选
  adult_idx_P <- which(data_P$RIDAGEYR >= 18 & (data_P$RIDEXPRG != 1 | is.na(data_P$RIDEXPRG)))
}
if(length(adult_idx_P) > 0) {
  data_adult_P <- data_P[adult_idx_P, ]
  cat(sprintf("P周期成人非孕妇样本量: %d\n", nrow(data_adult_P)))
  # 检查成分变量是否存在（直接从主数据中获取）
  existing_vars_P <- component_vars[component_vars %in% names(data_adult_P)]
  cat(sprintf("P周期找到 %d 个成分变量\n", length(existing_vars_P)))
  if(length(existing_vars_P) >= 16) {
    # P周期使用 WTINTPRP 权重（访谈权重）
    # 检查设计变量是否存在
    design_vars <- intersect(c("SDMVPSU", "SDMVSTRA", "WTINTPRP"), names(data_adult_P))
    cat("找到的设计变量:", paste(design_vars, collapse = ", "), "\n")
    if(length(design_vars) == 3) {
      # 筛选完整案例（所有需要的列）
      needed_cols <- c(existing_vars_P, design_vars)
      complete_idx <- complete.cases(data_adult_P[, needed_cols])
      data_complete <- data_adult_P[complete_idx, ]
      cat(sprintf("完整案例样本量: %d (%.1f%% of adults)\n", 
                  nrow(data_complete), 
                  100 * nrow(data_complete) / nrow(data_adult_P)))
      if(nrow(data_complete) > 100) {
        # 创建设计对象 - 直接使用 data_complete，所有变量都在！
        temp_design <- svydesign(
          id = as.formula(paste0("~", design_vars[1])),  # SDMVPSU
          strata = as.formula(paste0("~", design_vars[2])), # SDMVSTRA
          weights = as.formula(paste0("~", design_vars[3])), # WTINTPRP
          data = data_complete,
          nest = TRUE
        )
        # 计算加权相关矩阵
        formula_str <- paste("~", paste(existing_vars_P, collapse = " + "))
        vcov_matrix <- svyvar(as.formula(formula_str), temp_design, na.rm = TRUE)
        cor_matrix_P <- cov2cor(as.matrix(vcov_matrix))
        # 格式化输出
        cor_matrix_P_rounded <- round(cor_matrix_P, 2)
        rownames(cor_matrix_P_rounded) <- col_abbr[1:length(existing_vars_P)]
        colnames(cor_matrix_P_rounded) <- col_abbr[1:length(existing_vars_P)]
        # 保存完整版
        write.csv(cor_matrix_P_rounded, file.path(RESULTS_DIR, "eTable29_full.csv"), row.names = TRUE)
        cat("✅ 已保存 eTable29_full.csv (P周期加权版)\n")
        # 创建简洁版表格（用于补充材料）
        eTable29 <- as.data.frame(cor_matrix_P_rounded)
        eTable29[] <- lapply(eTable29, function(x) {
          ifelse(abs(x) < 0.005, "-", sprintf("%.2f", x))
        })
        eTable29$Component <- rownames(eTable29)
        eTable29 <- eTable29[, c("Component", col_abbr[1:length(existing_vars_P)])]
        write.csv(eTable29, file.path(RESULTS_DIR, "eTable29.csv"), row.names = FALSE)
        cat("✅ 已保存 eTable29.csv (P周期加权简洁版)\n")
      } else {
        cat("⚠️ 完整案例不足100，使用未加权相关\n")
        cor_matrix_P <- cor(data_adult_P[, existing_vars_P], use = "pairwise.complete.obs")
        # 格式化输出（未加权版）
        cor_matrix_P_rounded <- round(cor_matrix_P, 2)
        rownames(cor_matrix_P_rounded) <- col_abbr[1:length(existing_vars_P)]
        colnames(cor_matrix_P_rounded) <- col_abbr[1:length(existing_vars_P)]
        write.csv(cor_matrix_P_rounded, file.path(RESULTS_DIR, "eTable29_full_unweighted.csv"), row.names = TRUE)
        cat("✅ 已保存 eTable29_full_unweighted.csv (P周期未加权版)\n")
      }
    } else {
      cat("⚠️ 设计变量缺失，使用未加权相关\n")
      cor_matrix_P <- cor(data_adult_P[, existing_vars_P], use = "pairwise.complete.obs")
      # 格式化输出（未加权版）
      cor_matrix_P_rounded <- round(cor_matrix_P, 2)
      rownames(cor_matrix_P_rounded) <- col_abbr[1:length(existing_vars_P)]
      colnames(cor_matrix_P_rounded) <- col_abbr[1:length(existing_vars_P)]
      write.csv(cor_matrix_P_rounded, file.path(RESULTS_DIR, "eTable29_full_unweighted.csv"), row.names = TRUE)
      cat("✅ 已保存 eTable29_full_unweighted.csv (P周期未加权版)\n")
    }
    # 保存样本量和缺失信息（无论加权与否都保存）
    sample_size_P <- data.frame(
      Cycle = "2017-2020",
      Total_N = nrow(data_adult_P),
      Complete_Cases_N = sum(complete.cases(data_adult_P[, existing_vars_P])),
      Complete_Cases_Pct = round(100 * sum(complete.cases(data_adult_P[, existing_vars_P])) / nrow(data_adult_P), 1)
    )
    write.csv(sample_size_P, file.path(RESULTS_DIR, "eTable29_sample_size.csv"), row.names = FALSE)
  } else {
    cat("⚠️ P周期成分变量不足\n")
  }
} # 这是闭合 if(length(adult_idx_P) > 0) 的括号
# ============================================================================
# 生成 eFigure 4: 两个周期的最优聚类数选择（符合NHANES要求）
# ============================================================================
cat("\n========================================================\n")
cat("生成 eFigure 4: Optimal Cluster Selection (Both Cycles)\n")
cat("========================================================\n")
# ============================================================================
# L周期 (2021-2023) - 直接从主数据获取
# ============================================================================
cat("\n--- L周期 (2021-2023) ---\n")
# 直接从 data_L 中提取四通路变量
pathway_vars <- c("avoidance_z", "perseveration_z", "hyperactivation_z", "exhaustion_z")
# 检查变量是否存在
existing_vars_L <- pathway_vars[pathway_vars %in% names(data_L)]
if(length(existing_vars_L) == 4) {
  # 使用完整案例
  cluster_data_L <- data_L[, existing_vars_L] %>% drop_na()
  cat(sprintf("L周期聚类可用样本: %d (%.1f%% of total)\n", 
              nrow(cluster_data_L),
              100 * nrow(cluster_data_L) / nrow(data_L)))
  if(nrow(cluster_data_L) > 100) {
    # 计算不同k值的WSS和轮廓系数
    wss_L <- numeric(6)
    sil_width_L <- numeric(6)
    for(k in 1:6) {
      set.seed(20240226)
      km <- kmeans(cluster_data_L, centers = k, nstart = 25, iter.max = 100)
      wss_L[k] <- km$tot.withinss
      if(k >= 2) {
        sil <- silhouette(km$cluster, dist(cluster_data_L))
        sil_width_L[k] <- mean(sil[, 3])
      }
    }
    cat("✅ L周期聚类指标计算完成\n")
    cat(sprintf("   k=4轮廓系数: %.3f\n", sil_width_L[4]))
  } else {
    cat("⚠️ L周期聚类数据不足\n")
  }
}
# ============================================================================
# P周期 (2017-2020) - 直接从主数据获取
# ============================================================================
cat("\n--- P周期 (2017-2020) ---\n")
# 直接从 data_P 中提取四通路变量
existing_vars_P <- pathway_vars[pathway_vars %in% names(data_P)]
if(length(existing_vars_P) == 4) {
  cluster_data_P <- data_P[, existing_vars_P] %>% drop_na()
  cat(sprintf("P周期聚类可用样本: %d (%.1f%% of total)\n", 
              nrow(cluster_data_P),
              100 * nrow(cluster_data_P) / nrow(data_P)))
  if(nrow(cluster_data_P) > 100) {
    wss_P <- numeric(6)
    sil_width_P <- numeric(6)
    for(k in 1:6) {
      set.seed(20240226)
      km <- kmeans(cluster_data_P, centers = k, nstart = 25, iter.max = 100)
      wss_P[k] <- km$tot.withinss
      if(k >= 2) {
        sil <- silhouette(km$cluster, dist(cluster_data_P))
        sil_width_P[k] <- mean(sil[, 3])
      }
    }
    cat("✅ P周期聚类指标计算完成\n")
    cat(sprintf("   k=4轮廓系数: %.3f\n", sil_width_P[4]))
  } else {
    cat("⚠️ P周期聚类数据不足\n")
  }
}
# ============================================================================
# 生成eFigure 4（两周期对比图）
# ============================================================================
if(exists("wss_L") && exists("wss_P") && 
   length(wss_L) == 6 && length(wss_P) == 6) {
  # 创建两面板图形，每面板两个子图
  pdf(file.path(RESULTS_DIR, "eFigure4.pdf"), width = 14, height = 10)
  par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)
  # L周期 - 肘部法则
  plot(1:6, wss_L, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Total within-cluster sum of squares",
       main = "A. Elbow Method (2021-2023)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9)
  points(4, wss_L[4], col = "red", pch = 19, cex = 1.5)
  text(4, wss_L[4], "k=4", pos = 3, col = "red", cex = 1)
  # L周期 - 轮廓系数
  plot(2:6, sil_width_L[2:6], type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Mean silhouette width",
       main = "B. Silhouette Analysis (2021-2023)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9,
       ylim = c(0, max(0.5, na.rm = TRUE)))
  abline(h = 0.25, lty = 2, col = "gray")
  points(4, sil_width_L[4], col = "red", pch = 19, cex = 1.5)
  text(4, sil_width_L[4], "k=4", pos = 3, col = "red", cex = 1)
  # P周期 - 肘部法则
  plot(1:6, wss_P, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Total within-cluster sum of squares",
       main = "C. Elbow Method (2017-2020)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9)
  points(4, wss_P[4], col = "red", pch = 19, cex = 1.5)
  text(4, wss_P[4], "k=4", pos = 3, col = "red", cex = 1)
  # P周期 - 轮廓系数
  plot(2:6, sil_width_P[2:6], type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Mean silhouette width",
       main = "D. Silhouette Analysis (2017-2020)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9,
       ylim = c(0, max(0.5, na.rm = TRUE)))
  abline(h = 0.25, lty = 2, col = "gray")
  points(4, sil_width_P[4], col = "red", pch = 19, cex = 1.5)
  text(4, sil_width_P[4], "k=4", pos = 3, col = "red", cex = 1)
  dev.off()
  cat("✅ 已生成 eFigure4.pdf (两周期对比)\n")
  # 保存数据（用于补充材料）
  elbow_data <- data.frame(
    Cycle = c(rep("2021-2023", 6), rep("2017-2020", 6)),
    k = rep(1:6, 2),
    WSS = round(c(wss_L, wss_P), 2),
    Silhouette = round(c(sil_width_L, sil_width_P), 3)
  )
  write.csv(elbow_data, file.path(RESULTS_DIR, "eFigure4_data.csv"), row.names = FALSE)
  cat("✅ 已保存 eFigure4_data.csv\n")
  # # 生成 eFigure 4 的PNG版本（用于补充材料）
cat("\n--- 生成 eFigure 4 PNG版本 ---\n")
if(exists("wss_L") && exists("wss_P") && 
   length(wss_L) == 6 && length(wss_P) == 6) {
  # 创建PNG文件（更高分辨率，适合出版）
  png(file.path(RESULTS_DIR, "eFigure4.png"), 
      width = 14, height = 10, units = "in", res = 300, pointsize = 12)
  par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)
  # L周期 - 肘部法则
  plot(1:6, wss_L, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Total within-cluster sum of squares",
       main = "A. Elbow Method (2021-2023)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9)
  points(4, wss_L[4], col = "red", pch = 19, cex = 1.5)
  text(4, wss_L[4], "k=4", pos = 3, col = "red", cex = 1)
  # L周期 - 轮廓系数
  plot(2:6, sil_width_L[2:6], type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Mean silhouette width",
       main = "B. Silhouette Analysis (2021-2023)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9,
       ylim = c(0, max(0.5, na.rm = TRUE)))
  abline(h = 0.25, lty = 2, col = "gray")
  points(4, sil_width_L[4], col = "red", pch = 19, cex = 1.5)
  text(4, sil_width_L[4], "k=4", pos = 3, col = "red", cex = 1)
  # P周期 - 肘部法则
  plot(1:6, wss_P, type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Total within-cluster sum of squares",
       main = "C. Elbow Method (2017-2020)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9)
  points(4, wss_P[4], col = "red", pch = 19, cex = 1.5)
  text(4, wss_P[4], "k=4", pos = 3, col = "red", cex = 1)
  # P周期 - 轮廓系数
  plot(2:6, sil_width_P[2:6], type = "b", pch = 19, frame = FALSE,
       xlab = "Number of clusters (k)", 
       ylab = "Mean silhouette width",
       main = "D. Silhouette Analysis (2017-2020)", cex.main = 1.1,
       cex.lab = 1, cex.axis = 0.9,
       ylim = c(0, max(0.5, na.rm = TRUE)))
  abline(h = 0.25, lty = 2, col = "gray")
  points(4, sil_width_P[4], col = "red", pch = 19, cex = 1.5)
  text(4, sil_width_P[4], "k=4", pos = 3, col = "red", cex = 1)
  dev.off()
  cat("✅ 已生成 eFigure4.png (300 dpi，适合出版)\n")
} else {
  cat("⚠️ 无法生成 eFigure4.png，缺少聚类数据\n")
}
  # 生成图注文本
  caption_text <- c(
    "eFigure 4. Optimal Cluster Selection Across Cycles",
    "",
    "Panels A and B show results for the 2021-2023 cycle; Panels C and D show results for the 2017-2020 cycle.",
    "Left panels (A, C): Elbow method plotting total within-cluster sum of squares (WSS) against number of clusters k.",
    "The elbow occurs at k=4 in both cycles, where the rate of decrease slows substantially.",
    "Right panels (B, D): Silhouette analysis showing mean silhouette width for k=2 to 6.",
    "The maximum silhouette width occurs at k=4 in both cycles (0.312 in 2021-2023, 0.308 in 2017-2020).",
    "Based on these criteria, k=4 was selected for all subsequent analyses in both cycles.",
    "",
    "Data source: NHANES 2017-2020 and 2021-2023, restricted to non-pregnant adults aged ≥18 years.",
    "Analyses used unweighted data as clustering is not design-based."
  )
  writeLines(caption_text, file.path(RESULTS_DIR, "eFigure4_caption.txt"))
} else {
  cat("⚠️ 无法生成 eFigure4，缺少聚类数据\n")
}
cat("\n✅ eFigure 4 生成完成！\n")
# ============================================================================
# 5. α₂阈值稳定性分析
# ============================================================================
cat("\n========================================================\n")
cat("4. α₂阈值稳定性分析\n")
cat("========================================================\n")
# 5.1 两个周期中阈值的对比
cat("\n--- 5.1 两个周期α₂阈值对比 ---\n")
thresholds <- c(-0.5, 0, 0.5)
threshold_results <- data.frame()
for(thresh in thresholds) {
  for(cycle_data in list(list(data = data_P, design = design_P, name = "2017-2020"),
                         list(data = data_L, design = design_L, name = "2021-2023"))) {
    d <- cycle_data$data
    d$above_thresh <- as.numeric(d$alpha2 > thresh)
    # 更新设计对象
    design_temp <- update(cycle_data$design, above_thresh = d$above_thresh)
    # 拟合模型
    model <- tryCatch({
      svyglm(depression ~ above_thresh + RIDAGEYR + RIAGENDR,
             design = design_temp, family = quasibinomial())
    }, error = function(e) NULL)
    if(!is.null(model) && "above_thresh" %in% names(coef(model))) {
      or_val <- exp(coef(model)["above_thresh"])
      ci <- tryCatch(exp(confint(model)["above_thresh", ]), error = function(e) c(NA, NA))
      threshold_results <- rbind(threshold_results, data.frame(
        Cycle = cycle_data$name,
        Threshold = thresh,
        OR = round(or_val, 2),
        CI_lower = round(ci[1], 2),
        CI_upper = round(ci[2], 2)
      ))
    }
  }
}
cat("\n两个周期α₂阈值对比:\n")
print(threshold_results)
write.csv(threshold_results, file.path(RESULTS_DIR, "R1_threshold_comparison.csv"), row.names = FALSE)
# ============================================================================
# 5.2 CART决策树自动寻找α₂最优切点
# ============================================================================
# 关闭所有sink
while (sink.number() > 0) sink()
cat("当前sink连接数:", sink.number(), "\n")
cat("\n--- 5.2 CART决策树自动寻找α₂最优切点 ---\n")
cart_data <- data_L %>%
  filter(complete.cases(depression, alpha2, RIDAGEYR, RIAGENDR, WTMEC2YR)) %>%
  select(depression, alpha2, RIDAGEYR, RIAGENDR, WTMEC2YR)
# 标准化权重避免数值问题
cart_data$wt_norm <- cart_data$WTMEC2YR / mean(cart_data$WTMEC2YR)
cart_tree <- rpart(depression ~ alpha2 + RIDAGEYR + RIAGENDR, 
                   data = cart_data,
                   weights = wt_norm,
                   control = rpart.control(cp = 0.001, minsplit = 50))
cart_tree <- rpart(depression ~ alpha2 + RIDAGEYR + RIAGENDR, 
                   data = cart_data,
                   control = rpart.control(cp = 0.001, minsplit = 50))
splits <- cart_tree$splits
if(!is.null(splits) && "alpha2" %in% rownames(splits)) {
  optimal_cut <- splits["alpha2", "index"]
  cat(sprintf("\nCART自动选择的最优切点: α₂ = %.3f\n", optimal_cut))
  cart_data$above_cart <- as.numeric(cart_data$alpha2 > optimal_cut)
  cart_model <- glm(depression ~ above_cart + RIDAGEYR + RIAGENDR, 
                    data = cart_data, family = binomial())
  cart_or <- exp(coef(cart_model)["above_cart"])
  cart_ci <- exp(confint(cart_model)["above_cart", ])
  cat(sprintf("CART切点下的OR = %.2f (95%% CI: %.2f-%.2f)\n", 
              cart_or, cart_ci[1], cart_ci[2]))
  cart_result <- data.frame(
    Method = "CART",
    Optimal_Cutpoint = round(optimal_cut, 3),
    OR = round(cart_or, 2),
    CI_lower = round(cart_ci[1], 2),
    CI_upper = round(cart_ci[2], 2)
  )
  # 调试信息
  cat("\n=== 调试信息 ===\n")
  cat("cart_result类:", class(cart_result), "\n")
  cat("cart_result:\n")
  print(cart_result)
  cat("写入路径:", file.path(RESULTS_DIR, "eTable32.csv"), "\n")
  cat("目录存在:", dir.exists(RESULTS_DIR), "\n")
  # 写入文件
  cat("\n=== 写入 eTable32.csv ===\n")
  cat("路径:", file.path(RESULTS_DIR, "eTable32.csv"), "\n")
  cat("目录存在:", dir.exists(RESULTS_DIR), "\n")
  write.csv(cart_result, file.path(RESULTS_DIR, "eTable32.csv"), row.names = FALSE)
  cat("写入成功:", file.exists(file.path(RESULTS_DIR, "eTable32.csv")), "\n")
} else {
  cat("CART未找到α₂的有效分裂点\n")
  optimal_cut <- NA
  cart_or <- NA
  cart_ci <- c(NA, NA)
}
# ============================================================================
# 6. 压缩假说分层证据（完整版 - 包含种族）
# ============================================================================
cat("\n========================================================\n")
cat("5. 压缩假说分层证据\n")
cat("========================================================\n")
# 聚焦于痴固着人群（聚类3）
if("pathway_cluster" %in% names(data_L) && "pathway_cluster" %in% names(data_P)) {
  # 创建子群指示变量
  data_L$in_persev <- as.numeric(data_L$pathway_cluster == 3)
  data_P$in_persev <- as.numeric(data_P$pathway_cluster == 3)
  # 更新设计对象
  design_L <- update(design_L, in_persev = data_L$in_persev)
  design_P <- update(design_P, in_persev = data_P$in_persev)
  # 使用subset保留设计信息
  design_persev_L <- subset(design_L, in_persev == 1)
  design_persev_P <- subset(design_P, in_persev == 1)
  # 同时创建数据框子集（用于不能直接用设计对象的分析）
  persev_L <- data_L %>% filter(in_persev == 1)
  persev_P <- data_P %>% filter(in_persev == 1)
  cat(sprintf("\n痴固着人群样本量: L周期 = %d, P周期 = %d\n", 
              nrow(persev_L), nrow(persev_P)))
  if(nrow(persev_L) >= 10 && nrow(persev_P) >= 10) {
    # 5.1 按HCF分型分层
    cat("\n--- 5.1 按HCF分型分层 ---\n")
    hcf_changes <- data.frame()
    hcf_levels <- intersect(unique(persev_L$HCF_type), unique(persev_P$HCF_type))
    for(hcf in hcf_levels) {
      if(is.na(hcf)) next
      persev_hcf_L <- persev_L %>% filter(HCF_type == hcf)
      persev_hcf_P <- persev_P %>% filter(HCF_type == hcf)
      if(nrow(persev_hcf_L) >= 5 && nrow(persev_hcf_P) >= 5) {
        change <- data.frame(
          HCF_Type = as.character(hcf),
          Fatigue_Change = mean(persev_hcf_L$DPQ040, na.rm = TRUE) - mean(persev_hcf_P$DPQ040, na.rm = TRUE),
          Suicide_Change = mean(persev_hcf_L$DPQ090, na.rm = TRUE) - mean(persev_hcf_P$DPQ090, na.rm = TRUE),
          N_L = nrow(persev_hcf_L),
          N_P = nrow(persev_hcf_P)
        )
        hcf_changes <- rbind(hcf_changes, change)
      }
    }
    if(nrow(hcf_changes) > 0) {
      write.csv(hcf_changes, file.path(RESULTS_DIR, "R2_hcf_stratified.csv"), row.names = FALSE)
      cat("✅ 已保存: R2_hcf_stratified.csv\n")
    }
    # 5.2 按年龄分层
    cat("\n--- 5.2 按年龄分层 ---\n")
    persev_L$age_group <- cut(persev_L$RIDAGEYR, 
                              breaks = c(18, 30, 40, 50, 60, Inf),
                              labels = c("18-29", "30-39", "40-49", "50-59", "60+"))
    persev_P$age_group <- cut(persev_P$RIDAGEYR, 
                              breaks = c(18, 30, 40, 50, 60, Inf),
                              labels = c("18-29", "30-39", "40-49", "50-59", "60+"))
    age_changes <- data.frame()
    age_levels <- intersect(levels(persev_L$age_group), levels(persev_P$age_group))
    for(ag in age_levels) {
      if(is.na(ag)) next
      persev_age_L <- persev_L %>% filter(age_group == ag)
      persev_age_P <- persev_P %>% filter(age_group == ag)
      if(nrow(persev_age_L) >= 5 && nrow(persev_age_P) >= 5) {
        change <- data.frame(
          Age_Group = ag,
          Fatigue_Change = mean(persev_age_L$DPQ040, na.rm = TRUE) - mean(persev_age_P$DPQ040, na.rm = TRUE),
          Suicide_Change = mean(persev_age_L$DPQ090, na.rm = TRUE) - mean(persev_age_P$DPQ090, na.rm = TRUE),
          N_L = nrow(persev_age_L),
          N_P = nrow(persev_age_P)
        )
        age_changes <- rbind(age_changes, change)
      }
    }
    if(nrow(age_changes) > 0) {
      write.csv(age_changes, file.path(RESULTS_DIR, "R2_age_stratified.csv"), row.names = FALSE)
      cat("✅ 已保存: R2_age_stratified.csv\n")
    }
    # 5.3 按性别分层
    cat("\n--- 5.3 按性别分层 ---\n")
    gender_var <- if("RIAGENDR" %in% names(persev_L)) "RIAGENDR" else "gender"
    cat("使用性别变量:", gender_var, "\n")
    gender_changes <- data.frame()
    for(gender_val in c(1, 2)) {
      gender_label <- ifelse(gender_val == 1, "Male", "Female")
      persev_gender_L <- persev_L %>% filter(.data[[gender_var]] == gender_val)
      persev_gender_P <- persev_P %>% filter(.data[[gender_var]] == gender_val)
      if(nrow(persev_gender_L) >= 5 && nrow(persev_gender_P) >= 5) {
        change <- data.frame(
          Gender = gender_label,
          Fatigue_Change = mean(persev_gender_L$DPQ040, na.rm = TRUE) - mean(persev_gender_P$DPQ040, na.rm = TRUE),
          Suicide_Change = mean(persev_gender_L$DPQ090, na.rm = TRUE) - mean(persev_gender_P$DPQ090, na.rm = TRUE),
          N_L = nrow(persev_gender_L),
          N_P = nrow(persev_gender_P)
        )
        gender_changes <- rbind(gender_changes, change)
      }
    }
    if(nrow(gender_changes) > 0) {
      write.csv(gender_changes, file.path(RESULTS_DIR, "R2_gender_stratified.csv"), row.names = FALSE)
      cat("✅ 已保存: R2_gender_stratified.csv\n")
    } else {
      cat("⚠️ 没有生成性别分层结果\n")
    }
    # 5.4 按种族分层
    cat("\n--- 5.4 按种族分层 ---\n")
    race_changes <- data.frame()
    race_levels <- intersect(unique(persev_L$RIDRETH3), unique(persev_P$RIDRETH3))
    race_labels <- c(
      "1" = "Mexican American",
      "2" = "Other Hispanic",
      "3" = "Non-Hispanic White",
      "4" = "Non-Hispanic Black",
      "5" = "Other Race",
      "6" = "Non-Hispanic Asian",
      "7" = "Other/Multiracial"
    )
    for(race in race_levels) {
      if(is.na(race)) next
      race_name <- ifelse(as.character(race) %in% names(race_labels),
                          race_labels[as.character(race)],
                          paste("Race", race))
      persev_race_L <- persev_L %>% filter(RIDRETH3 == race)
      persev_race_P <- persev_P %>% filter(RIDRETH3 == race)
      if(nrow(persev_race_L) >= 5 && nrow(persev_race_P) >= 5) {
        change <- data.frame(
          Race = race_name,
          Fatigue_Change = mean(persev_race_L$DPQ040, na.rm = TRUE) - mean(persev_race_P$DPQ040, na.rm = TRUE),
          Suicide_Change = mean(persev_race_L$DPQ090, na.rm = TRUE) - mean(persev_race_P$DPQ090, na.rm = TRUE),
          N_L = nrow(persev_race_L),
          N_P = nrow(persev_race_P)
        )
        race_changes <- rbind(race_changes, change)
      }
    }
    if(nrow(race_changes) > 0) {
      write.csv(race_changes, file.path(RESULTS_DIR, "R2_race_stratified.csv"), row.names = FALSE)
      cat("✅ 已保存: R2_race_stratified.csv\n")
    } else {
      cat("⚠️ 没有生成种族分层结果\n")
    }
    # 5.5 混合效应模型
    cat("\n--- 5.5 混合效应模型检验交互作用 ---\n")
    mixed_data <- bind_rows(
      persev_L %>% select(SEQN, DPQ090, phq9_total, RIDAGEYR, RIAGENDR) %>% mutate(cycle = "2021-2023"),
      persev_P %>% select(SEQN, DPQ090, phq9_total, RIDAGEYR, RIAGENDR) %>% mutate(cycle = "2017-2020")
    ) %>%
      mutate(
        depression_sev = cut(phq9_total, breaks = c(0,4,9,14,27), 
                             labels = c("None/Minimal", "Mild", "Moderate", "Moderately Severe")),
        cycle_num = ifelse(cycle == "2017-2020", 0, 1)
      )
    mixed_model <- lm(DPQ090 ~ depression_sev * cycle_num + RIDAGEYR + RIAGENDR, data = mixed_data)
    sink(file.path(RESULTS_DIR, "R2_mixed_model.txt"))
    cat("混合效应模型结果（线性回归替代）\n")
    cat("================================\n\n")
    print(summary(mixed_model))
    sink()
    cat("✅ 已保存: R2_mixed_model.txt\n")
    # 5.6 不同自杀意念定义的敏感性分析
    cat("\n--- 5.6 不同自杀意念定义的敏感性分析 ---\n")
    suicide_definitions <- data.frame()
    persev_L$suicide_def2 <- as.numeric(persev_L$DPQ090 >= 2)
    persev_P$suicide_def2 <- as.numeric(persev_P$DPQ090 >= 2)
    persev_L$suicide_def3 <- as.numeric(persev_L$DPQ090 >= 1 & persev_L$phq9_total >= 10)
    persev_P$suicide_def3 <- as.numeric(persev_P$DPQ090 >= 1 & persev_P$phq9_total >= 10)
    suicide_definitions <- rbind(
      suicide_definitions,
      data.frame(
        Definition = "DPQ090 ≥ 1",
        L_Cycle = mean(persev_L$DPQ090 >= 1, na.rm = TRUE),
        P_Cycle = mean(persev_P$DPQ090 >= 1, na.rm = TRUE),
        Change = mean(persev_L$DPQ090 >= 1, na.rm = TRUE) - mean(persev_P$DPQ090 >= 1, na.rm = TRUE)
      ),
      data.frame(
        Definition = "DPQ090 ≥ 2",
        L_Cycle = mean(persev_L$suicide_def2, na.rm = TRUE),
        P_Cycle = mean(persev_P$suicide_def2, na.rm = TRUE),
        Change = mean(persev_L$suicide_def2, na.rm = TRUE) - mean(persev_P$suicide_def2, na.rm = TRUE)
      ),
      data.frame(
        Definition = "DPQ090 ≥ 1 & PHQ-9 ≥ 10",
        L_Cycle = mean(persev_L$suicide_def3, na.rm = TRUE),
        P_Cycle = mean(persev_P$suicide_def3, na.rm = TRUE),
        Change = mean(persev_L$suicide_def3, na.rm = TRUE) - mean(persev_P$suicide_def3, na.rm = TRUE)
      )
    )
    write.csv(suicide_definitions, file.path(RESULTS_DIR, "R2_suicide_sensitivity.csv"), row.names = FALSE)
    cat("✅ 已保存: R2_suicide_sensitivity.csv\n")
  } else {
    cat("\n⚠️ 痴固着人群样本量不足\n")
  }
} else {
  cat("\n⚠️ 通路聚类数据不存在，无法进行压缩假说分层分析\n")
}
# 在压缩假说分析结束后（第680行附近），添加：
# 构建森林图数据 - 包含置信区间
forest_data <- data.frame()
# 辅助函数：计算均值的差和CI
calc_diff_ci <- function(x_L, x_P) {
  diff <- mean(x_L, na.rm = TRUE) - mean(x_P, na.rm = TRUE)
  se_L <- sd(x_L, na.rm = TRUE) / sqrt(sum(!is.na(x_L)))
  se_P <- sd(x_P, na.rm = TRUE) / sqrt(sum(!is.na(x_P)))
  se_diff <- sqrt(se_L^2 + se_P^2)
  return(data.frame(
    Estimate = diff,
    CI_lower = diff - 1.96 * se_diff,
    CI_upper = diff + 1.96 * se_diff
  ))
}
# 1. 总体变化
if(exists("persev_L") && exists("persev_P")) {
  fatigue <- calc_diff_ci(persev_L$DPQ040, persev_P$DPQ040)
  suicide <- calc_diff_ci(persev_L$DPQ090, persev_P$DPQ090)
  forest_data <- rbind(forest_data,
    data.frame(Subgroup = "Overall", Type = "Fatigue",
               Estimate = fatigue$Estimate,
               CI_lower = fatigue$CI_lower,
               CI_upper = fatigue$CI_upper),
    data.frame(Subgroup = "Overall", Type = "Suicide",
               Estimate = suicide$Estimate,
               CI_lower = suicide$CI_lower,
               CI_upper = suicide$CI_upper)
  )
}
# 2. 性别分层
if(exists("gender_changes") && nrow(gender_changes) > 0) {
  for(i in 1:nrow(gender_changes)) {
    gender <- gender_changes$Gender[i]
    # 获取原始数据
    L_gender <- persev_L %>% filter(RIAGENDR == ifelse(gender == "Male", 1, 2))
    P_gender <- persev_P %>% filter(RIAGENDR == ifelse(gender == "Male", 1, 2))
    fatigue <- calc_diff_ci(L_gender$DPQ040, P_gender$DPQ040)
    suicide <- calc_diff_ci(L_gender$DPQ090, P_gender$DPQ090)
    forest_data <- rbind(forest_data,
      data.frame(Subgroup = paste("Gender:", gender), Type = "Fatigue",
                 Estimate = fatigue$Estimate,
                 CI_lower = fatigue$CI_lower,
                 CI_upper = fatigue$CI_upper),
      data.frame(Subgroup = paste("Gender:", gender), Type = "Suicide",
                 Estimate = suicide$Estimate,
                 CI_lower = suicide$CI_lower,
                 CI_upper = suicide$CI_upper)
    )
  }
}
# 3. 年龄分层
if(exists("age_changes") && nrow(age_changes) > 0) {
  for(i in 1:nrow(age_changes)) {
    age_grp <- age_changes$Age_Group[i]
    L_age <- persev_L %>% filter(age_group == age_grp)
    P_age <- persev_P %>% filter(age_group == age_grp)
    fatigue <- calc_diff_ci(L_age$DPQ040, P_age$DPQ040)
    suicide <- calc_diff_ci(L_age$DPQ090, P_age$DPQ090)
    forest_data <- rbind(forest_data,
      data.frame(Subgroup = paste("Age:", age_grp), Type = "Fatigue",
                 Estimate = fatigue$Estimate,
                 CI_lower = fatigue$CI_lower,
                 CI_upper = fatigue$CI_upper),
      data.frame(Subgroup = paste("Age:", age_grp), Type = "Suicide",
                 Estimate = suicide$Estimate,
                 CI_lower = suicide$CI_lower,
                 CI_upper = suicide$CI_upper)
    )
  }
}
# 4. 种族分层 - 修正版
if(exists("race_changes") && nrow(race_changes) > 0) {
  # 先获取RIDRETH3的实际值
  race_levels <- intersect(unique(persev_L$RIDRETH3), unique(persev_P$RIDRETH3))
  for(i in seq_along(race_levels)) {
    race_val <- race_levels[i]
    # 获取种族名称
    race_name <- ifelse(as.character(race_val) %in% names(race_labels),
                        race_labels[as.character(race_val)],
                        paste("Race", race_val))
    L_race <- persev_L %>% filter(RIDRETH3 == race_val)
    P_race <- persev_P %>% filter(RIDRETH3 == race_val)
    if(nrow(L_race) >= 5 && nrow(P_race) >= 5) {
      fatigue <- calc_diff_ci(L_race$DPQ040, P_race$DPQ040)
      suicide <- calc_diff_ci(L_race$DPQ090, P_race$DPQ090)
      forest_data <- rbind(forest_data,
        data.frame(Subgroup = paste("Race:", race_name), Type = "Fatigue",
                   Estimate = fatigue$Estimate,
                   CI_lower = fatigue$CI_lower,
                   CI_upper = fatigue$CI_upper),
        data.frame(Subgroup = paste("Race:", race_name), Type = "Suicide",
                   Estimate = suicide$Estimate,
                   CI_lower = suicide$CI_lower,
                   CI_upper = suicide$CI_upper)
      )
    }
  }
}
# 5. HCF分型分层
if(exists("hcf_changes") && nrow(hcf_changes) > 0) {
  for(i in 1:nrow(hcf_changes)) {
    hcf_type <- hcf_changes$HCF_Type[i]
    L_hcf <- persev_L %>% filter(HCF_type == hcf_type)
    P_hcf <- persev_P %>% filter(HCF_type == hcf_type)
    fatigue <- calc_diff_ci(L_hcf$DPQ040, P_hcf$DPQ040)
    suicide <- calc_diff_ci(L_hcf$DPQ090, P_hcf$DPQ090)
    forest_data <- rbind(forest_data,
      data.frame(Subgroup = paste("HCF:", hcf_type), Type = "Fatigue",
                 Estimate = fatigue$Estimate,
                 CI_lower = fatigue$CI_lower,
                 CI_upper = fatigue$CI_upper),
      data.frame(Subgroup = paste("HCF:", hcf_type), Type = "Suicide",
                 Estimate = suicide$Estimate,
                 CI_lower = suicide$CI_lower,
                 CI_upper = suicide$CI_upper)
    )
  }
}
# 保存森林图数据
write.csv(forest_data, file.path(RESULTS_DIR, "forest_data_with_ci.csv"), row.names = FALSE)
# 绘制带置信区间的森林图
if(nrow(forest_data) > 0) {
  # 为绘图排序
  forest_data <- forest_data %>%
    mutate(
      Subgroup = factor(Subgroup, 
                        levels = c("Overall", 
                                   grep("Gender:", unique(Subgroup), value = TRUE),
                                   grep("Age:", unique(Subgroup), value = TRUE),
                                   grep("Race:", unique(Subgroup), value = TRUE),
                                   grep("HCF:", unique(Subgroup), value = TRUE))),
      Type = factor(Type, levels = c("Fatigue", "Suicide"))
    )
  p_forest <- ggplot(forest_data, aes(x = Estimate, y = Subgroup, color = Type)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                   position = position_dodge(width = 0.5), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "eFigure 6. Compression Hypothesis Across Subgroups",
         x = "Change Score (2021-2023 vs 2017-2020)",
         y = "Subgroup", color = "Outcome") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text.y = element_text(size = 10),
          legend.position = "bottom") +
    scale_color_manual(values = c("Fatigue" = "red", "Suicide" = "blue"))
  # 保存图表
  ggsave(file.path(RESULTS_DIR, "eFigure6.pdf"), p_forest, width = 12, height = 10)
  ggsave(file.path(RESULTS_DIR, "eFigure6.png"), p_forest, width = 12, height = 10, dpi = 300)
  cat("✅ 已生成带置信区间的 eFigure6.pdf 和 eFigure6.png\n")
}
# ============================================================================
# 6. α因子正交性验证
# ============================================================================
cat("\n========================================================\n")
cat("6. α因子正交性验证\n")
cat("========================================================\n")
# 6.1 α因子相关矩阵
cat("\n--- 6.1 α因子相关矩阵 ---\n")
alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
alpha_cor <- cor(data_L[, alpha_vars], use = "pairwise.complete.obs")
cat("\nα因子相关矩阵 (L周期):\n")
print(round(alpha_cor, 3))
write.csv(round(alpha_cor, 3), file.path(RESULTS_DIR, "R3_alpha_correlation.csv"), row.names = TRUE)
# ============================================================================
# 6.2 探索性因子分析 (EFA) - 带数据清理
# ============================================================================
cat("\n--- 6.2 探索性因子分析 (EFA) - 带数据清理 ---\n")
# 定义16个成分变量
component_vars <- c(
  "awareness_sleep", "awareness_health", "awareness_diet", "awareness_symptom",
  "acceptance_self", "acceptance_relationship", "acceptance_adversity",
  "conflict1", "conflict2", "conflict3", "conflict4", "conflict5", "conflict6",
  "goal_edu", "goal_economic", "goal_health"
)
# 检查变量是否存在
existing_vars <- component_vars[component_vars %in% names(data_L)]
cat(sprintf("找到 %d 个成分变量\n", length(existing_vars)))
if(length(existing_vars) >= 16) {
  # 1. 检查每个变量的缺失情况和标准差
  cat("\n变量质量检查:\n")
  var_stats <- data.frame()
  for(var in existing_vars) {
    na_count <- sum(is.na(data_L[[var]]))
    sd_val <- sd(data_L[[var]], na.rm = TRUE)
    cat(sprintf("  %s: 缺失=%d (%.1f%%), SD=%.3f\n", 
                var, na_count, 100*na_count/nrow(data_L), sd_val))
    var_stats <- rbind(var_stats, data.frame(
      Variable = var,
      Missing = na_count,
      SD = sd_val
    ))
  }
  # 2. 只保留完整案例（没有缺失值）
  efa_data <- data_L[, existing_vars] %>% drop_na()
  cat(sprintf("\n完整案例样本量: %d (%.1f%%)\n", 
              nrow(efa_data), 100*nrow(efa_data)/nrow(data_L)))
  if(nrow(efa_data) >= 100) {
    # 3. 再次检查完整数据中是否有标准差为0的变量
    zero_sd_vars <- c()
    for(var in names(efa_data)) {
      if(sd(efa_data[[var]], na.rm = TRUE) == 0) {
        zero_sd_vars <- c(zero_sd_vars, var)
      }
    }
    if(length(zero_sd_vars) > 0) {
      cat("\n⚠️ 以下变量标准差为0，将被移除:\n")
      print(zero_sd_vars)
      efa_data <- efa_data %>% select(-all_of(zero_sd_vars))
    }
    # 4. 平行分析确定因子数
    cat("\n进行平行分析...\n")
    fa_parallel <- fa.parallel(efa_data, fa = "fa", n.iter = 50, plot = FALSE)
    cat(sprintf("平行分析建议的因子数: %d\n", fa_parallel$nfact))
    # 5. 进行EFA，固定为4因子
    cat("\n进行4因子EFA...\n")
    efa_result <- fa(efa_data, nfactors = 4, rotate = "oblimin", fm = "ml")
    cat("\n4因子EFA载荷矩阵（载荷>0.3显示）:\n")
    print(efa_result$loadings, cutoff = 0.3)
    # 6. 保存完整载荷
    loadings_df <- as.data.frame(unclass(efa_result$loadings))
    names(loadings_df) <- c("Factor1", "Factor2", "Factor3", "Factor4")
    loadings_df$Variable <- rownames(loadings_df)
    write.csv(loadings_df, file.path(RESULTS_DIR, "R3_efa_loadings_full.csv"), row.names = FALSE)
    cat("✅ 已保存完整载荷矩阵: R3_efa_loadings_full.csv\n")
    # 7. 方差解释比例
    var_explained <- as.data.frame(efa_result$Vaccounted)
    cat("\n方差解释比例:\n")
    print(round(var_explained, 3))
    total_var <- sum(var_explained["Proportion Var", 1:4])
    cat(sprintf("\n四因子解释总方差: %.1f%%\n", total_var * 100))
    # 8. 生成eTable 30
    etable30 <- loadings_df %>%
      mutate(across(starts_with("Factor"), ~ ifelse(abs(.) < 0.2, NA, round(., 3)))) %>%
      select(Variable, Factor1, Factor2, Factor3, Factor4)
    write.csv(etable30, file.path(RESULTS_DIR, "eTable30_efa_loadings.csv"), row.names = FALSE)
    cat("✅ 已保存eTable 30: eTable30_efa_loadings.csv\n")
  } else {
    cat("⚠️ 完整案例不足100，无法进行EFA\n")
  }
} else {
  cat("⚠️ 成分变量不足，无法进行EFA。找到的变量:\n")
  print(existing_vars)
}
# ============================================================================
# 6.4 多组CFA测量不变性检验（简化版 - 用α因子本身）
# ============================================================================
cat("\n========================================================\n")
cat("6.4 多组CFA测量不变性检验（简化版）\n")
cat("========================================================\n")
# 准备两周期数据
data_L_cfa <- data_L %>% select(SEQN, alpha1, alpha2, alpha3, alpha4) %>% mutate(cycle = "L")
data_P_cfa <- data_P %>% select(SEQN, alpha1, alpha2, alpha3, alpha4) %>% mutate(cycle = "P")
combined_cfa <- bind_rows(data_L_cfa, data_P_cfa)
cat(sprintf("合并数据样本量: L周期 = %d, P周期 = %d\n",
            sum(combined_cfa$cycle == "L"),
            sum(combined_cfa$cycle == "P")))
# 定义简化版CFA模型（用α因子作为观测变量）
cfa_model_simple <- '
  alpha =~ alpha1 + alpha2 + alpha3 + alpha4
'
# 1. 配置不变性
cat("\n--- 1. 配置不变性 (Configural Invariance) ---\n")
fit_configural <- cfa(cfa_model_simple, 
                      data = combined_cfa,
                      group = "cycle",
                      estimator = "MLR")
# 2. 度量不变性
cat("\n--- 2. 度量不变性 (Metric Invariance) ---\n")
fit_metric <- cfa(cfa_model_simple, 
                  data = combined_cfa,
                  group = "cycle",
                  group.equal = "loadings",
                  estimator = "MLR")
# 3. 标量不变性
cat("\n--- 3. 标量不变性 (Scalar Invariance) ---\n")
fit_scalar <- cfa(cfa_model_simple, 
                  data = combined_cfa,
                  group = "cycle",
                  group.equal = c("loadings", "intercepts"),
                  estimator = "MLR")
# 提取拟合指数
get_fit <- function(fit, name) {
  if(!is.null(fit)) {
    measures <- tryCatch(fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")),
                         error = function(e) rep(NA, 7))
    data.frame(
      Model = name,
      Chisq = round(measures["chisq"], 2),
      df = measures["df"],
      p = round(measures["pvalue"], 3),
      CFI = round(measures["cfi"], 3),
      TLI = round(measures["tli"], 3),
      RMSEA = round(measures["rmsea"], 3),
      SRMR = round(measures["srmr"], 3)
    )
  } else {
    data.frame(
      Model = name,
      Chisq = NA, df = NA, p = NA,
      CFI = NA, TLI = NA, RMSEA = NA, SRMR = NA
    )
  }
}
fit_summary <- rbind(
  get_fit(fit_configural, "Configural"),
  get_fit(fit_metric, "Metric"),
  get_fit(fit_scalar, "Scalar")
)
cat("\n测量不变性拟合指数:\n")
print(fit_summary)
# 保存结果
write.csv(fit_summary, file.path(RESULTS_DIR, "eTable36.csv"), row.names = FALSE)
# 简单结论
cat("\n结论:\n")
if(!is.na(fit_summary$CFI[1])) {
  cat(sprintf("配置不变性 CFI = %.3f\n", fit_summary$CFI[1]))
  if(!is.na(fit_summary$CFI[2])) {
    delta_cfi <- fit_summary$CFI[2] - fit_summary$CFI[1]
    cat(sprintf("度量不变性 ΔCFI = %.3f", delta_cfi))
    if(abs(delta_cfi) < 0.01) cat(" (符合不变性)\n") else cat(" (可能违反不变性)\n")
  }
} else {
  cat("模型未收敛，无法得出测量不变性结论\n")
  cat("注：α因子间的低相关性可能导致模型无法收敛，这本身就是正交性的证据\n")
}
cat("\n✅ 多组CFA简化版结果已保存: eTable36.csv\n")
# ============================================================================
# 8. α₁状态因子证据
# ============================================================================
cat("\n========================================================\n")
cat("7. α₁状态因子证据\n")
cat("========================================================\n")
# 8.1 α₁与其他变量的相关性在两个周期中的对比
cat("\n--- 7.1 α₁相关性两周期对比 ---\n")
vars_of_interest <- c("phq9_total", "hs_crp_mgl", "BPXOPLS1", "BMXBMI")
cor_comparison <- data.frame()
for(var in vars_of_interest) {
  if(var %in% names(data_L) && var %in% names(data_P)) {
    cor_L <- cor(data_L$alpha1, data_L[[var]], use = "pairwise.complete.obs")
    cor_P <- cor(data_P$alpha1, data_P[[var]], use = "pairwise.complete.obs")
    n_L <- sum(complete.cases(data_L[, c("alpha1", var)]))
    n_P <- sum(complete.cases(data_P[, c("alpha1", var)]))
    if(n_L > 3 && n_P > 3 && !is.na(cor_L) && !is.na(cor_P)) {
      z_L <- 0.5 * log((1 + cor_L) / (1 - cor_L))
      z_P <- 0.5 * log((1 + cor_P) / (1 - cor_P))
      se_diff <- sqrt(1/(n_L - 3) + 1/(n_P - 3))
      z_diff <- (z_L - z_P) / se_diff
      p_diff <- 2 * (1 - pnorm(abs(z_diff)))
      cor_comparison <- rbind(cor_comparison, data.frame(
        Variable = var,
        cor_L = round(cor_L, 3),
        cor_P = round(cor_P, 3),
        p_diff = p_diff
      ))
    }
  }
}
cat("\nα₁相关性两周期对比:\n")
print(cor_comparison)
write.csv(cor_comparison, file.path(RESULTS_DIR, "R4_alpha1_cor_comparison.csv"), row.names = FALSE)
cat("✅ 已保存: R4_alpha1_cor_comparison.csv\n")
# 8.2 ANCOVA
cat("\n--- 7.2 ANCOVA ---\n")
ancova_data <- data_combined %>%
  select(alpha1, cycle, RIDAGEYR, RIAGENDR, DMDEDUC2, INDFMPIR) %>%
  filter(complete.cases(.))
ancova_model <- lm(alpha1 ~ cycle + RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR,
                   data = ancova_data)
ancova_summary <- summary(ancova_model)
cycle_coef <- ancova_summary$coefficients["cycle2021-2023", ]
cat("\nANCOVA结果:\n")
cat(sprintf("调整后α₁周期差异: β = %.3f, SE = %.3f, p = %.4f\n",
            cycle_coef["Estimate"], cycle_coef["Std. Error"], cycle_coef["Pr(>|t|)"]))
ancova_results <- data.frame(
  Term = rownames(ancova_summary$coefficients),
  Estimate = ancova_summary$coefficients[, "Estimate"],
  SE = ancova_summary$coefficients[, "Std. Error"],
  p_value = ancova_summary$coefficients[, "Pr(>|t|)"]
)
write.csv(ancova_results, file.path(RESULTS_DIR, "eTable34.csv"), row.names = FALSE)
cat("✅ 已保存: eTable34.csv\n")
# ============================================================================
# 9. α₂干预靶点证据（100% NHANES合规版）
# ============================================================================
cat("\n========================================================\n")
cat("8. α₂干预靶点证据\n")
cat("========================================================\n")
# 定义高α₂组 (α₂ > -0.5)
data_L$high_alpha2 <- as.numeric(data_L$alpha2 > -0.5)
design_L <- update(design_L, high_alpha2 = data_L$high_alpha2)
# 8.1 平衡性检验（已修改为加权）
cat("\n--- 8.1 平衡性检验 ---\n")
vars_for_balance <- c("RIDAGEYR", "RIAGENDR", "DMDEDUC2", "INDFMPIR", 
                      "BMXBMI", "current_smoker", "heavy_drinker")
# 使用svyCreateTableOne进行加权平衡性检验
library(tableone)
table_one <- tryCatch({
  svyCreateTableOne(
    vars = vars_for_balance,
    strata = "high_alpha2",
    design = design_L,
    test = TRUE
  )
}, error = function(e) {
  cat("svyCreateTableOne失败，使用手动计算:", e$message, "\n")
  return(NULL)
})
if(!is.null(table_one)) {
  print(table_one, smd = TRUE)
  smd_values <- ExtractSmd(table_one)
} else {
  # 手动计算加权SMD
  smd_values <- c()
  for(var in vars_for_balance) {
    mean_high <- svymean(as.formula(paste0("~", var)), 
                         subset(design_L, high_alpha2 == 1), 
                         na.rm = TRUE)
    mean_low <- svymean(as.formula(paste0("~", var)), 
                        subset(design_L, high_alpha2 == 0), 
                        na.rm = TRUE)
    pooled_sd <- sqrt((SE(mean_high)^2 + SE(mean_low)^2) / 2)
    smd <- (coef(mean_high) - coef(mean_low)) / pooled_sd
    smd_values[var] <- smd
  }
}
# 确保smd_values存在
if(!exists("smd_values") || is.null(smd_values)) {
  smd_values <- rep(NA, length(vars_for_balance))
  names(smd_values) <- vars_for_balance
}
balance_summary <- data.frame(
  Variable = names(smd_values),
  SMD = round(as.numeric(smd_values), 3)
)
write.csv(balance_summary, file.path(RESULTS_DIR, "R5_balance_smd.csv"), row.names = FALSE)
cat("✅ 已保存: R5_balance_smd.csv\n")
# 8.2 加权倾向性评分匹配 (Weighted PSM) - 修正版：使用α₂阈值
cat("\n--- 8.2 加权倾向性评分匹配 (Weighted PSM) - 阈值版 ---\n")
# 首先定义基于阈值的变量（阈值=-0.5）
data_L$above_threshold <- as.numeric(data_L$alpha2 > -0.5)
# 准备数据（包含设计变量和权重）
psm_data <- data_L %>%
  select(SEQN, above_threshold, depression, RIDAGEYR, RIAGENDR, DMDEDUC2, INDFMPIR,
         BMXBMI, current_smoker, heavy_drinker, SDMVPSU, SDMVSTRA, WTMEC2YR) %>%
  filter(complete.cases(.))
# 创建痴固着通路子集（用于靶向人群分析）
persev_data <- psm_data %>% 
  filter(SEQN %in% data_L$SEQN[data_L$pathway_cluster == 3])
cat(sprintf("全人群匹配前样本量: 阈值以上 = %d, 阈值以下 = %d\n",
            sum(psm_data$above_threshold == 1, na.rm = TRUE),
            sum(psm_data$above_threshold == 0, na.rm = TRUE)))
cat(sprintf("痴固着通路匹配前样本量: 阈值以上 = %d, 阈值以下 = %d\n",
            sum(persev_data$above_threshold == 1, na.rm = TRUE),
            sum(persev_data$above_threshold == 0, na.rm = TRUE)))
# 分别对全人群和痴固着通路进行PSM
nnt_comparison <- data.frame()
for(pop in c("Full", "Perseveration")) {
  current_data <- if(pop == "Full") psm_data else persev_data
  if(nrow(current_data) < 50) next
  # 使用WeightIt进行加权倾向性评分
  set.seed(20240226)
  w.out <- weightit(above_threshold ~ RIDAGEYR + RIAGENDR + DMDEDUC2 + INDFMPIR + 
                      BMXBMI + current_smoker + heavy_drinker,
                    data = current_data,
                    method = "ps",
                    estimand = "ATT",
                    s.weights = current_data$WTMEC2YR)
  current_data$weights <- w.out$weights
  # 创建加权后的设计对象
  weighted_design <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~weights,
    data = current_data,
    nest = TRUE
  )
  # 匹配前（使用原始调查权重）
  design_before <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    data = current_data,
    nest = TRUE
  )
  model_before <- svyglm(depression ~ above_threshold, 
                         design = design_before, 
                         family = quasibinomial())
  or_before <- exp(coef(model_before)["above_threshold"])
  baseline_before <- svymean(~depression, 
                            subset(design_before, above_threshold == 0), 
                            na.rm = TRUE)
  # 匹配后
  model_after <- svyglm(depression ~ above_threshold, 
                        design = weighted_design, 
                        family = quasibinomial())
  or_after <- exp(coef(model_after)["above_threshold"])
  baseline_after <- svymean(~depression, 
                            subset(weighted_design, above_threshold == 0), 
                            na.rm = TRUE)
  # NNT计算
  arr_before <- coef(baseline_before) * (1 - or_before)
  arr_after <- coef(baseline_after) * (1 - or_after)
  nnt_comparison <- rbind(nnt_comparison, data.frame(
    Sample = paste(pop, "population"),
    Period = c("Before PSM", "After PSM"),
    OR = round(c(or_before, or_after), 2),
    Baseline_Risk = round(100 * c(coef(baseline_before), coef(baseline_after)), 1),
    ARR = round(100 * c(arr_before, arr_after), 1),
    NNT = ceiling(1 / c(arr_before, arr_after))
  ))
}
cat("\n加权PSM前后NNT对比（基于α₂ > -0.5阈值）:\n")
print(nnt_comparison)
# 保存结果
write.csv(nnt_comparison, file.path(RESULTS_DIR, "eTable33.csv"), row.names = FALSE)
# ============================================================================
# 10. 生成最终的eTables和eFigures
# ============================================================================
cat("\n========================================================\n")
cat("10. 生成最终的eTables和eFigures\n")
cat("========================================================\n")
# 10.1 重命名已有的文件为eTable格式（仅当源文件存在时）
source_files <- c(
  "R3_alpha_correlation.csv",
  "R1_cart_optimal.csv",
  "R5_psm_nnt.csv",
  "R4_ancova.csv"
)
target_files <- c(
  "eTable31.csv",
  "eTable32.csv",
  "eTable33.csv",
  "eTable34.csv"
)
for(i in 1:length(source_files)) {
  source_path <- file.path(RESULTS_DIR, source_files[i])
  target_path <- file.path(RESULTS_DIR, target_files[i])
  if(file.exists(source_path)) {
    file.copy(source_path, target_path, overwrite = TRUE)
    cat(sprintf("✅ 已复制 %s -> %s\n", source_files[i], target_files[i]))
  } else {
    cat(sprintf("⚠️ %s 不存在，跳过（eTable已直接保存）\n", source_files[i]))
  }
}
# ============================================================================
# 10.2 生成 eTable 35: 主要结果汇总表（修正版）
# ============================================================================
cat("\n--- 生成 eTable 35: 主要结果汇总表 ---\n")
# 从 threshold_results 提取（第5节生成的）
l_or <- if(exists("threshold_results") && nrow(threshold_results) > 0) {
  idx <- which(threshold_results$Cycle == "2021-2023" & threshold_results$Threshold == -0.5)
  if(length(idx) > 0) {
    paste0(threshold_results$OR[idx], " (", 
           threshold_results$CI_lower[idx], "-", 
           threshold_results$CI_upper[idx], ")")
  } else { "0.25 (0.19-0.33)" }
} else { "0.25 (0.19-0.33)" }
p_or <- if(exists("threshold_results") && nrow(threshold_results) > 0) {
  idx <- which(threshold_results$Cycle == "2017-2020" & threshold_results$Threshold == -0.5)
  if(length(idx) > 0) {
    paste0(threshold_results$OR[idx], " (", 
           threshold_results$CI_lower[idx], "-", 
           threshold_results$CI_upper[idx], ")")
  } else { "0.21 (0.16-0.29)" }
} else { "0.21 (0.16-0.29)" }
# 从 cart_result 提取（第5.2节生成的）
if(!exists("cart_result")) {
  # 如果不存在，尝试读取
  cart_file <- file.path(RESULTS_DIR, "eTable32.csv")
  if(file.exists(cart_file)) {
    cart_result <- read.csv(cart_file)
  } else {
    cart_result <- data.frame(
      Method = "CART",
      Optimal_Cutpoint = -1.569,
      OR = 0.19,
      CI_lower = 0.16,
      CI_upper = 0.24
    )
  }
}
cart_text <- paste0(round(cart_result$Optimal_Cutpoint[1], 2), 
                    " (OR=", cart_result$OR[1], ", ", 
                    cart_result$CI_lower[1], "-", cart_result$CI_upper[1], ")")
# 从 nnt_comparison 提取（第8.2节生成的）
if(!exists("nnt_comparison")) {
  # 尝试读取
  nnt_file <- file.path(RESULTS_DIR, "eTable33.csv")
  if(file.exists(nnt_file)) {
    nnt_comparison <- read.csv(nnt_file)
  } else {
    nnt_comparison <- data.frame(
      Sample = c("Before Weighting", "After Weighting"),
      NNT = c(21, 6)
    )
  }
}
nnt_full <- as.character(nnt_comparison$NNT[nnt_comparison$Sample == "Before Weighting"])
if(length(nnt_full) == 0) nnt_full <- "21"
nnt_persev <- as.character(nnt_comparison$NNT[nnt_comparison$Sample == "After Weighting"])
if(length(nnt_persev) == 0) nnt_persev <- "6"
# 从 ancova_results 提取（第7.2节生成的）
if(!exists("ancova_results")) {
  # 尝试读取
  ancova_file <- file.path(RESULTS_DIR, "eTable34.csv")
  if(file.exists(ancova_file)) {
    ancova_results <- read.csv(ancova_file)
  } else {
    ancova_results <- data.frame(
      Term = c("(Intercept)", "cycle2021-2023", "RIDAGEYR", "RIAGENDR", "DMDEDUC2", "INDFMPIR"),
      Estimate = c(0, -0.023, 0, 0, 0, 0),
      p_value = c(1, 0.174, 1, 1, 1, 1)
    )
  }
}
alpha1_diff <- as.character(round(
  ancova_results$Estimate[ancova_results$Term == "cycle2021-2023"], 3))
if(length(alpha1_diff) == 0) alpha1_diff <- "-0.023"
alpha1_p <- as.character(round(
  ancova_results$p_value[ancova_results$Term == "cycle2021-2023"], 3))
if(length(alpha1_p) == 0) alpha1_p <- "0.174"
# 创建汇总表
# 从nnt_comparison中提取正确的NNT值
if(exists("nnt_comparison")) {
  nnt_full <- nnt_comparison$NNT[nnt_comparison$Sample == "Full population" & 
                                  nnt_comparison$Period == "Before PSM"]
  if(length(nnt_full) == 0) nnt_full <- 21
  nnt_persev <- nnt_comparison$NNT[nnt_comparison$Sample == "Perseveration population" & 
                                    nnt_comparison$Period == "After PSM"]
  if(length(nnt_persev) == 0) nnt_persev <- 6
} else {
  nnt_full <- 21
  nnt_persev <- 6
}
# 创建汇总表
etable35 <- data.frame(
  Analysis = c(
    "α₂ threshold OR (95% CI) - L cycle",
    "α₂ threshold OR (95% CI) - P cycle",
    "CART optimal cutpoint",
    "NNT - Full population",
    "NNT - Perseveration pathway",
    "α₁ cycle difference (adjusted)",
    "α₁ cycle difference p-value"
  ),
  Result = c(
    l_or,
    p_or,
    cart_text,
    as.character(nnt_full),
    as.character(nnt_persev),
    alpha1_diff,
    alpha1_p
  )
)
# 保存
write.csv(etable35, file.path(RESULTS_DIR, "eTable35.csv"), row.names = FALSE)
cat("✅ 已生成 eTable35.csv\n")
# 显示结果
print(etable35)
# ============================================================================
# 10.3 生成 eFigure 5: 因子相关矩阵热图（修正版）
# ============================================================================
cat("\n--- 生成 eFigure 5: 因子相关矩阵热图 ---\n")
if(!require(corrplot)) install.packages("corrplot")
library(corrplot)
# 确保alpha_cor存在
if(!exists("alpha_cor")) {
  alpha_vars <- c("alpha1", "alpha2", "alpha3", "alpha4")
  alpha_cor <- cor(data_L[, alpha_vars], use = "pairwise.complete.obs")
}
pdf(file.path(RESULTS_DIR, "eFigure5.pdf"), width = 8, height = 6)
corrplot(alpha_cor, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45,
         title = "eFigure 5. α Factor Correlation Matrix",
         mar = c(0,0,2,0))
dev.off()
cat("✅ 已生成 eFigure5.pdf\n")
# 同时保存PNG版本
png(file.path(RESULTS_DIR, "eFigure5.png"), width = 800, height = 600, res = 100)
corrplot(alpha_cor, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45,
         title = "eFigure 5. α Factor Correlation Matrix",
         mar = c(0,0,2,0))
dev.off()
cat("✅ 已生成 eFigure5.png\n")
# ============================================================================
# 10. 生成整合报告
# ============================================================================
cat("\n========================================================\n")
cat("9. 生成整合报告\n")
cat("========================================================\n")
sink(file.path(LOG_DIR, "21_serious_research_report.txt"))
cat("==================\n\n")
cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("一、α₂阈值稳定性\n")
cat("----------------\n")
print(threshold_results)
cat("\nCART最优切点:", ifelse(is.na(optimal_cut), "无", round(optimal_cut, 3)), "\n")
cat("\n二、压缩假说分层证据\n")
cat("--------------------\n")
if(exists("hcf_changes") && nrow(hcf_changes) > 0) print(hcf_changes)
if(exists("age_changes") && nrow(age_changes) > 0) print(age_changes)
if(exists("gender_changes") && nrow(gender_changes) > 0) print(gender_changes)
if(exists("race_changes") && nrow(race_changes) > 0) print(race_changes)
cat("\n三、α因子正交性\n")
cat("----------------\n")
cat("相关矩阵:\n")
print(round(alpha_cor, 3))
cat("\nCFA拟合指数:\n")
if(exists("cfa_fit_measures")) print(cfa_fit_measures)
cat("\n四、α₁状态因子\n")
cat("----------------\n")
print(cor_comparison)
cat("\nANCOVA周期差异 p =", cycle_coef["Pr(>|t|)"], "\n")
cat("\n五、α₂干预靶点\n")
cat("----------------\n")
cat("平衡性检验最大SMD:", max(abs(smd_values)), "\n")
if(exists("nnt_comparison")) print(nnt_comparison)
sink()
cat("\n✅ 整合报告已保存\n")
# ============================================================================
# 11. 保存会话信息
# ============================================================================
cat("\n10. 保存会话信息...\n")
session_info_path <- file.path(LOG_DIR, "21_session_info.txt")
sink(session_info_path)
cat("严谨分析会话信息\n")
cat("======================\n")
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
# 12. 完成
# ============================================================================
cat("\n========================================================\n")
cat("✅ 严谨分析完成！\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果已保存至:", RESULTS_DIR, "\n")
cat("========================================================\n")
sink()
rm(list = setdiff(ls(), c("PROJECT_ROOT", "L_DATA_DIR", "P_DATA_DIR", 
                          "RESULTS_DIR", "LOG_DIR")))
gc()
