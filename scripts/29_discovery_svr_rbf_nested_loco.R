#!/usr/bin/env Rscript

# ============================================================================
# MMER — ITALY + MANITOBA DISCOVERY COHORT
# NESTED LEAVE-ONE-COW-OUT RBF SUPPORT VECTOR REGRESSION
# ============================================================================
#
# PURPOSE
# -------
# Test whether a nonlinear radial-basis-function support vector regression
# (RBF-SVR) can predict ecological resistance from baseline ecological
# architecture better than the current linear/ridge benchmark.
#
# BIOLOGICAL COHORT
# -----------------
# Frozen discovery cohort only:
#   Italy    : 20 quarter trajectories / 5 cows
#   Manitoba : 29 quarter trajectories / 9 cows
#   Total    : 49 quarter trajectories / 14 cows
#
# OUTCOME
# -------
# ecological_resistance = 1 - Bray-Curtis(T1, T2)
#
# FEATURES
# --------
# The script uses the already-frozen Italy+Manitoba baseline PCA coordinate
# system and compares three pre-defined architecture representations:
#   1. PC1-PC2   : frozen low-dimensional discovery representation
#   2. PC1-PC6   : best-performing ridge dimensionality sensitivity
#   3. PC1-PC31  : >=90% cumulative baseline PCA variance
#
# Study is included as a fixed binary input feature:
#   Italy = 0; Manitoba = 1.
#
# IMPORTANT INTERPRETATION
# ------------------------
# This is a PREDICTIVE SENSITIVITY ANALYSIS, not a replacement for the frozen
# CR2 inferential discovery model. Hyperparameters are tuned only within the
# training cows of each outer fold. The outer test cow is never used for
# hyperparameter selection.
#
# The PCA coordinate system itself is frozen from the discovery analysis and
# therefore is not re-estimated inside CV folds. Thus this analysis estimates
# prediction conditional on the frozen discovery PCA representation; it is not
# a fully nested feature-engineering pipeline.
#
# PRIMARY PERFORMANCE METRICS
# ---------------------------
#   RMSE
#   MAE
#   Pearson correlation(observed, predicted)
#   Spearman correlation(observed, predicted)
#   cross-validated R2
#
# MODEL
# -----
# epsilon-SVR with radial basis kernel:
#
#   K(x_i, x_j) = exp[-gamma * ||x_i - x_j||^2]
#
# Hyperparameters:
#   cost    : margin-violation penalty
#   gamma   : RBF width / locality
#   epsilon : width of epsilon-insensitive loss tube
#
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

if (!requireNamespace("e1071", quietly = TRUE)) {
  stop(
    "Package 'e1071' is required. Install in the mmer environment, e.g.\n",
    "  conda install -c conda-forge r-e1071\n",
    "then rerun this script."
  )
}

set.seed(290821)

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"

PCA_FILE <- file.path(
  ROOT,
  "results/primary_healthy_twoAims/models/baseline_withinStudy_CLR_PCA.rds"
)

ANALYSIS_FILE <- file.path(
  ROOT,
  "results/primary_healthy_twoAims/tables/PRIMARY_analysis_table.csv"
)

OUT <- file.path(
  ROOT,
  "results/primary_healthy_twoAims/SVR_RBF_nested_LOCO"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================================
# 1. LOAD FROZEN DISCOVERY DATA
# ============================================================================

pca <- readRDS(PCA_FILE)

dat0 <- read_csv(
  ANALYSIS_FILE,
  show_col_types = FALSE
)

stopifnot(
  nrow(dat0) == 49,
  n_distinct(dat0$cow_global) == 14,
  all(c(
    "trajectory_global",
    "cow_global",
    "study",
    "ecological_resistance"
  ) %in% names(dat0))
)

# ============================================================================
# 2. PCA VARIANCE COVERAGE
# ============================================================================

ve <- pca$sdev^2 / sum(pca$sdev^2)
cum <- cumsum(ve)

K90 <- which(cum >= 0.90)[1]

cat("\n====================================\n")
cat("DISCOVERY PCA COVERAGE\n")
cat("====================================\n")
cat("PC1-PC2 :", round(100 * cum[2], 2), "%\n")
cat("PC1-PC6 :", round(100 * cum[6], 2), "%\n")
cat("PC1-PC31:", round(100 * cum[31], 2), "%\n")
cat("First K reaching >=90%:", K90, "\n")

stopifnot(K90 == 31)

# ============================================================================
# 3. EXTRACT AND STANDARDIZE FROZEN PCA SCORES
# ============================================================================

scores <- as.data.frame(
  pca$x[, 1:31, drop = FALSE]
) %>%
  rownames_to_column("trajectory_global")

for (k in 1:31) {
  pc <- paste0("PC", k)
  scores[[paste0("z_PC", k)]] <- as.numeric(scale(scores[[pc]]))
}

dat <- dat0 %>%
  select(-any_of(c("PC1", "PC2", "z_PC1", "z_PC2"))) %>%
  left_join(scores, by = "trajectory_global") %>%
  mutate(
    cow_global = as.character(cow_global),
    study = as.character(study)
  )

needed <- c(
  "ecological_resistance",
  "cow_global",
  "study",
  paste0("z_PC", 1:31)
)

stopifnot(
  nrow(dat) == 49,
  n_distinct(dat$cow_global) == 14,
  !anyNA(dat[, needed])
)

# ============================================================================
# 4. DESIGN MATRIX
# ============================================================================

make_matrix <- function(data, k) {

  pcs <- paste0("z_PC", seq_len(k))

  Xpc <- as.matrix(
    data[, pcs, drop = FALSE]
  )

  storage.mode(Xpc) <- "double"

  study_manitoba <- as.numeric(
    data$study == "Manitoba"
  )

  X <- cbind(
    Xpc,
    studyManitoba = study_manitoba
  )

  storage.mode(X) <- "double"
  X
}

# ============================================================================
# 5. PERFORMANCE METRICS
# ============================================================================

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

mae <- function(obs, pred) {
  mean(abs(obs - pred))
}

cv_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

safe_cor <- function(x, y, method = "pearson") {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = method)
}

# ============================================================================
# 6. HYPERPARAMETER GRID
#
# Cost spans weak to strong penalty.
# Gamma is scaled relative to input dimensionality so that PC2, PC6 and PC31
# models are searched on comparable kernel-width scales.
# ============================================================================

cost_grid <- c(0.25, 1, 4, 16, 64, 256)
gamma_multiplier_grid <- c(0.25, 0.5, 1, 2, 4)
epsilon_grid <- c(0.025, 0.05, 0.10)

make_grid <- function(k) {

  p <- k + 1  # PCs + study indicator

  expand.grid(
    cost = cost_grid,
    gamma_multiplier = gamma_multiplier_grid,
    epsilon = epsilon_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      gamma = gamma_multiplier / p
    )
}

# ============================================================================
# 7. INNER LEAVE-ONE-COW-OUT TUNING
# ============================================================================

choose_hyperparameters_inner <- function(train_data, k) {

  inner_cows <- unique(train_data$cow_global)
  grid <- make_grid(k)

  score_one <- function(cost_value, gamma_value, epsilon_value) {

    pred <- rep(NA_real_, nrow(train_data))

    for (held_cow in inner_cows) {

      idx_val <- which(train_data$cow_global == held_cow)
      idx_train <- which(train_data$cow_global != held_cow)

      tr <- train_data[idx_train, , drop = FALSE]
      va <- train_data[idx_val, , drop = FALSE]

      Xtr <- make_matrix(tr, k)
      Xva <- make_matrix(va, k)

      fit <- e1071::svm(
        x = Xtr,
        y = tr$ecological_resistance,
        type = "eps-regression",
        kernel = "radial",
        cost = cost_value,
        gamma = gamma_value,
        epsilon = epsilon_value,
        scale = FALSE
      )

      pred[idx_val] <- as.numeric(
        predict(fit, Xva)
      )
    }

    tibble(
      cost = cost_value,
      gamma = gamma_value,
      epsilon = epsilon_value,
      inner_RMSE = rmse(train_data$ecological_resistance, pred),
      inner_MAE = mae(train_data$ecological_resistance, pred),
      inner_CV_R2 = cv_r2(train_data$ecological_resistance, pred)
    )
  }

  results <- bind_rows(
    lapply(seq_len(nrow(grid)), function(i) {
      score_one(
        grid$cost[i],
        grid$gamma[i],
        grid$epsilon[i]
      )
    })
  )

  # Primary tuning criterion = minimum inner cow-held-out RMSE.
  # Ties are resolved toward the simpler/more regularized solution:
  # lower cost, lower gamma, larger epsilon.
  best <- results %>%
    arrange(
      inner_RMSE,
      cost,
      gamma,
      desc(epsilon)
    ) %>%
    slice(1)

  list(best = best, table = results)
}

# ============================================================================
# 8. OUTER LEAVE-ONE-COW-OUT EVALUATION
# ============================================================================

run_nested_loco <- function(k) {

  cows <- unique(dat$cow_global)
  predictions <- list()
  tuning_records <- list()

  cat("\n====================================\n")
  cat("RBF-SVR PC1-PC", k, "\n", sep = "")
  cat("====================================\n")

  for (held_cow in cows) {

    cat("Outer held-out cow:", held_cow, "\n")

    train <- dat %>%
      filter(cow_global != held_cow)

    test <- dat %>%
      filter(cow_global == held_cow)

    tuned <- choose_hyperparameters_inner(train, k)
    best <- tuned$best

    Xtrain <- make_matrix(train, k)
    Xtest <- make_matrix(test, k)

    fit <- e1071::svm(
      x = Xtrain,
      y = train$ecological_resistance,
      type = "eps-regression",
      kernel = "radial",
      cost = best$cost,
      gamma = best$gamma,
      epsilon = best$epsilon,
      scale = FALSE
    )

    pred <- as.numeric(
      predict(fit, Xtest)
    )

    predictions[[held_cow]] <- tibble(
      model = paste0("PC1-PC", k, " RBF-SVR"),
      n_PCs = k,
      trajectory_global = test$trajectory_global,
      cow_global = test$cow_global,
      study = test$study,
      observed = test$ecological_resistance,
      predicted = pred
    )

    tuning_records[[held_cow]] <- tibble(
      model = paste0("PC1-PC", k, " RBF-SVR"),
      n_PCs = k,
      held_out_cow = held_cow,
      cost = best$cost,
      gamma = best$gamma,
      epsilon = best$epsilon,
      inner_RMSE = best$inner_RMSE,
      inner_MAE = best$inner_MAE,
      inner_CV_R2 = best$inner_CV_R2
    )
  }

  pred <- bind_rows(predictions)
  tune <- bind_rows(tuning_records)

  performance <- tibble(
    model = paste0("PC1-PC", k, " RBF-SVR"),
    n_PCs = k,
    cumulative_PCA_variance_percent = 100 * cum[k],
    RMSE = rmse(pred$observed, pred$predicted),
    MAE = mae(pred$observed, pred$predicted),
    Pearson_r = safe_cor(pred$observed, pred$predicted, "pearson"),
    Spearman_rho = safe_cor(pred$observed, pred$predicted, "spearman"),
    CV_R2 = cv_r2(pred$observed, pred$predicted)
  )

  list(
    performance = performance,
    predictions = pred,
    tuning = tune
  )
}

# ============================================================================
# 9. RUN PRE-DEFINED REPRESENTATIONS
# ============================================================================

res2 <- run_nested_loco(2)
res6 <- run_nested_loco(6)
res31 <- run_nested_loco(31)

performance <- bind_rows(
  res2$performance,
  res6$performance,
  res31$performance
)

predictions <- bind_rows(
  res2$predictions,
  res6$predictions,
  res31$predictions
)

tuning <- bind_rows(
  res2$tuning,
  res6$tuning,
  res31$tuning
)

cat("\n====================================\n")
cat("RBF-SVR NESTED LOCO PERFORMANCE\n")
cat("====================================\n")
print(performance, width = Inf)

cat("\nSelected hyperparameters by outer fold:\n")
print(tuning, width = Inf)

# ============================================================================
# 10. FIT FINAL FULL-DATA MODELS FOR REPRODUCIBILITY ONLY
#
# These fits are NOT used for reported generalization performance.
# ============================================================================

fit_final <- function(k) {

  tuned <- choose_hyperparameters_inner(dat, k)
  best <- tuned$best
  X <- make_matrix(dat, k)

  fit <- e1071::svm(
    x = X,
    y = dat$ecological_resistance,
    type = "eps-regression",
    kernel = "radial",
    cost = best$cost,
    gamma = best$gamma,
    epsilon = best$epsilon,
    scale = FALSE
  )

  list(
    fit = fit,
    best = best,
    tuning = tuned$table
  )
}

final2 <- fit_final(2)
final6 <- fit_final(6)
final31 <- fit_final(31)

final_parameters <- bind_rows(
  final2$best %>% mutate(model = "PC1-PC2 RBF-SVR", .before = 1),
  final6$best %>% mutate(model = "PC1-PC6 RBF-SVR", .before = 1),
  final31$best %>% mutate(model = "PC1-PC31 RBF-SVR", .before = 1)
)

cat("\n====================================\n")
cat("FINAL FULL-DATA TUNED PARAMETERS\n")
cat("====================================\n")
print(final_parameters, width = Inf)

# ============================================================================
# 11. SAVE OUTPUTS
# ============================================================================

write_csv(
  performance,
  file.path(OUT, "SVR_RBF_nested_LOCO_performance.csv")
)

write_csv(
  predictions,
  file.path(OUT, "SVR_RBF_nested_LOCO_predictions.csv")
)

write_csv(
  tuning,
  file.path(OUT, "SVR_RBF_outer_fold_selected_hyperparameters.csv")
)

write_csv(
  final_parameters,
  file.path(OUT, "SVR_RBF_final_full_data_hyperparameters.csv")
)

write_csv(
  final2$tuning,
  file.path(OUT, "SVR_RBF_PC2_full_data_inner_grid.csv")
)

write_csv(
  final6$tuning,
  file.path(OUT, "SVR_RBF_PC6_full_data_inner_grid.csv")
)

write_csv(
  final31$tuning,
  file.path(OUT, "SVR_RBF_PC31_full_data_inner_grid.csv")
)

saveRDS(
  final2$fit,
  file.path(OUT, "SVR_RBF_PC2_final_model.rds")
)

saveRDS(
  final6$fit,
  file.path(OUT, "SVR_RBF_PC6_final_model.rds")
)

saveRDS(
  final31$fit,
  file.path(OUT, "SVR_RBF_PC31_final_model.rds")
)

variance_table <- tibble(
  PC = paste0("PC", seq_along(ve)),
  variance_percent = 100 * ve,
  cumulative_percent = 100 * cum
)

write_csv(
  variance_table,
  file.path(OUT, "PCA_variance_full.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "sessionInfo.txt")
)

cat("\n====================================\n")
cat("[PASS] DISCOVERY RBF-SVR NESTED LOCO COMPLETE\n")
cat("====================================\n")
cat("Trajectories: 49\n")
cat("Independent cows: 14\n")
cat("Studies: Italy + Manitoba\n")
cat("Models: PC1-PC2, PC1-PC6, PC1-PC31\n")
cat("Outer validation: leave-one-cow-out\n")
cat("Inner tuning: leave-one-cow-out within training cows\n")
cat("Kernel: radial basis function\n")
cat("Output:", OUT, "\n")
