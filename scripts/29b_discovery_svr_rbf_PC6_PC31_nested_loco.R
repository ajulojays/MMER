#!/usr/bin/env Rscript

# =============================================================================
# MMER — ITALY + MANITOBA DISCOVERY COHORT
# RBF-SVR BENCHMARK: PC1-PC6 vs PC1-PC31
# =============================================================================
#
# PURPOSE
# -------
# Compare two fixed baseline-architecture representations for prediction of
# ecological resistance in the frozen Italy + Manitoba discovery cohort:
#
#   Model A: PC1-PC6  (42.71% cumulative baseline PCA variance)
#   Model B: PC1-PC31 (90.93% cumulative baseline PCA variance)
#
# MODEL
# -----
# epsilon-support vector regression with radial basis kernel:
#
#   K(x_i, x_j) = exp[-gamma * ||x_i - x_j||^2]
#
# VALIDATION
# ----------
# Nested leave-one-cow-out cross-validation (LOCO):
#   * outer loop: one cow held out completely for final prediction;
#   * inner loop: LOCO among training cows only for cost/gamma/epsilon tuning;
#   * all PC scaling is learned from the relevant training split only and then
#     applied to validation/test observations;
#   * the held-out cow is never used for hyperparameter tuning or scaling.
#
# IMPORTANT
# ---------
# The PCA coordinate system itself is the frozen discovery PCA fitted before
# this predictive sensitivity analysis. Therefore this is prediction
# conditional on the frozen PCA representation, not a fully nested PCA pipeline.
# The original PC1+PC2 CR2 model remains the primary inferential discovery model.
#
# COHORT
# ------
# Italy    : 20 quarter trajectories / 5 cows
# Manitoba : 29 quarter trajectories / 9 cows
# Total    : 49 trajectories / 14 cows
#
# OUTCOME
# -------
# ecological_resistance = 1 - Bray-Curtis(T1, T2)
#
# OUTPUTS
# -------
# results/primary_healthy_twoAims/SVR_RBF_PC6_PC31_nested_LOCO/
#   SVR_PC6_PC31_nested_LOCO_performance.csv
#   SVR_PC6_PC31_nested_LOCO_predictions.csv
#   SVR_PC6_PC31_outer_fold_tuning.csv
#   SVR_PC6_PC31_final_hyperparameters.csv
#   SVR_PC6_PC31_full_grid_final_training.csv
#   sessionInfo.txt
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

if (!requireNamespace("e1071", quietly = TRUE)) {
  stop(
    "Package 'e1071' is required. Install it with:\n",
    "  conda install -c conda-forge r-e1071\n"
  )
}

set.seed(290822)

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
  "results/primary_healthy_twoAims/SVR_RBF_PC6_PC31_nested_LOCO"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(PCA_FILE), file.exists(ANALYSIS_FILE))

# =============================================================================
# 1. LOAD FROZEN DISCOVERY DATA
# =============================================================================

pca <- readRDS(PCA_FILE)

dat0 <- read_csv(
  ANALYSIS_FILE,
  show_col_types = FALSE
)

stopifnot(
  nrow(dat0) == 49,
  dplyr::n_distinct(dat0$cow_global) == 14,
  all(c(
    "trajectory_global",
    "cow_global",
    "study",
    "ecological_resistance"
  ) %in% names(dat0))
)

# =============================================================================
# 2. PCA COVERAGE + RAW FROZEN SCORES
# =============================================================================

ve <- pca$sdev^2 / sum(pca$sdev^2)
cum <- cumsum(ve)
K90 <- which(cum >= 0.90)[1]

cat("\n====================================\n")
cat("DISCOVERY PCA COVERAGE\n")
cat("====================================\n")
cat("PC1-PC6 :", round(100 * cum[6], 2), "%\n")
cat("PC1-PC31:", round(100 * cum[31], 2), "%\n")
cat("First K reaching >=90%:", K90, "\n")

stopifnot(K90 == 31)

scores <- as.data.frame(
  pca$x[, 1:31, drop = FALSE]
) %>%
  rownames_to_column("trajectory_global")

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
  paste0("PC", 1:31)
)

stopifnot(
  nrow(dat) == 49,
  dplyr::n_distinct(dat$cow_global) == 14,
  !anyNA(dat[, needed])
)

# =============================================================================
# 3. HELPERS
# =============================================================================

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
mae  <- function(obs, pred) mean(abs(obs - pred))

cv_r2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

safe_cor <- function(x, y, method = "pearson") {
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = method)
}

# Learn PC scaling ONLY from training observations.
fit_scaler <- function(train_data, k) {
  pcs <- paste0("PC", seq_len(k))
  mu <- vapply(train_data[, pcs, drop = FALSE], mean, numeric(1))
  ss <- vapply(train_data[, pcs, drop = FALSE], sd, numeric(1))
  if (any(!is.finite(ss)) || any(ss <= 0)) {
    stop("Non-finite or zero PC SD in training fold.")
  }
  list(mean = mu, sd = ss)
}

make_matrix <- function(data, k, scaler) {
  pcs <- paste0("PC", seq_len(k))
  X <- as.matrix(data[, pcs, drop = FALSE])
  storage.mode(X) <- "double"

  X <- sweep(X, 2, scaler$mean, "-")
  X <- sweep(X, 2, scaler$sd, "/")

  # Fixed study encoding avoids one-level factor failures in held-out subsets.
  study_manitoba <- as.numeric(data$study == "Manitoba")
  X <- cbind(X, studyManitoba = study_manitoba)
  storage.mode(X) <- "double"
  X
}

# =============================================================================
# 4. HYPERPARAMETER GRID
# =============================================================================

# Deliberately compact but broad grid suitable for n=14 independent cows.
cost_grid <- c(0.25, 1, 4, 16, 64, 256)
gamma_multiplier_grid <- c(0.25, 0.5, 1, 2, 4)
epsilon_grid <- c(0.025, 0.05, 0.10)

make_grid <- function(k) {
  p <- k + 1L  # PCs + study indicator
  expand.grid(
    cost = cost_grid,
    gamma_multiplier = gamma_multiplier_grid,
    epsilon = epsilon_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(gamma = gamma_multiplier / p)
}

# =============================================================================
# 5. INNER LOCO HYPERPARAMETER TUNING
# =============================================================================

choose_hyperparameters_inner <- function(train_data, k) {

  inner_cows <- unique(train_data$cow_global)
  grid <- make_grid(k)

  score_one <- function(cost_value, gamma_value, epsilon_value) {

    pred <- rep(NA_real_, nrow(train_data))

    for (held_cow in inner_cows) {

      idx_val <- which(train_data$cow_global == held_cow)
      idx_fit <- which(train_data$cow_global != held_cow)

      tr <- train_data[idx_fit, , drop = FALSE]
      va <- train_data[idx_val, , drop = FALSE]

      scaler <- fit_scaler(tr, k)
      Xtr <- make_matrix(tr, k, scaler)
      Xva <- make_matrix(va, k, scaler)

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

      pred[idx_val] <- as.numeric(predict(fit, Xva))
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
  # Tie-breakers favor stronger regularization/smoother solutions.
  best <- results %>%
    arrange(inner_RMSE, cost, gamma, desc(epsilon)) %>%
    slice(1)

  list(best = best, table = results)
}

# =============================================================================
# 6. OUTER NESTED LOCO
# =============================================================================

run_nested_loco <- function(k) {

  cows <- unique(dat$cow_global)
  predictions <- list()
  tuning <- list()

  cat("\n====================================\n")
  cat("RBF-SVR PC1-PC", k, "\n", sep = "")
  cat("====================================\n")

  for (held_cow in cows) {

    cat("Outer held-out cow:", held_cow, "\n")

    train <- dat %>% filter(cow_global != held_cow)
    test  <- dat %>% filter(cow_global == held_cow)

    tuned <- choose_hyperparameters_inner(train, k)
    best <- tuned$best

    # Outer-fold scaling uses OUTER TRAINING DATA ONLY.
    scaler <- fit_scaler(train, k)
    Xtrain <- make_matrix(train, k, scaler)
    Xtest  <- make_matrix(test,  k, scaler)

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

    yhat <- as.numeric(predict(fit, Xtest))

    predictions[[held_cow]] <- tibble(
      model = paste0("PC1-PC", k, " RBF-SVR"),
      n_PCs = k,
      trajectory_global = test$trajectory_global,
      cow_global = test$cow_global,
      study = test$study,
      observed = test$ecological_resistance,
      predicted = yhat
    )

    tuning[[held_cow]] <- tibble(
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
  tune <- bind_rows(tuning)

  stopifnot(
    nrow(pred) == 49,
    !anyNA(pred$predicted)
  )

  perf <- tibble(
    model = paste0("PC1-PC", k, " RBF-SVR"),
    n_PCs = k,
    cumulative_PCA_variance_percent = 100 * cum[k],
    RMSE = rmse(pred$observed, pred$predicted),
    MAE = mae(pred$observed, pred$predicted),
    Pearson_r = safe_cor(pred$observed, pred$predicted, "pearson"),
    Spearman_rho = safe_cor(pred$observed, pred$predicted, "spearman"),
    CV_R2 = cv_r2(pred$observed, pred$predicted)
  )

  list(performance = perf, predictions = pred, tuning = tune)
}

# =============================================================================
# 7. RUN THE TWO PRESPECIFIED REPRESENTATIONS
# =============================================================================

res6  <- run_nested_loco(6)
res31 <- run_nested_loco(31)

performance <- bind_rows(
  res6$performance,
  res31$performance
) %>%
  arrange(RMSE)

predictions <- bind_rows(
  res6$predictions,
  res31$predictions
)

outer_tuning <- bind_rows(
  res6$tuning,
  res31$tuning
)

cat("\n====================================\n")
cat("RBF-SVR PC6 vs PC31 NESTED LOCO PERFORMANCE\n")
cat("====================================\n")
print(performance, width = Inf)

cat("\nCurrent reference benchmark:\n")
cat("PC1-PC6 ridge nested LOCO: CV_R2 = 0.364; RMSE = 0.118\n")

# =============================================================================
# 8. FINAL FULL-DISCOVERY MODELS FOR HYPERPARAMETER RECORD ONLY
#
# These are NOT used to estimate generalization performance.
# =============================================================================

final_records <- list()
final_grid_records <- list()

for (k in c(6, 31)) {

  tuned <- choose_hyperparameters_inner(dat, k)
  best <- tuned$best

  final_records[[as.character(k)]] <- tibble(
    model = paste0("PC1-PC", k, " RBF-SVR"),
    n_PCs = k,
    cumulative_PCA_variance_percent = 100 * cum[k],
    cost = best$cost,
    gamma = best$gamma,
    epsilon = best$epsilon,
    full_data_inner_LOCO_RMSE = best$inner_RMSE,
    full_data_inner_LOCO_MAE = best$inner_MAE,
    full_data_inner_LOCO_CV_R2 = best$inner_CV_R2
  )

  final_grid_records[[as.character(k)]] <- tuned$table %>%
    mutate(
      model = paste0("PC1-PC", k, " RBF-SVR"),
      n_PCs = k,
      .before = 1
    )
}

final_hyperparameters <- bind_rows(final_records)
final_grid <- bind_rows(final_grid_records)

# =============================================================================
# 9. SAVE
# =============================================================================

write_csv(
  performance,
  file.path(OUT, "SVR_PC6_PC31_nested_LOCO_performance.csv")
)

write_csv(
  predictions,
  file.path(OUT, "SVR_PC6_PC31_nested_LOCO_predictions.csv")
)

write_csv(
  outer_tuning,
  file.path(OUT, "SVR_PC6_PC31_outer_fold_tuning.csv")
)

write_csv(
  final_hyperparameters,
  file.path(OUT, "SVR_PC6_PC31_final_hyperparameters.csv")
)

write_csv(
  final_grid,
  file.path(OUT, "SVR_PC6_PC31_full_grid_final_training.csv")
)

capture.output(
  sessionInfo(),
  file = file.path(OUT, "sessionInfo.txt")
)

cat("\n====================================\n")
cat("[PASS] PC1-PC6 / PC1-PC31 RBF-SVR COMPLETE\n")
cat("====================================\n")
cat("Discovery trajectories: 49\n")
cat("Independent cows: 14\n")
cat("Outer validation: leave-one-cow-out\n")
cat("Inner tuning: leave-one-cow-out among training cows\n")
cat("Fold-local scaling: YES\n")
cat("PC1-PC6 variance:", round(100 * cum[6], 2), "%\n")
cat("PC1-PC31 variance:", round(100 * cum[31], 2), "%\n")
cat("Output:", OUT, "\n")
