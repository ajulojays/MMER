#!/usr/bin/env Rscript

# MMER Study 1 perturbability model
# --------------------------------
# Discovery target: mean Bray-Curtis displacement across T2/T3 for each
# independent cow-quarter trajectory (n = 20).
# Predictors: T1 richness, Shannon diversity, evenness, and dominance.
#
# Design principles:
#   * quarter trajectory, not post-baseline row, is the modeling unit
#   * all four quarters from a cow stay together in resampling folds
#   * ridge regularization handles correlated ecological predictors
#   * nested leave-one-cow-out CV estimates discovery performance
#   * cow-level bootstrap quantifies coefficient uncertainty
#   * final frozen model is intended for external validation, not clinical use

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

ROOT <- normalizePath(".", mustWork = TRUE)
INFILE <- file.path(ROOT, "results/production/questions/tables/Q3_baseline_predictor_dataset.csv")
OUTDIR <- file.path(ROOT, "results/production/perturbability_model")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

set.seed(20260812)
B <- as.integer(Sys.getenv("MMER_BOOTSTRAPS", "5000"))
if (!is.finite(B) || B < 100) stop("MMER_BOOTSTRAPS must be >=100")

FEATURES <- c(
  "baseline_richness",
  "baseline_shannon",
  "baseline_evenness",
  "baseline_dominance"
)
LAMBDA_GRID <- 10^seq(-4, 3, length.out = 80)

cat("\n============================================================\n")
cat(" MMER — STUDY 1 PERTURBABILITY MODEL\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Data assembly: collapse T2/T3 to one outcome per cow-quarter trajectory.
# -----------------------------------------------------------------------------
q3 <- read_csv(INFILE, show_col_types = FALSE, col_types = cols(cow_id = col_character()))
required <- c("cow_id", "quarter", "treatment", "timepoint", "bray_displacement", FEATURES)
missing <- setdiff(required, names(q3))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

trajectory <- q3 %>%
  filter(timepoint %in% c("T2", "T3")) %>%
  group_by(cow_id, quarter, treatment) %>%
  summarise(
    mean_bray_displacement = mean(bray_displacement, na.rm = TRUE),
    across(all_of(FEATURES), ~ first(.x)),
    n_post = n(),
    .groups = "drop"
  ) %>%
  mutate(cow_quarter = paste(cow_id, quarter, sep = "_"))

stopifnot(nrow(trajectory) == 20)
stopifnot(n_distinct(trajectory$cow_id) == 5)
stopifnot(all(trajectory$n_post == 2))
stopifnot(!anyNA(trajectory[, c("mean_bray_displacement", FEATURES)]))

write_csv(trajectory, file.path(OUTDIR, "study1_20trajectory_training_table.csv"))
cat("[PASS] 20 independent cow-quarter trajectories assembled\n")

X <- as.matrix(trajectory[, FEATURES])
y <- trajectory$mean_bray_displacement
cows <- trajectory$cow_id

# -----------------------------------------------------------------------------
# Minimal ridge implementation: intercept is unpenalized; predictors are
# standardized using training data only.
# -----------------------------------------------------------------------------
ridge_fit <- function(X, y, lambda) {
  mu <- colMeans(X)
  sig <- apply(X, 2, sd)
  if (any(!is.finite(sig) | sig <= 0)) stop("Zero/invalid predictor variance")
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
  ybar <- mean(y)
  p <- ncol(Z)
  beta_z <- solve(crossprod(Z) + lambda * diag(p), crossprod(Z, y - ybar))
  beta_z <- as.numeric(beta_z)
  names(beta_z) <- colnames(X)
  beta_raw <- beta_z / sig
  intercept_raw <- ybar - sum(beta_raw * mu)
  list(
    lambda = lambda,
    mean = mu,
    sd = sig,
    intercept_z = ybar,
    beta_z = beta_z,
    intercept_raw = intercept_raw,
    beta_raw = beta_raw
  )
}

ridge_predict <- function(fit, Xnew) {
  Z <- sweep(sweep(Xnew, 2, fit$mean, "-"), 2, fit$sd, "/")
  as.numeric(fit$intercept_z + Z %*% fit$beta_z)
}

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
mae <- function(obs, pred) mean(abs(obs - pred))

# -----------------------------------------------------------------------------
# Inner cow-blocked CV for lambda selection within an outer training set.
# -----------------------------------------------------------------------------
select_lambda_loco <- function(X, y, cows, grid) {
  ucows <- unique(cows)
  score <- lapply(grid, function(lambda) {
    pred <- rep(NA_real_, length(y))
    for (held in ucows) {
      tr <- cows != held
      te <- cows == held
      fit <- ridge_fit(X[tr, , drop = FALSE], y[tr], lambda)
      pred[te] <- ridge_predict(fit, X[te, , drop = FALSE])
    }
    tibble(lambda = lambda, RMSE = rmse(y, pred), MAE = mae(y, pred))
  }) %>% bind_rows()
  best <- score %>% arrange(RMSE, MAE, lambda) %>% slice(1)
  list(lambda = best$lambda[[1]], table = score)
}

# -----------------------------------------------------------------------------
# Nested leave-one-cow-out performance estimate.
# Outer cow is never used to select lambda.
# -----------------------------------------------------------------------------
cat("[STEP] Nested leave-one-cow-out validation\n")
outer_pred <- rep(NA_real_, length(y))
outer_lambda <- rep(NA_real_, length(y))

for (held in unique(cows)) {
  tr <- cows != held
  te <- cows == held
  inner <- select_lambda_loco(X[tr, , drop = FALSE], y[tr], cows[tr], LAMBDA_GRID)
  fit <- ridge_fit(X[tr, , drop = FALSE], y[tr], inner$lambda)
  outer_pred[te] <- ridge_predict(fit, X[te, , drop = FALSE])
  outer_lambda[te] <- inner$lambda
}

nested_predictions <- trajectory %>%
  transmute(
    cow_id, quarter, treatment, cow_quarter,
    observed_mean_bray = mean_bray_displacement,
    predicted_mean_bray = outer_pred,
    residual = observed_mean_bray - predicted_mean_bray,
    outer_fold_lambda = outer_lambda
  )
write_csv(nested_predictions, file.path(OUTDIR, "nested_LOCO_predictions.csv"))

nested_metrics <- tibble(
  n_trajectories = length(y),
  n_cows = n_distinct(cows),
  RMSE = rmse(y, outer_pred),
  MAE = mae(y, outer_pred),
  Pearson_r = suppressWarnings(cor(y, outer_pred, method = "pearson")),
  Spearman_rho = suppressWarnings(cor(y, outer_pred, method = "spearman"))
)
write_csv(nested_metrics, file.path(OUTDIR, "nested_LOCO_metrics.csv"))
cat("[PASS] nested LOCO predictions generated\n")

# -----------------------------------------------------------------------------
# Final frozen model. Lambda is selected by cow-blocked CV using the complete
# discovery cohort. This model is for external application after internal
# performance has already been estimated by the nested procedure above.
# -----------------------------------------------------------------------------
cat("[STEP] Fit frozen discovery model\n")
full_cv <- select_lambda_loco(X, y, cows, LAMBDA_GRID)
final_lambda <- full_cv$lambda
final_fit <- ridge_fit(X, y, final_lambda)
write_csv(full_cv$table, file.path(OUTDIR, "final_lambda_LOCO_grid.csv"))

coef_table <- tibble(
  term = c("(Intercept)", FEATURES),
  coefficient_raw = c(final_fit$intercept_raw, unname(final_fit$beta_raw)),
  coefficient_standardized = c(final_fit$intercept_z, unname(final_fit$beta_z))
)
write_csv(coef_table, file.path(OUTDIR, "frozen_model_coefficients.csv"))

scaling <- tibble(
  predictor = FEATURES,
  discovery_mean = unname(final_fit$mean),
  discovery_sd = unname(final_fit$sd)
)
write_csv(scaling, file.path(OUTDIR, "frozen_model_scaling.csv"))

model <- list(
  model_name = "MMER_Study1_Perturbability_Ridge_v1",
  target = "mean Bray-Curtis displacement across T2 and T3 from each quarter's T1 baseline",
  predictors = FEATURES,
  biological_unit = "cow-quarter trajectory",
  n_trajectories = nrow(trajectory),
  n_cows = n_distinct(cows),
  lambda = final_lambda,
  training_mean = final_fit$mean,
  training_sd = final_fit$sd,
  intercept_z = final_fit$intercept_z,
  beta_z = final_fit$beta_z,
  intercept_raw = final_fit$intercept_raw,
  beta_raw = final_fit$beta_raw,
  preprocessing = "predictors measured at T1; no outcome-derived feature engineering",
  intended_use = "external ecological perturbability validation only; not a mastitis-risk or clinical decision model"
)
saveRDS(model, file.path(OUTDIR, "MMER_Study1_Perturbability_Ridge_v1.rds"))

# -----------------------------------------------------------------------------
# Cow-level bootstrap uncertainty. All four quarter trajectories from a sampled
# cow are resampled together. Store raw-scale coefficient draws so an external
# cohort can propagate discovery-model uncertainty without outcome leakage.
# -----------------------------------------------------------------------------
cat("[STEP] Cow-level bootstrap (B = ", B, ")\n", sep = "")
ucows <- unique(cows)
boot <- vector("list", B)

for (b in seq_len(B)) {
  sampled <- sample(ucows, length(ucows), replace = TRUE)
  idx <- unlist(lapply(sampled, function(cc) which(cows == cc)), use.names = FALSE)
  fit_b <- tryCatch(ridge_fit(X[idx, , drop = FALSE], y[idx], final_lambda), error = function(e) NULL)
  if (is.null(fit_b)) next
  boot[[b]] <- tibble(
    bootstrap = b,
    intercept = fit_b$intercept_raw,
    baseline_richness = fit_b$beta_raw[["baseline_richness"]],
    baseline_shannon = fit_b$beta_raw[["baseline_shannon"]],
    baseline_evenness = fit_b$beta_raw[["baseline_evenness"]],
    baseline_dominance = fit_b$beta_raw[["baseline_dominance"]]
  )
}
boot_coef <- bind_rows(boot)
if (nrow(boot_coef) < 0.95 * B) warning("More than 5% bootstrap fits failed")
write_csv(boot_coef, file.path(OUTDIR, "cow_bootstrap_raw_coefficients.csv"))

boot_summary <- boot_coef %>%
  summarise(across(-bootstrap, list(
    median = ~ median(.x, na.rm = TRUE),
    lo95 = ~ quantile(.x, 0.025, na.rm = TRUE),
    hi95 = ~ quantile(.x, 0.975, na.rm = TRUE)
  )))
write_csv(boot_summary, file.path(OUTDIR, "cow_bootstrap_coefficient_summary.csv"))

# -----------------------------------------------------------------------------
# Compact manifest for external validation preregistration/model freezing.
# -----------------------------------------------------------------------------
manifest <- tibble(
  field = c(
    "model_name", "target", "n_cows", "n_trajectories", "final_lambda",
    "predictors", "internal_validation", "bootstrap_unit", "external_primary_use"
  ),
  value = c(
    model$model_name,
    model$target,
    as.character(model$n_cows),
    as.character(model$n_trajectories),
    format(model$lambda, scientific = TRUE),
    paste(FEATURES, collapse = ";"),
    "nested leave-one-cow-out CV",
    "cow (all four quarters retained together)",
    "rank/calibration of observed ecological displacement in an independent cohort"
  )
)
write_csv(manifest, file.path(OUTDIR, "model_manifest.csv"))

cat("\n============================================================\n")
cat(" PERTURBABILITY MODEL COMPLETE\n")
cat("============================================================\n\n")
print(nested_metrics)
cat("\nFinal lambda:", final_lambda, "\n")
cat("Bootstrap fits retained:", nrow(boot_coef), "/", B, "\n")
cat("\n[PASS] 20-trajectory target\n")
cat("[PASS] nested cow-blocked validation\n")
cat("[PASS] ridge regularization\n")
cat("[PASS] cow-level bootstrap uncertainty\n")
cat("[PASS] frozen external-validation model\n")
cat("[IMPORTANT] This predicts ecological perturbability, not mastitis risk.\n")
