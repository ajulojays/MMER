#!/usr/bin/env Rscript

# Apply frozen MMER Study 1 perturbability model to an external baseline cohort.
# Usage:
#   Rscript scripts/07b_apply_perturbability_model.R external_baseline.csv out_prefix
#
# Required columns:
#   sample_id, baseline_richness, baseline_shannon, baseline_evenness,
#   baseline_dominance
#
# The script does not use external outcomes. Predictions can therefore be frozen
# before observed follow-up displacement is calculated/unblinded.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/07b_apply_perturbability_model.R <external_baseline.csv> <out_prefix>")
}

ROOT <- normalizePath(".", mustWork = TRUE)
INPUT <- args[[1]]
PREFIX <- args[[2]]
MODEL_DIR <- file.path(ROOT, "results/production/perturbability_model")
MODEL_FILE <- file.path(MODEL_DIR, "MMER_Study1_Perturbability_Ridge_v1.rds")
BOOT_FILE <- file.path(MODEL_DIR, "cow_bootstrap_raw_coefficients.csv")

model <- readRDS(MODEL_FILE)
dat <- read_csv(INPUT, show_col_types = FALSE)
required <- c("sample_id", model$predictors)
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
if (anyNA(dat[, model$predictors])) stop("Predictor columns contain missing values")

X <- as.matrix(dat[, model$predictors])
Z <- sweep(sweep(X, 2, model$training_mean, "-"), 2, model$training_sd, "/")
pred <- as.numeric(model$intercept_z + Z %*% model$beta_z)

out <- dat %>%
  mutate(
    predicted_mean_bray_displacement = pred,
    predicted_perturbability_rank = rank(pred, ties.method = "average") / n()
  )

# Propagate discovery-model uncertainty with frozen raw-scale bootstrap draws.
if (file.exists(BOOT_FILE)) {
  boot <- read_csv(BOOT_FILE, show_col_types = FALSE)
  predmat <- matrix(NA_real_, nrow = nrow(dat), ncol = nrow(boot))
  for (b in seq_len(nrow(boot))) {
    predmat[, b] <- boot$intercept[b]
    for (f in model$predictors) predmat[, b] <- predmat[, b] + dat[[f]] * boot[[f]][b]
  }
  out$prediction_bootstrap_median <- apply(predmat, 1, median, na.rm = TRUE)
  out$prediction_bootstrap_lo95 <- apply(predmat, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$prediction_bootstrap_hi95 <- apply(predmat, 1, quantile, probs = 0.975, na.rm = TRUE)
}

outfile <- paste0(PREFIX, "_MMER_perturbability_predictions.csv")
write_csv(out, outfile)

cat("[PASS] Frozen MMER perturbability model applied\n")
cat("[INFO] External samples:", nrow(out), "\n")
cat("[INFO] No external outcome was used to generate predictions\n")
cat("[INFO] Output:", outfile, "\n")
