#!/usr/bin/env Rscript

# MMER — Van Beeck external validation only
# Purpose: independent phenomenon-level validation of the Italy + Manitoba
# baseline-architecture -> ecological-resistance hypothesis.
# No subgroup, SCC, treatment-interaction, or temporal-endpoint hypothesis tests.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(vegan)
  library(sandwich)
})

set.seed(260817)

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"
META <- file.path(ROOT, "results/vanbeeck/analysis_3k/metadata_60_complete_trajectories_locked.csv")
SEQ  <- file.path(ROOT, "results/vanbeeck/analysis_3k/seqtab_bacterial_rarefied_3000.rds")
TAX  <- file.path(ROOT, "results/vanbeeck/taxonomy/taxonomy_bacterial.rds")
OUT  <- file.path(ROOT, "results/vanbeeck/external_validation_only")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(META), file.exists(SEQ), file.exists(TAX))

# -----------------------------------------------------------------------------
# 1. Locked cohort
# -----------------------------------------------------------------------------
meta <- read_csv(META, show_col_types = FALSE) %>%
  mutate(
    sample_name = as.character(sample_name),
    core = as.character(core),
    timepoint = as.character(timepoint),
    dairy = factor(dairy, levels = c("Dairy 1", "Dairy 2", "Dairy 3")),
    treatment = factor(
      treatment,
      levels = c("control low SCC", "control high SCC", "CB", "CH")
    )
  )

design <- meta %>% distinct(core, dairy, treatment)
required_timepoints <- c("Baseline", "7 Days", "55-75DIM")

stopifnot(
  nrow(meta) == 180,
  n_distinct(meta$core) == 60,
  nrow(design) == 60,
  all(table(meta$core) == 3),
  !anyNA(meta$dairy),
  !anyNA(meta$treatment)
)

for (id in unique(meta$core)) {
  z <- meta %>% filter(core == id)
  stopifnot(setequal(z$timepoint, required_timepoints))
}

expected_counts <- c(
  "control low SCC" = 14,
  "control high SCC" = 12,
  "CB" = 18,
  "CH" = 16
)
observed_counts <- design %>% count(treatment)
observed_named <- setNames(observed_counts$n, as.character(observed_counts$treatment))
stopifnot(all(observed_named[names(expected_counts)] == expected_counts))

cat("\n====================================\n")
cat("VAN BEECK LOCKED EXTERNAL COHORT\n")
cat("====================================\n")
cat("Samples:", nrow(meta), "\n")
cat("Trajectories:", n_distinct(meta$core), "\n")
print(design %>% count(dairy, treatment), width = Inf)
write_csv(design, file.path(OUT, "locked_60_trajectory_design.csv"))

# -----------------------------------------------------------------------------
# 2. ASV table and taxonomy
# -----------------------------------------------------------------------------
seqtab_all <- readRDS(SEQ)
tax_all <- readRDS(TAX)
stopifnot(
  all(meta$sample_name %in% rownames(seqtab_all)),
  all(colnames(seqtab_all) %in% rownames(tax_all))
)
seqtab <- seqtab_all[meta$sample_name, , drop = FALSE]
tax <- tax_all[colnames(seqtab), , drop = FALSE]
stopifnot(nrow(seqtab) == 180, all(rowSums(seqtab) == 3000))

# Remove organellar features before family aggregation.
tax_text <- apply(tax, 1, paste, collapse = " ")
organellar <- grepl("mitochond|chloroplast|plastid|organelle", tax_text, ignore.case = TRUE)
seqtab <- seqtab[, !organellar, drop = FALSE]
tax <- tax[!organellar, , drop = FALSE]
stopifnot(ncol(seqtab) > 0, all(rowSums(seqtab) > 0))

# -----------------------------------------------------------------------------
# 3. Aggregate to bacterial family
# -----------------------------------------------------------------------------
required_tax <- c("Phylum", "Class", "Order", "Family")
stopifnot(all(required_tax %in% colnames(tax)))
family_label <- as.character(tax[, "Family"])

for (i in seq_along(family_label)) {
  if (is.na(family_label[i]) || trimws(family_label[i]) == "") {
    ord <- as.character(tax[i, "Order"])
    cls <- as.character(tax[i, "Class"])
    phy <- as.character(tax[i, "Phylum"])
    if (!is.na(ord) && trimws(ord) != "") {
      family_label[i] <- paste0("Unclassified_Order__", ord)
    } else if (!is.na(cls) && trimws(cls) != "") {
      family_label[i] <- paste0("Unclassified_Class__", cls)
    } else if (!is.na(phy) && trimws(phy) != "") {
      family_label[i] <- paste0("Unclassified_Phylum__", phy)
    } else {
      family_label[i] <- "Unclassified_Bacteria"
    }
  }
}

families <- sort(unique(family_label))
fam_counts <- matrix(
  0,
  nrow = nrow(seqtab),
  ncol = length(families),
  dimnames = list(rownames(seqtab), families)
)
for (f in families) {
  fam_counts[, f] <- rowSums(seqtab[, family_label == f, drop = FALSE])
}
stopifnot(all(rowSums(fam_counts) > 0))
fam_ra <- fam_counts / rowSums(fam_counts)
stopifnot(max(abs(rowSums(fam_ra) - 1)) < 1e-10)

cat("ASVs after organellar removal:", ncol(seqtab), "\n")
cat("Family matrix:", nrow(fam_ra), "samples x", ncol(fam_ra), "families\n")

# -----------------------------------------------------------------------------
# 4. Primary ecological-resistance phenotype
# -----------------------------------------------------------------------------
# early_R = 1 - Bray(Baseline, 7 Days)
# late_R  = 1 - Bray(Baseline, 55-75DIM)
# overall_R = mean(early_R, late_R)
# The components construct the frozen primary outcome; they are not separately
# hypothesis-tested in this external-validation-only script.
trajectory_resistance <- lapply(design$core, function(id) {
  m <- meta %>% filter(core == id)
  b <- m$sample_name[m$timepoint == "Baseline"]
  e <- m$sample_name[m$timepoint == "7 Days"]
  l <- m$sample_name[m$timepoint == "55-75DIM"]
  stopifnot(length(b) == 1, length(e) == 1, length(l) == 1)
  dm <- as.matrix(vegdist(fam_ra[c(b, e, l), , drop = FALSE], method = "bray"))
  d_early <- unname(dm[b, e])
  d_late <- unname(dm[b, l])
  early_R <- 1 - d_early
  late_R <- 1 - d_late
  tibble(
    core = id,
    baseline_sample = b,
    day7_sample = e,
    late_sample = l,
    displacement_early = d_early,
    displacement_late = d_late,
    early_resistance = early_R,
    longterm_resistance = late_R,
    overall_resistance = mean(c(early_R, late_R))
  )
}) %>% bind_rows() %>% left_join(design, by = "core")

stopifnot(
  nrow(trajectory_resistance) == 60,
  !anyNA(trajectory_resistance$overall_resistance),
  all(trajectory_resistance$overall_resistance >= 0),
  all(trajectory_resistance$overall_resistance <= 1)
)

write_csv(
  trajectory_resistance,
  file.path(OUT, "vanbeeck_family_ecological_resistance_60.csv")
)

resistance_summary <- trajectory_resistance %>%
  summarise(
    n = n(),
    mean = mean(overall_resistance),
    sd = sd(overall_resistance),
    median = median(overall_resistance),
    IQR = IQR(overall_resistance),
    min = min(overall_resistance),
    max = max(overall_resistance)
  )
write_csv(resistance_summary, file.path(OUT, "vanbeeck_primary_resistance_summary.csv"))

# -----------------------------------------------------------------------------
# 5. Baseline family architecture
# -----------------------------------------------------------------------------
baseline_meta <- meta %>% filter(timepoint == "Baseline") %>% arrange(core)
stopifnot(nrow(baseline_meta) == 60, n_distinct(baseline_meta$core) == 60)
baseline_ra <- fam_ra[baseline_meta$sample_name, , drop = FALSE]

# Prevalence >=10% among baseline trajectories.
prevalence <- colMeans(baseline_ra > 0)
family_ids <- colnames(baseline_ra)
stopifnot(length(family_ids) == length(prevalence))
keep <- family_ids[prevalence >= 0.10]
write_csv(
  tibble(
    family = family_ids,
    prevalence = as.numeric(prevalence),
    retained = family_ids %in% keep
  ),
  file.path(OUT, "vanbeeck_baseline_family_prevalence.csv")
)
baseline_ra <- baseline_ra[, keep, drop = FALSE]
stopifnot(ncol(baseline_ra) >= 2)

# Zero replacement: half the smallest positive relative abundance, zeros only.
positive <- baseline_ra[baseline_ra > 0]
pseudo <- 0.5 * min(positive)
baseline_replaced <- baseline_ra
baseline_replaced[baseline_replaced == 0] <- pseudo
baseline_replaced <- baseline_replaced / rowSums(baseline_replaced)

# CLR transform.
log_base <- log(baseline_replaced)
clr <- log_base - rowMeans(log_base)
rownames(clr) <- rownames(baseline_replaced)
colnames(clr) <- colnames(baseline_replaced)
stopifnot(all(is.finite(clr)))

# Center CLR features within dairy before PCA. Do NOT center by treatment.
clr_centered <- clr
for (dd in levels(droplevels(baseline_meta$dairy))) {
  ii <- which(baseline_meta$dairy == dd)
  clr_centered[ii, ] <- sweep(
    clr[ii, , drop = FALSE],
    2,
    colMeans(clr[ii, , drop = FALSE]),
    "-"
  )
}
feature_sd <- apply(clr_centered, 2, sd)
clr_centered <- clr_centered[, is.finite(feature_sd) & feature_sd > 0, drop = FALSE]
stopifnot(ncol(clr_centered) >= 2)

# Independent Van Beeck PCA: this validates the phenomenon, not fixed axis transport.
pca <- prcomp(clr_centered, center = TRUE, scale. = FALSE)
variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pca_variance <- tibble(
  component = paste0("PC", seq_along(variance_explained)),
  variance_explained = variance_explained
)

scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  rownames_to_column("sample_name") %>%
  left_join(
    baseline_meta %>% select(sample_name, core, dairy, treatment),
    by = "sample_name"
  ) %>%
  mutate(
    z_PC1 = as.numeric(scale(PC1)),
    z_PC2 = as.numeric(scale(PC2))
  )

loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE]) %>%
  rownames_to_column("family")

write_csv(pca_variance, file.path(OUT, "vanbeeck_PCA_variance_explained.csv"))
write_csv(scores, file.path(OUT, "vanbeeck_baseline_PCA_scores.csv"))
write_csv(loadings, file.path(OUT, "vanbeeck_baseline_PCA_loadings.csv"))
write_csv(
  loadings %>% mutate(abs_loading = abs(PC1)) %>% arrange(desc(abs_loading)) %>% slice_head(n = 25),
  file.path(OUT, "vanbeeck_PC1_top25_family_loadings.csv")
)
write_csv(
  loadings %>% mutate(abs_loading = abs(PC2)) %>% arrange(desc(abs_loading)) %>% slice_head(n = 25),
  file.path(OUT, "vanbeeck_PC2_top25_family_loadings.csv")
)

cat("Families retained >=10%:", length(keep), "\n")
cat("CLR pseudocount:", format(pseudo, digits = 10), "\n")
cat("Families entering PCA:", ncol(clr_centered), "\n")
cat("PC1 variance explained:", variance_explained[1], "\n")
cat("PC2 variance explained:", variance_explained[2], "\n")

# -----------------------------------------------------------------------------
# 6. Frozen primary external-validation model
# -----------------------------------------------------------------------------
analysis <- trajectory_resistance %>%
  left_join(scores %>% select(core, z_PC1, z_PC2), by = "core") %>%
  mutate(dairy = droplevels(dairy), treatment = droplevels(treatment))
stopifnot(nrow(analysis) == 60, !anyNA(analysis$z_PC1), !anyNA(analysis$z_PC2))
write_csv(analysis, file.path(OUT, "vanbeeck_external_validation_analysis_table.csv"))

fit_reduced <- lm(overall_resistance ~ treatment + dairy, data = analysis)
fit_full <- lm(overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy, data = analysis)
stopifnot(!anyNA(coef(fit_full)))

# Classical nested partial F.
nested <- anova(fit_reduced, fit_full)
classical_F <- nested$F[2]
classical_p <- nested$`Pr(>F)`[2]
classical_df_num <- nested$Df[2]
classical_df_den <- df.residual(fit_full)

# Partial R2 for adding PC1 + PC2 beyond treatment + dairy.
rss0 <- sum(residuals(fit_reduced)^2)
rss1 <- sum(residuals(fit_full)^2)
partial_R2 <- (rss0 - rss1) / rss0

# HC3 coefficients and joint Wald test.
V3 <- vcovHC(fit_full, type = "HC3")
beta <- coef(fit_full)
se3 <- sqrt(diag(V3))
t3 <- beta / se3
dfr <- df.residual(fit_full)
p3 <- 2 * pt(abs(t3), df = dfr, lower.tail = FALSE)
HC3_coefficients <- tibble(
  term = names(beta),
  estimate = as.numeric(beta),
  HC3_SE = as.numeric(se3),
  t = as.numeric(t3),
  df = dfr,
  p = as.numeric(p3)
)

idx <- match(c("z_PC1", "z_PC2"), names(beta))
stopifnot(!anyNA(idx))
bpc <- beta[idx]
Vpc <- V3[idx, idx, drop = FALSE]
HC3_Wald <- as.numeric(t(bpc) %*% solve(Vpc) %*% bpc)
HC3_joint_p <- pchisq(HC3_Wald, df = 2, lower.tail = FALSE)

# -----------------------------------------------------------------------------
# 7. 9,999 Freedman-Lane permutations within dairy x treatment
# -----------------------------------------------------------------------------
N_PERM <- 9999
perm_strata <- interaction(analysis$dairy, analysis$treatment, drop = TRUE, lex.order = TRUE)
yhat0 <- fitted(fit_reduced)
e0 <- residuals(fit_reduced)
perm_F <- rep(NA_real_, N_PERM)

for (b in seq_len(N_PERM)) {
  ep <- e0
  for (ss in levels(perm_strata)) {
    ii <- which(perm_strata == ss)
    if (length(ii) > 1) ep[ii] <- sample(e0[ii], length(ii), replace = FALSE)
  }
  pd <- analysis
  pd$.perm_y <- yhat0 + ep
  pf0 <- lm(.perm_y ~ treatment + dairy, data = pd)
  pf1 <- lm(.perm_y ~ z_PC1 + z_PC2 + treatment + dairy, data = pd)
  perm_F[b] <- anova(pf0, pf1)$F[2]
}
permutation_p <- (1 + sum(perm_F >= classical_F, na.rm = TRUE)) / (N_PERM + 1)

primary_result <- tibble(
  study = "Van Beeck",
  role = "Independent external validation",
  taxonomic_level = "Family",
  n = nrow(analysis),
  outcome = "overall_resistance",
  model = "overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy",
  R2 = summary(fit_full)$r.squared,
  adjusted_R2 = summary(fit_full)$adj.r.squared,
  partial_R2_PC1_PC2 = partial_R2,
  classical_F = classical_F,
  classical_df_num = classical_df_num,
  classical_df_den = classical_df_den,
  classical_p = classical_p,
  HC3_Wald_chisq = HC3_Wald,
  HC3_df = 2,
  HC3_joint_p = HC3_joint_p,
  permutation_n = N_PERM,
  permutation_p = permutation_p
)

write_csv(HC3_coefficients, file.path(OUT, "PRIMARY_HC3_coefficients.csv"))
write_csv(primary_result, file.path(OUT, "PRIMARY_external_validation_joint_PC1_PC2.csv"))
write_csv(tibble(replicate = seq_len(N_PERM), partial_F = perm_F), file.path(OUT, "PERMUTATION_9999_partial_F.csv"))
capture.output(summary(fit_full), file = file.path(OUT, "PRIMARY_OLS_model_summary.txt"))

# -----------------------------------------------------------------------------
# 8. 5,000 complete-trajectory stratified bootstrap replicates
# -----------------------------------------------------------------------------
# PCA scores are frozen. This quantifies coefficient/effect-size uncertainty
# conditional on the independently learned Van Beeck ecological coordinate system.
N_BOOT <- 5000
boot_strata <- interaction(analysis$dairy, analysis$treatment, drop = TRUE, lex.order = TRUE)
stratum_indices <- split(seq_len(nrow(analysis)), boot_strata)

bootstrap_one <- function(b) {
  ids <- unlist(lapply(stratum_indices, function(ii) sample(ii, length(ii), replace = TRUE)), use.names = FALSE)
  z <- droplevels(analysis[ids, , drop = FALSE])
  f0 <- try(lm(overall_resistance ~ treatment + dairy, data = z), silent = TRUE)
  f1 <- try(lm(overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy, data = z), silent = TRUE)
  if (inherits(f0, "try-error") || inherits(f1, "try-error")) {
    return(tibble(replicate = b, beta_PC1 = NA_real_, beta_PC2 = NA_real_, partial_R2_PC1_PC2 = NA_real_, R2 = NA_real_))
  }
  bb <- coef(f1)
  if (!all(c("z_PC1", "z_PC2") %in% names(bb)) || anyNA(bb[c("z_PC1", "z_PC2")])) {
    return(tibble(replicate = b, beta_PC1 = NA_real_, beta_PC2 = NA_real_, partial_R2_PC1_PC2 = NA_real_, R2 = NA_real_))
  }
  r0 <- sum(residuals(f0)^2)
  r1 <- sum(residuals(f1)^2)
  tibble(
    replicate = b,
    beta_PC1 = unname(bb["z_PC1"]),
    beta_PC2 = unname(bb["z_PC2"]),
    partial_R2_PC1_PC2 = ifelse(r0 > 0, (r0 - r1) / r0, NA_real_),
    R2 = summary(f1)$r.squared
  )
}

cat("Running", N_BOOT, "stratified bootstrap replicates...\n")
boot <- bind_rows(lapply(seq_len(N_BOOT), bootstrap_one))
valid <- boot %>% filter(complete.cases(beta_PC1, beta_PC2, partial_R2_PC1_PC2, R2))
stopifnot(nrow(valid) > 0)

boot_summary_one <- function(x, original, statistic) {
  tibble(
    statistic = statistic,
    original = as.numeric(original),
    bootstrap_mean = mean(x),
    bootstrap_median = median(x),
    bootstrap_sd = sd(x),
    CI_2.5 = quantile(x, 0.025, names = FALSE),
    CI_97.5 = quantile(x, 0.975, names = FALSE),
    proportion_negative = mean(x < 0),
    proportion_positive = mean(x > 0),
    valid_replicates = length(x),
    requested_replicates = N_BOOT
  )
}

boot_summary <- bind_rows(
  boot_summary_one(valid$beta_PC1, beta["z_PC1"], "beta_PC1"),
  boot_summary_one(valid$beta_PC2, beta["z_PC2"], "beta_PC2"),
  boot_summary_one(valid$partial_R2_PC1_PC2, partial_R2, "partial_R2_PC1_PC2")
)

write_csv(boot, file.path(OUT, "BOOTSTRAP_5000_replicates.csv"))
write_csv(boot_summary, file.path(OUT, "BOOTSTRAP_5000_summary.csv"))

# -----------------------------------------------------------------------------
# 9. Frozen analysis record
# -----------------------------------------------------------------------------
freeze <- tibble(
  field = c(
    "study", "role", "primary_cohort", "primary_n", "taxonomic_level",
    "resistance_definition", "primary_outcome", "baseline_feature_space",
    "prevalence_filter", "zero_handling", "ordination", "ordination_centering",
    "primary_model", "primary_hypothesis", "design_covariates", "robust_inference",
    "permutation", "bootstrap", "subgroup_inference", "interaction_inference",
    "cross_study_PC_sign_comparison"
  ),
  value = c(
    "Van Beeck et al.; PRJEB63336",
    "Independent external validation of Italy + Manitoba ecological-resistance hypothesis",
    "All locked complete milk trajectories", "60", "Family",
    "R = 1 - Bray-Curtis(baseline, follow-up)",
    "Mean of Baseline->7 Days and Baseline->55-75DIM resistance",
    "Baseline family relative abundance -> CLR",
    ">=10% of baseline trajectories",
    "Half the smallest positive baseline relative abundance; zeros only; re-closure before CLR",
    "Independent Van Beeck PCA",
    "CLR features centered within dairy before PCA",
    "overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy",
    "Joint PC1 + PC2", "Treatment + dairy", "HC3",
    "9999 Freedman-Lane permutations within dairy x treatment",
    "5000 complete-trajectory resamples within dairy x treatment; frozen PCA scores",
    "None", "None",
    "Not permitted: independent PCA supports phenomenon-level replication, not fixed-axis transport"
  )
)
write_csv(freeze, file.path(OUT, "VANBEECK_EXTERNAL_VALIDATION_FREEZE.csv"))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

# -----------------------------------------------------------------------------
# 10. Terminal summary
# -----------------------------------------------------------------------------
cat("\n====================================\n")
cat("VAN BEECK EXTERNAL VALIDATION\n")
cat("====================================\n")
cat("N:", nrow(analysis), "\n")
cat("Taxonomic level: FAMILY\n")
cat("Primary model: overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy\n\n")
cat("HC3 COEFFICIENTS\n")
print(HC3_coefficients, width = Inf)
cat("\nJOINT PC1 + PC2 EXTERNAL-VALIDATION TEST\n")
print(primary_result, width = Inf)
cat("\nBOOTSTRAP UNCERTAINTY\n")
print(boot_summary, width = Inf)
cat("\n[PASS] VAN BEECK EXTERNAL VALIDATION COMPLETE\n")
cat("Output:", OUT, "\n")
