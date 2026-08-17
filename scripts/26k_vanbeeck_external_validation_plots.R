#!/usr/bin/env Rscript

# Publication-ready figures for the canonical Van Beeck external-validation run.
# Reads only frozen outputs from scripts/26j_vanbeeck_external_validation_only.R.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"
IN <- file.path(ROOT, "results/vanbeeck/external_validation_only")
OUT <- file.path(ROOT, "figures/vanbeeck_external_validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

analysis <- read_csv(file.path(IN, "vanbeeck_external_validation_analysis_table.csv"), show_col_types = FALSE) %>%
  mutate(
    dairy = factor(dairy),
    treatment = factor(treatment, levels = c("control low SCC", "control high SCC", "CB", "CH"))
  )
primary <- read_csv(file.path(IN, "PRIMARY_external_validation_joint_PC1_PC2.csv"), show_col_types = FALSE)
boot <- read_csv(file.path(IN, "BOOTSTRAP_5000_replicates.csv"), show_col_types = FALSE)
boot_summary <- read_csv(file.path(IN, "BOOTSTRAP_5000_summary.csv"), show_col_types = FALSE)
pca_scores <- read_csv(file.path(IN, "vanbeeck_baseline_PCA_scores.csv"), show_col_types = FALSE)
pca_var <- read_csv(file.path(IN, "vanbeeck_PCA_variance_explained.csv"), show_col_types = FALSE)

stopifnot(nrow(analysis) == 60, nrow(primary) == 1)

save_plot <- function(p, stem, width = 7.2, height = 5.2) {
  ggsave(file.path(OUT, paste0(stem, ".png")), p, width = width, height = height, dpi = 600)
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = width, height = height)
}

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, hjust = 0),
    legend.title = element_text(face = "bold")
  )

# Figure VB1. Distribution of the frozen primary phenotype.
p1 <- ggplot(analysis, aes(x = overall_resistance)) +
  geom_histogram(bins = 12, boundary = 0, closed = "left") +
  geom_vline(xintercept = mean(analysis$overall_resistance), linetype = 2) +
  labs(
    title = "Van Beeck ecological resistance",
    subtitle = "All 60 locked longitudinal milk trajectories",
    x = "Overall ecological resistance",
    y = "Number of trajectories",
    caption = "Overall resistance = mean[1 - Bray(Baseline, 7 Days), 1 - Bray(Baseline, 55–75 DIM)]. Dashed line: cohort mean."
  ) +
  base_theme
save_plot(p1, "Fig_VB1_overall_resistance_distribution", 7.2, 4.8)

# Figures VB2A/B. Added-variable visualizations of the frozen model coefficients.
make_added_variable <- function(pc, other_pc, label) {
  y_formula <- as.formula(paste("overall_resistance ~", other_pc, "+ treatment + dairy"))
  x_formula <- as.formula(paste(pc, "~", other_pc, "+ treatment + dairy"))
  z <- analysis %>%
    mutate(
      outcome_residual = resid(lm(y_formula, data = analysis)),
      pc_residual = resid(lm(x_formula, data = analysis))
    )

  ggplot(z, aes(x = pc_residual, y = outcome_residual)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = paste0("Adjusted ", label, " association with ecological resistance"),
      subtitle = "Added-variable visualization of the frozen external-validation model",
      x = paste0(label, " residual after adjustment"),
      y = "Resistance residual after adjustment",
      caption = "Adjustment includes the other PC, treatment group and dairy. The inferential target remains the joint PC1 + PC2 test."
    ) +
    base_theme
}

p2a <- make_added_variable("z_PC1", "z_PC2", "PC1")
p2b <- make_added_variable("z_PC2", "z_PC1", "PC2")
save_plot(p2a, "Fig_VB2A_PC1_added_variable", 7.2, 5.2)
save_plot(p2b, "Fig_VB2B_PC2_added_variable", 7.2, 5.2)

# Figure VB3. Bootstrap coefficient intervals.
forest <- boot_summary %>%
  filter(statistic %in% c("beta_PC1", "beta_PC2")) %>%
  mutate(
    axis = recode(statistic, beta_PC1 = "PC1", beta_PC2 = "PC2"),
    axis = factor(axis, levels = c("PC2", "PC1"))
  )

p3 <- ggplot(forest, aes(y = axis, x = original)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = CI_2.5, xmax = CI_97.5),
    orientation = "y",
    width = 0.15
  ) +
  geom_point(size = 2.5) +
  labs(
    title = "Bootstrap stability of Van Beeck baseline-architecture coefficients",
    subtitle = "5,000 complete-trajectory resamples within dairy × treatment strata",
    x = "Coefficient for ecological resistance",
    y = NULL,
    caption = "Points are observed coefficients; horizontal bars are percentile 95% bootstrap intervals. PCA scores were frozen during resampling."
  ) +
  base_theme
save_plot(p3, "Fig_VB3_bootstrap_PC_coefficients", 7.2, 4.2)

# Figure VB4. Bootstrap distribution of the joint architecture effect size.
observed_pr2 <- primary$partial_R2_PC1_PC2[[1]]
p4 <- ggplot(boot, aes(x = partial_R2_PC1_PC2)) +
  geom_histogram(bins = 35) +
  geom_vline(xintercept = observed_pr2, linetype = 2) +
  labs(
    title = "Bootstrap uncertainty in the joint architecture effect size",
    subtitle = paste0("Observed partial R² = ", sprintf("%.3f", observed_pr2)),
    x = "Partial R² for adding PC1 + PC2 beyond treatment + dairy",
    y = "Bootstrap replicates",
    caption = "The partial-R² bootstrap is an effect-size uncertainty analysis, not an additional null-hypothesis test."
  ) +
  base_theme
save_plot(p4, "Fig_VB4_bootstrap_partial_R2", 7.2, 4.8)

# Figure VB5. Concordance of the three joint hypothesis tests.
inference <- tibble(
  method = c("Classical partial F", "HC3 joint Wald", "Freedman–Lane permutation"),
  p = c(primary$classical_p[[1]], primary$HC3_joint_p[[1]], primary$permutation_p[[1]])
) %>%
  mutate(
    minus_log10_p = -log10(p),
    method = factor(method, levels = rev(method))
  )

p5 <- ggplot(inference, aes(x = minus_log10_p, y = method)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2) +
  geom_segment(aes(x = 0, xend = minus_log10_p, yend = method), linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(aes(label = paste0("p = ", signif(p, 3))), hjust = -0.08, size = 3.4) +
  expand_limits(x = max(inference$minus_log10_p) * 1.25) +
  labs(
    title = "Concordant support for the joint PC1 + PC2 validation test",
    x = expression(-log[10](italic(p))),
    y = NULL,
    caption = "Dashed line marks p = 0.05. The permutation test preserves the dairy × treatment design strata."
  ) +
  base_theme
save_plot(p5, "Fig_VB5_joint_test_concordance", 7.2, 4.5)

# Figure VB6. Baseline PC space, shown only as ordination/QC context.
pc1_pct <- 100 * pca_var$variance_explained[pca_var$component == "PC1"][[1]]
pc2_pct <- 100 * pca_var$variance_explained[pca_var$component == "PC2"][[1]]
p6 <- ggplot(pca_scores, aes(x = PC1, y = PC2, shape = dairy)) +
  geom_point(size = 2.4, alpha = 0.8) +
  labs(
    title = "Van Beeck baseline family-level ecological architecture",
    subtitle = "Independent PCA after within-dairy centering of CLR family profiles",
    x = paste0("PC1 (", sprintf("%.1f", pc1_pct), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", pc2_pct), "%)"),
    shape = "Dairy",
    caption = "The Van Beeck axes are cohort-specific. PC signs/loadings must not be interpreted as identical to independently fitted Italy/Manitoba or Porcellato axes."
  ) +
  base_theme
save_plot(p6, "Fig_VB6_baseline_PCA", 7.2, 5.2)

cat("[PASS] Van Beeck external-validation plots written to:", OUT, "\n")
