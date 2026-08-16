#!/usr/bin/env Rscript

# ============================================================================
# MMER
#
# PC2 ECOLOGICAL-STATE INTERPRETATION
#
# DESCRIPTIVE ONLY
#
# Primary cohorts:
#   Italy
#   Manitoba
#
# Purpose:
# Characterize the baseline family-level ecological configuration represented
# by the frozen PC2 axis from the primary healthy two-aim analysis.
#
# NO NEW HYPOTHESIS TESTS ARE PERFORMED.
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"

PRIMARY_DIR <- file.path(
  ROOT,
  "results/primary_healthy_twoAims"
)

FAMILY_FILE <- file.path(
  ROOT,
  "results/harmonized/three_study_family_long.csv"
)

ANALYSIS_FILE <- file.path(
  PRIMARY_DIR,
  "tables/PRIMARY_analysis_table.csv"
)

LOADINGS_FILE <- file.path(
  PRIMARY_DIR,
  "tables/baseline_PC_family_loadings.csv"
)

OUT <- file.path(
  ROOT,
  "results/PC2_ecological_state_interpretation"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(OUT, "tables"),
  showWarnings = FALSE
)

dir.create(
  file.path(OUT, "figures"),
  showWarnings = FALSE
)

PRIMARY_STUDIES <- c(
  "Italy",
  "Manitoba"
)

TOP_N <- 10


# ============================================================================
# READ
# ============================================================================

analysis <- read_csv(
  ANALYSIS_FILE,
  show_col_types = FALSE
)

loadings <- read_csv(
  LOADINGS_FILE,
  show_col_types = FALSE
)

fam <- read_csv(
  FAMILY_FILE,
  show_col_types = FALSE
)


cat("\n============================================================\n")
cat("PC2 ECOLOGICAL-STATE INTERPRETATION\n")
cat("============================================================\n")


# ============================================================================
# DEFINE PC2 TERTILES
#
# Global tertiles using the frozen z_PC2 scores.
# Study centering has already been performed upstream.
# ============================================================================

analysis <- analysis %>%
  mutate(
    PC2_tertile_num =
      ntile(
        z_PC2,
        3
      ),

    PC2_tertile =
      factor(
        PC2_tertile_num,
        levels = 1:3,
        labels = c(
          "Low_PC2",
          "Mid_PC2",
          "High_PC2"
        )
      )
  )


tertile_manifest <- analysis %>%
  select(
    study,
    trajectory_global,
    cow_global,
    cow_id,
    quarter,
    z_PC2,
    ecological_resistance,
    bray_displacement,
    PC2_tertile
  )


write_csv(
  tertile_manifest,
  file.path(
    OUT,
    "tables",
    "PC2_tertile_manifest.csv"
  )
)


cat("\nPC2 tertile distribution:\n")

print(
  analysis %>%
    count(
      study,
      PC2_tertile
    )
)


# ============================================================================
# IDENTIFY TOP POSITIVE / NEGATIVE PC2 FAMILIES
# ============================================================================

top_positive <- loadings %>%
  arrange(
    desc(
      PC2_loading
    )
  ) %>%
  slice_head(
    n = TOP_N
  ) %>%
  mutate(
    PC2_side =
      "Positive_PC2"
  )


top_negative <- loadings %>%
  arrange(
    PC2_loading
  ) %>%
  slice_head(
    n = TOP_N
  ) %>%
  mutate(
    PC2_side =
      "Negative_PC2"
  )


selected_families <- bind_rows(
  top_positive,
  top_negative
) %>%
  distinct(
    family,
    .keep_all = TRUE
  )


cat("\nSelected PC2 families:\n")
print(
  selected_families,
  n = Inf
)


write_csv(
  selected_families,
  file.path(
    OUT,
    "tables",
    "selected_PC2_families.csv"
  )
)


# ============================================================================
# BASELINE T1 FAMILY DATA
# ============================================================================

baseline <- fam %>%
  filter(
    study %in% PRIMARY_STUDIES,
    harmonized_timepoint == "T1",
    trajectory_global %in%
      analysis$trajectory_global,
    family %in%
      selected_families$family
  ) %>%

  group_by(
    study,
    trajectory_global,
    family
  ) %>%

  summarise(
    count_3k =
      mean(
        count_3k
      ),
    .groups = "drop"
  )


# ============================================================================
# COMPLETE ZERO-FILLED FAMILY TABLE
#
# Needed for correct prevalence / relative abundance summaries.
# ============================================================================

baseline_complete <- expand_grid(
  trajectory_global =
    analysis$trajectory_global,
  family =
    selected_families$family
) %>%

  left_join(
    analysis %>%
      select(
        trajectory_global,
        study,
        PC2_tertile,
        z_PC2,
        ecological_resistance
      ),
    by = "trajectory_global"
  ) %>%

  left_join(
    baseline,
    by = c(
      "study",
      "trajectory_global",
      "family"
    )
  ) %>%

  mutate(
    count_3k =
      replace_na(
        count_3k,
        0
      ),

    relative_abundance =
      count_3k /
      3000
  ) %>%

  left_join(
    selected_families %>%
      select(
        family,
        PC2_loading,
        PC2_side
      ),
    by = "family"
  )


# ============================================================================
# FAMILY × STUDY × TERTILE SUMMARY
# ============================================================================

family_summary <- baseline_complete %>%

  group_by(
    study,
    PC2_tertile,
    family,
    PC2_loading,
    PC2_side
  ) %>%

  summarise(
    n_quarters =
      n(),

    prevalence =
      mean(
        count_3k > 0
      ),

    median_relative_abundance =
      median(
        relative_abundance
      ),

    mean_relative_abundance =
      mean(
        relative_abundance
      ),

    q25_relative_abundance =
      quantile(
        relative_abundance,
        0.25
      ),

    q75_relative_abundance =
      quantile(
        relative_abundance,
        0.75
      ),

    .groups =
      "drop"
  )


write_csv(
  family_summary,
  file.path(
    OUT,
    "tables",
    "PC2_family_by_study_tertile_summary.csv"
  )
)


# ============================================================================
# POOLED TERTILE SUMMARY
#
# Descriptive across Italy + Manitoba.
# ============================================================================

pooled_family_summary <- baseline_complete %>%

  group_by(
    PC2_tertile,
    family,
    PC2_loading,
    PC2_side
  ) %>%

  summarise(
    n_quarters =
      n(),

    prevalence =
      mean(
        count_3k > 0
      ),

    median_relative_abundance =
      median(
        relative_abundance
      ),

    mean_relative_abundance =
      mean(
        relative_abundance
      ),

    .groups =
      "drop"
  )


write_csv(
  pooled_family_summary,
  file.path(
    OUT,
    "tables",
    "PC2_family_pooled_tertile_summary.csv"
  )
)


# ============================================================================
# HIGH vs LOW PC2 DESCRIPTIVE CONTRAST
# ============================================================================

high_low <- pooled_family_summary %>%

  filter(
    PC2_tertile %in%
      c(
        "Low_PC2",
        "High_PC2"
      )
  ) %>%

  select(
    PC2_tertile,
    family,
    PC2_loading,
    PC2_side,
    prevalence,
    median_relative_abundance,
    mean_relative_abundance
  ) %>%

  pivot_wider(
    names_from =
      PC2_tertile,

    values_from =
      c(
        prevalence,
        median_relative_abundance,
        mean_relative_abundance
      )
  ) %>%

  mutate(
    prevalence_difference =
      prevalence_High_PC2 -
      prevalence_Low_PC2,

    mean_RA_difference =
      mean_relative_abundance_High_PC2 -
      mean_relative_abundance_Low_PC2
  ) %>%

  arrange(
    desc(
      abs(
        mean_RA_difference
      )
    )
  )


write_csv(
  high_low,
  file.path(
    OUT,
    "tables",
    "PC2_high_vs_low_descriptive_contrast.csv"
  )
)


cat("\nHigh vs low PC2 descriptive contrast:\n")

print(
  high_low,
  n = Inf
)


# ============================================================================
# RESISTANCE SUMMARY BY PC2 TERTILE
# ============================================================================

resistance_summary <- analysis %>%

  group_by(
    study,
    PC2_tertile
  ) %>%

  summarise(
    n =
      n(),

    mean_resistance =
      mean(
        ecological_resistance
      ),

    sd_resistance =
      sd(
        ecological_resistance
      ),

    median_resistance =
      median(
        ecological_resistance
      ),

    mean_displacement =
      mean(
        bray_displacement
      ),

    .groups =
      "drop"
  )


write_csv(
  resistance_summary,
  file.path(
    OUT,
    "tables",
    "PC2_tertile_resistance_summary.csv"
  )
)


# ============================================================================
# FAMILY ORDER FOR FIGURES
#
# Negative PC2 at bottom -> positive PC2 at top.
# ============================================================================

family_order <- selected_families %>%
  arrange(
    PC2_loading
  ) %>%
  pull(
    family
  )


# ============================================================================
# FIGURE 1
# MEDIAN BASELINE RELATIVE ABUNDANCE HEATMAP
# ============================================================================

heat_ra <- family_summary %>%
  mutate(
    family =
      factor(
        family,
        levels = family_order
      ),

    group_label =
      paste(
        study,
        PC2_tertile,
        sep = "\n"
      ),

    median_RA_percent =
      100 *
      median_relative_abundance
  )


p1 <- ggplot(
  heat_ra,
  aes(
    x = group_label,
    y = family,
    fill = median_RA_percent
  )
) +

  geom_tile() +

  labs(
    x = NULL,
    y = NULL,
    fill =
      "Median baseline\nrelative abundance (%)",
    title =
      "Baseline family composition across the PC2 ecological gradient"
  ) +

  theme_minimal(
    base_size = 11
  ) +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


ggsave(
  file.path(
    OUT,
    "figures",
    "Figure_PC2_family_median_RA_heatmap.pdf"
  ),
  p1,
  width = 8,
  height = 7
)


# ============================================================================
# FIGURE 2
# PREVALENCE HEATMAP
# ============================================================================

heat_prev <- family_summary %>%
  mutate(
    family =
      factor(
        family,
        levels = family_order
      ),

    group_label =
      paste(
        study,
        PC2_tertile,
        sep = "\n"
      ),

    prevalence_percent =
      100 *
      prevalence
  )


p2 <- ggplot(
  heat_prev,
  aes(
    x = group_label,
    y = family,
    fill = prevalence_percent
  )
) +

  geom_tile() +

  labs(
    x = NULL,
    y = NULL,
    fill =
      "Prevalence (%)",
    title =
      "Baseline prevalence across the PC2 ecological gradient"
  ) +

  theme_minimal(
    base_size = 11
  ) +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


ggsave(
  file.path(
    OUT,
    "figures",
    "Figure_PC2_family_prevalence_heatmap.pdf"
  ),
  p2,
  width = 8,
  height = 7
)


# ============================================================================
# FIGURE 3
# RESISTANCE ACROSS PC2 TERTILES
#
# Descriptive visualization only.
# ============================================================================

p3 <- ggplot(
  analysis,
  aes(
    x = PC2_tertile,
    y = ecological_resistance,
    shape = study
  )
) +

  geom_boxplot(
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.10,
    alpha = 0.8,
    size = 2
  ) +

  labs(
    x =
      "Baseline PC2 ecological state",

    y =
      "Ecological resistance (1 - Bray-Curtis)",

    title =
      "Ecological resistance across the baseline PC2 gradient"
  ) +

  theme_classic(
    base_size = 12
  )


ggsave(
  file.path(
    OUT,
    "figures",
    "Figure_PC2_tertile_resistance.pdf"
  ),
  p3,
  width = 6.5,
  height = 4.5
)


# ============================================================================
# HUMAN-READABLE SUMMARY
# ============================================================================

sink(
  file.path(
    OUT,
    "PC2_ECOLOGICAL_STATE_SUMMARY.txt"
  )
)


cat(
  "PC2 ECOLOGICAL-STATE INTERPRETATION\n"
)

cat(
  "===================================\n\n"
)


cat(
  "This analysis is descriptive only.\n"
)

cat(
  "No new hypothesis tests were performed.\n\n"
)


cat(
  "PC2 RESISTANCE SUMMARY\n"
)

print(
  resistance_summary
)


cat(
  "\n\nTOP POSITIVE PC2 FAMILIES\n"
)

print(
  top_positive,
  n = Inf
)


cat(
  "\n\nTOP NEGATIVE PC2 FAMILIES\n"
)

print(
  top_negative,
  n = Inf
)


cat(
  "\n\nHIGH vs LOW PC2 DESCRIPTIVE CONTRAST\n"
)

print(
  high_low,
  n = Inf
)


sink()


capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "sessionInfo.txt"
  )
)


cat("\n============================================================\n")
cat("[PASS] PC2 ECOLOGICAL-STATE INTERPRETATION COMPLETE\n")
cat("============================================================\n")

cat(
  "Output directory:",
  OUT,
  "\n"
)
