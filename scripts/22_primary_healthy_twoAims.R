#!/usr/bin/env Rscript

# ============================================================================
# MMER
#
# PRIMARY HEALTHY MAMMARY-QUARTER ECOLOGICAL RESISTANCE ANALYSIS
#
# PRIMARY DATASETS ONLY:
#   Italy
#   Manitoba
#
# Wisconsin is intentionally excluded.
#
# ---------------------------------------------------------------------------
# CENTRAL QUESTION
#
# Does the ecological state of a healthy mammary quarter at dry-off predict
# how strongly that same quarter changes during the transition into early
# postpartum?
#
# ---------------------------------------------------------------------------
# AIM 1 — ECOLOGICAL SUSCEPTIBILITY
#
# Test whether baseline multivariate ecological architecture predicts
# quarter-level ecological resistance across Italy + Manitoba.
#
# Primary outcome:
#
#     resistance = 1 - Bray-Curtis(T1, T2)
#
# Larger values = greater compositional retention / resistance.
# Smaller values = greater ecological displacement / susceptibility.
#
# Baseline architecture:
#
#     T1 family composition
#       -> prevalence filter
#       -> CLR
#       -> within-study centering
#       -> pooled PCA
#       -> PC1 + PC2
#
# Primary inferential test:
#
#     H0: beta_PC1 = beta_PC2 = 0
#
# Cow-clustered CR2 small-sample robust inference.
#
# ---------------------------------------------------------------------------
# AIM 2 — CROSS-COHORT REPRODUCIBILITY
#
# Test whether the baseline architecture -> resistance relationship differs
# between Italy and Manitoba.
#
#     H0:
#       PC1 x Study = 0
#       PC2 x Study = 0
#
# Again using cow-clustered CR2 inference.
#
# ---------------------------------------------------------------------------
# TAXONOMIC INTERPRETATION
#
# PC loadings are descriptive interpretation of the prespecified multivariate
# ecological axes.
#
# No family-by-family hypothesis fishing is performed.
# ============================================================================


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(clubSandwich)
})


set.seed(123)


# ============================================================================
# CONFIGURATION
# ============================================================================

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"

FAMILY_FILE <- file.path(
  ROOT,
  "results/harmonized/three_study_family_long.csv"
)

ELIG_FILE <- file.path(
  ROOT,
  "results/harmonized/three_study_trajectory_eligibility.csv"
)

OUT <- file.path(
  ROOT,
  "results/primary_healthy_twoAims"
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

dir.create(
  file.path(OUT, "models"),
  showWarnings = FALSE
)

dir.create(
  file.path(OUT, "qc"),
  showWarnings = FALSE
)


PRIMARY_STUDIES <- c(
  "Italy",
  "Manitoba"
)

PSEUDOCOUNT <- 0.5

PREVALENCE_THRESHOLD <- 0.10


# ============================================================================
# HELPERS
# ============================================================================

bray_pair <- function(a, b) {

  den <- sum(a + b)

  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }

  sum(
    abs(a - b)
  ) / den
}


clr_transform <- function(
  M,
  pseudocount = 0.5
) {

  L <- log(
    M + pseudocount
  )

  sweep(
    L,
    1,
    rowMeans(L),
    "-"
  )
}


is_organelle_family <- function(x) {

  grepl(
    "mitochond|chloroplast|plastid|organelle",
    x,
    ignore.case = TRUE
  )
}


cr2_joint <- function(
  fit,
  cluster,
  terms,
  label
) {

  coefficient_names <- names(
    coef(fit)
  )

  missing_terms <- setdiff(
    terms,
    coefficient_names
  )


  if (length(missing_terms) > 0) {

    stop(
      "Missing model coefficient(s): ",
      paste(
        missing_terms,
        collapse = ", "
      )
    )
  }


  V <- clubSandwich::vcovCR(
    fit,
    cluster = cluster,
    type = "CR2"
  )


  C <- clubSandwich::constrain_zero(
    terms,
    coefs = coef(fit)
  )


  W <- clubSandwich::Wald_test(
    fit,
    constraints = C,
    vcov = V,
    test = "HTZ"
  )


  z <- as.data.frame(W)


  tibble(
    hypothesis = label,
    Fstat = z$Fstat,
    df_num = z$df_num,
    df_denom = z$df_denom,
    p_value = z$p_val
  )
}


cr2_coefficients <- function(
  fit,
  cluster
) {

  V <- clubSandwich::vcovCR(
    fit,
    cluster = cluster,
    type = "CR2"
  )


  z <- clubSandwich::coef_test(
    fit,
    vcov = V,
    test = "Satterthwaite"
  )


  as.data.frame(z) %>%
    rownames_to_column(
      "term"
    ) %>%
    as_tibble()
}


# ============================================================================
# READ DATA
# ============================================================================

fam <- read_csv(
  FAMILY_FILE,
  show_col_types = FALSE
)


elig <- read_csv(
  ELIG_FILE,
  show_col_types = FALSE
)


cat("\n")
cat("============================================================\n")
cat("MMER — PRIMARY HEALTHY TWO-AIM ANALYSIS\n")
cat("============================================================\n")

cat(
  "Studies:",
  paste(
    PRIMARY_STUDIES,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Wisconsin: EXCLUDED\n\n"
)


# ============================================================================
# REMOVE ORGANELLE LABELS
# ============================================================================

organelle_rows <- fam %>%
  filter(
    is_organelle_family(
      family
    )
  )


write_csv(
  organelle_rows %>%
    distinct(
      study,
      family
    ),
  file.path(
    OUT,
    "qc",
    "removed_organelle_families.csv"
  )
)


cat(
  "Organelle-labelled rows removed:",
  nrow(organelle_rows),
  "\n"
)


fam <- fam %>%
  filter(
    !is_organelle_family(
      family
    )
  )


# ============================================================================
# RESTRICT TO PRIMARY HEALTHY STUDIES
# ============================================================================

fam <- fam %>%
  filter(
    study %in% PRIMARY_STUDIES
  )


elig_primary <- elig %>%
  filter(
    study %in% PRIMARY_STUDIES,
    primary_complete
  )


primary_ids <- unique(
  elig_primary$trajectory_global
)


fam <- fam %>%
  filter(
    trajectory_global %in%
      primary_ids
  )


# ============================================================================
# BUILD T1 / T2 QUARTER PROFILES
# ============================================================================

profile_long <- fam %>%
  filter(
    harmonized_timepoint %in%
      c(
        "T1",
        "T2"
      )
  ) %>%

  group_by(
    study,
    trajectory_global,
    cow_global,
    cow_id,
    quarter,
    harmonized_timepoint,
    family
  ) %>%

  summarise(
    count = mean(
      count_3k
    ),
    .groups = "drop"
  )


# Verify both harmonized states exist

complete <- profile_long %>%
  distinct(
    study,
    trajectory_global,
    harmonized_timepoint
  ) %>%

  count(
    study,
    trajectory_global,
    name = "n_timepoints"
  ) %>%

  filter(
    n_timepoints == 2
  )


profile_long <- profile_long %>%
  semi_join(
    complete,
    by = c(
      "study",
      "trajectory_global"
    )
  )


traj_meta <- profile_long %>%
  distinct(
    study,
    trajectory_global,
    cow_global,
    cow_id,
    quarter
  )


sample_summary <- traj_meta %>%
  group_by(
    study
  ) %>%

  summarise(
    trajectories = n(),
    cows = n_distinct(
      cow_global
    ),
    .groups = "drop"
  )


cat("\nAnalysis population:\n")
print(sample_summary)


write_csv(
  sample_summary,
  file.path(
    OUT,
    "tables",
    "analysis_population.csv"
  )
)


write_csv(
  traj_meta,
  file.path(
    OUT,
    "tables",
    "quarter_trajectory_manifest.csv"
  )
)


# ============================================================================
# CONSTRUCT FAMILY MATRICES
# ============================================================================

family_universe <- sort(
  unique(
    profile_long$family
  )
)


make_matrix <- function(tp) {

  w <- profile_long %>%

    filter(
      harmonized_timepoint == tp
    ) %>%

    select(
      trajectory_global,
      family,
      count
    ) %>%

    pivot_wider(
      names_from = family,
      values_from = count,
      values_fill = 0
    )


  ids <- w$trajectory_global


  M <- as.matrix(
    w[
      ,
      setdiff(
        names(w),
        "trajectory_global"
      ),
      drop = FALSE
    ]
  )


  rownames(M) <- ids


  missing_families <- setdiff(
    family_universe,
    colnames(M)
  )


  if (length(missing_families) > 0) {

    M <- cbind(
      M,
      matrix(
        0,
        nrow = nrow(M),
        ncol = length(missing_families),
        dimnames = list(
          NULL,
          missing_families
        )
      )
    )
  }


  M[
    ,
    family_universe,
    drop = FALSE
  ]
}


T1 <- make_matrix(
  "T1"
)

T2 <- make_matrix(
  "T2"
)


common_ids <- intersect(
  rownames(T1),
  rownames(T2)
)


T1 <- T1[
  common_ids,
  ,
  drop = FALSE
]


T2 <- T2[
  common_ids,
  ,
  drop = FALSE
]


traj_meta <- traj_meta %>%
  filter(
    trajectory_global %in%
      common_ids
  )


# ============================================================================
# DEFINE ECOLOGICAL RESPONSE PHENOTYPE
# ============================================================================

phenotype <- tibble(

  trajectory_global =
    common_ids,

  bray_displacement =
    vapply(
      common_ids,
      function(id) {

        bray_pair(
          T1[id, ],
          T2[id, ]
        )
      },
      numeric(1)
    )

) %>%

  mutate(

    ecological_resistance =
      1 -
      bray_displacement

  ) %>%

  left_join(
    traj_meta,
    by = "trajectory_global"
  )


if (
  any(
    !is.finite(
      phenotype$ecological_resistance
    )
  )
) {

  stop(
    "Non-finite ecological resistance values detected."
  )
}


write_csv(
  phenotype,
  file.path(
    OUT,
    "tables",
    "quarter_ecological_resistance.csv"
  )
)


phenotype_summary <- phenotype %>%

  group_by(
    study
  ) %>%

  summarise(

    trajectories =
      n(),

    cows =
      n_distinct(
        cow_global
      ),

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

    IQR_resistance =
      IQR(
        ecological_resistance
      ),

    mean_displacement =
      mean(
        bray_displacement
      ),

    .groups =
      "drop"
  )


cat("\nEcological resistance phenotype:\n")
print(phenotype_summary)


write_csv(
  phenotype_summary,
  file.path(
    OUT,
    "tables",
    "ecological_resistance_summary.csv"
  )
)


# ============================================================================
# BASELINE ECOLOGICAL ARCHITECTURE
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("BASELINE ECOLOGICAL ARCHITECTURE\n")
cat("============================================================\n")


# ---------------------------------------------------------------------------
# Prevalence filter based ONLY on baseline T1
# ---------------------------------------------------------------------------

prevalence <- colMeans(
  T1 > 0
)


prevalence_table <- tibble(

  family =
    names(prevalence),

  prevalence =
    as.numeric(
      prevalence
    ),

  retained =
    prevalence >=
    PREVALENCE_THRESHOLD

)


write_csv(
  prevalence_table,
  file.path(
    OUT,
    "tables",
    "baseline_family_prevalence.csv"
  )
)


keep_family <- prevalence >=
  PREVALENCE_THRESHOLD


X <- T1[
  ,
  keep_family,
  drop = FALSE
]


cat(
  "Baseline families before filter:",
  ncol(T1),
  "\n"
)

cat(
  "Baseline families retained:",
  ncol(X),
  "\n"
)


# ============================================================================
# CLR TRANSFORMATION
# ============================================================================

Xclr <- clr_transform(
  X,
  pseudocount =
    PSEUDOCOUNT
)


# ============================================================================
# REMOVE STUDY CENTROIDS BEFORE PCA
#
# This prevents the PCA from merely rediscovering Italy vs Manitoba.
# ============================================================================

meta_ordered <- phenotype[
  match(
    rownames(Xclr),
    phenotype$trajectory_global
  ),
]


if (
  !all(
    meta_ordered$trajectory_global ==
      rownames(Xclr)
  )
) {

  stop(
    "Metadata / matrix row-order mismatch."
  )
}


Xcenter <- Xclr


centroid_table <- list()


for (s in PRIMARY_STUDIES) {

  idx <- which(
    meta_ordered$study == s
  )


  centroid <- colMeans(
    Xclr[
      idx,
      ,
      drop = FALSE
    ]
  )


  centroid_table[[s]] <- tibble(

    study = s,

    family =
      colnames(Xclr),

    centroid =
      centroid

  )


  Xcenter[
    idx,
  ] <- sweep(
    Xclr[
      idx,
      ,
      drop = FALSE
    ],
    2,
    centroid,
    "-"
  )
}


centroid_table <- bind_rows(
  centroid_table
)


write_csv(
  centroid_table,
  file.path(
    OUT,
    "tables",
    "within_study_CLR_centroids.csv"
  )
)


# ---------------------------------------------------------------------------
# Numerical QC
# ---------------------------------------------------------------------------

centering_QC <- bind_rows(

  lapply(
    PRIMARY_STUDIES,
    function(s) {

      idx <- which(
        meta_ordered$study == s
      )


      means <- colMeans(
        Xcenter[
          idx,
          ,
          drop = FALSE
        ]
      )


      tibble(

        study = s,

        max_absolute_centered_family_mean =
          max(
            abs(
              means
            )
          )
      )
    }
  )
)


cat("\nWithin-study centering QC:\n")
print(centering_QC)


write_csv(
  centering_QC,
  file.path(
    OUT,
    "qc",
    "within_study_centering_QC.csv"
  )
)


# ============================================================================
# REMOVE ZERO-VARIANCE FEATURES
# ============================================================================

feature_sd <- apply(
  Xcenter,
  2,
  sd
)


keep_variable <- is.finite(
  feature_sd
) &
  feature_sd > 0


Xcenter <- Xcenter[
  ,
  keep_variable,
  drop = FALSE
]


cat(
  "Families entering PCA:",
  ncol(Xcenter),
  "\n"
)


# ============================================================================
# PCA
# ============================================================================

pca <- prcomp(
  Xcenter,
  center = FALSE,
  scale. = FALSE
)


saveRDS(
  pca,
  file.path(
    OUT,
    "models",
    "baseline_withinStudy_CLR_PCA.rds"
  )
)


variance <- pca$sdev^2 /
  sum(
    pca$sdev^2
  )


variance_table <- tibble(

  PC =
    paste0(
      "PC",
      seq_along(
        variance
      )
    ),

  variance_explained =
    variance,

  cumulative_variance =
    cumsum(
      variance
    )

)


write_csv(
  variance_table,
  file.path(
    OUT,
    "tables",
    "PCA_variance_explained.csv"
  )
)


cat(
  "\nPC1 variance:",
  round(
    variance[1] * 100,
    2
  ),
  "%\n"
)

cat(
  "PC2 variance:",
  round(
    variance[2] * 100,
    2
  ),
  "%\n"
)

cat(
  "PC1 + PC2:",
  round(
    sum(
      variance[1:2]
    ) * 100,
    2
  ),
  "%\n"
)


# ============================================================================
# PCA SCORES
# ============================================================================

scores <- tibble(

  trajectory_global =
    rownames(
      pca$x
    ),

  PC1 =
    pca$x[
      ,
      1
    ],

  PC2 =
    pca$x[
      ,
      2
    ]

)


# Standardize axes for interpretable beta per 1 SD

scores <- scores %>%
  mutate(

    z_PC1 =
      as.numeric(
        scale(
          PC1
        )
      ),

    z_PC2 =
      as.numeric(
        scale(
          PC2
        )
      )
  )


write_csv(
  scores,
  file.path(
    OUT,
    "tables",
    "baseline_PC_scores.csv"
  )
)


# ============================================================================
# PCA LOADINGS
# ============================================================================

loadings <- tibble(

  family =
    rownames(
      pca$rotation
    ),

  PC1_loading =
    pca$rotation[
      ,
      1
    ],

  PC2_loading =
    pca$rotation[
      ,
      2
    ]

)


write_csv(
  loadings,
  file.path(
    OUT,
    "tables",
    "baseline_PC_family_loadings.csv"
  )
)


# ============================================================================
# ANALYSIS TABLE
# ============================================================================

analysis <- phenotype %>%

  left_join(
    scores,
    by = "trajectory_global"
  ) %>%

  mutate(

    study =
      factor(
        study,
        levels =
          PRIMARY_STUDIES
      ),

    cow_global =
      factor(
        cow_global
      )
  )


if (
  any(
    !complete.cases(
      analysis[
        ,
        c(
          "ecological_resistance",
          "z_PC1",
          "z_PC2",
          "study",
          "cow_global"
        )
      ]
    )
  )
) {

  stop(
    "Missing values detected in primary analysis variables."
  )
}


write_csv(
  analysis,
  file.path(
    OUT,
    "tables",
    "PRIMARY_analysis_table.csv"
  )
)


cat("\nFinal primary sample:\n")

print(
  analysis %>%
    group_by(
      study
    ) %>%
    summarise(
      trajectories =
        n(),
      cows =
        n_distinct(
          cow_global
        ),
      .groups =
        "drop"
    )
)


# ============================================================================
# AIM 1
#
# DOES BASELINE ECOLOGICAL ARCHITECTURE PREDICT RESISTANCE?
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("AIM 1\n")
cat("BASELINE ARCHITECTURE -> ECOLOGICAL RESISTANCE\n")
cat("============================================================\n")


aim1_fit <- lm(

  ecological_resistance ~
    z_PC1 +
    z_PC2 +
    study,

  data =
    analysis

)


saveRDS(
  aim1_fit,
  file.path(
    OUT,
    "models",
    "Aim1_primary_model.rds"
  )
)


# ---------------------------------------------------------------------------
# Omnibus test: PC1 + PC2 jointly
# ---------------------------------------------------------------------------

aim1_omnibus <- cr2_joint(

  fit =
    aim1_fit,

  cluster =
    analysis$cow_global,

  terms =
    c(
      "z_PC1",
      "z_PC2"
    ),

  label =
    paste0(
      "Baseline ecological architecture predicts ",
      "quarter ecological resistance"
    )

)


cat("\nAim 1 omnibus CR2 test:\n")
print(aim1_omnibus)


write_csv(
  aim1_omnibus,
  file.path(
    OUT,
    "tables",
    "Aim1_omnibus_CR2.csv"
  )
)


# ---------------------------------------------------------------------------
# Individual robust coefficients
# ---------------------------------------------------------------------------

aim1_coefficients <- cr2_coefficients(

  fit =
    aim1_fit,

  cluster =
    analysis$cow_global

)


cat("\nAim 1 robust coefficients:\n")
print(aim1_coefficients)


write_csv(
  aim1_coefficients,
  file.path(
    OUT,
    "tables",
    "Aim1_CR2_coefficients.csv"
  )
)


# ============================================================================
# AIM 1 EFFECT-SIZE MODEL
#
# Report R2 descriptively — NOT as clustered inferential statistic.
# ============================================================================

aim1_model_summary <- summary(
  aim1_fit
)


aim1_effect_size <- tibble(

  n =
    nrow(
      analysis
    ),

  cows =
    n_distinct(
      analysis$cow_global
    ),

  R_squared =
    aim1_model_summary$r.squared,

  adjusted_R_squared =
    aim1_model_summary$adj.r.squared,

  residual_SD =
    sigma(
      aim1_fit
    )

)


cat("\nAim 1 model-level descriptive effect size:\n")
print(aim1_effect_size)


write_csv(
  aim1_effect_size,
  file.path(
    OUT,
    "tables",
    "Aim1_model_effect_size.csv"
  )
)


# ============================================================================
# AIM 2
#
# IS THE ASSOCIATION REPRODUCIBLE ACROSS ITALY AND MANITOBA?
#
# H0:
#   PC1 x study = 0
#   PC2 x study = 0
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("AIM 2\n")
cat("ITALY vs MANITOBA REPRODUCIBILITY\n")
cat("============================================================\n")


aim2_fit <- lm(

  ecological_resistance ~
    (
      z_PC1 +
      z_PC2
    ) *
    study,

  data =
    analysis

)


saveRDS(
  aim2_fit,
  file.path(
    OUT,
    "models",
    "Aim2_reproducibility_model.rds"
  )
)


coefficient_names <- names(
  coef(
    aim2_fit
  )
)


interaction_terms <- coefficient_names[
  grepl(
    "z_PC[12]:study|study.*:z_PC[12]",
    coefficient_names
  )
]


if (
  length(
    interaction_terms
  ) != 2
) {

  stop(
    "Expected exactly 2 PC x study interactions, found: ",
    paste(
      interaction_terms,
      collapse = ", "
    )
  )
}


aim2_omnibus <- cr2_joint(

  fit =
    aim2_fit,

  cluster =
    analysis$cow_global,

  terms =
    interaction_terms,

  label =
    paste0(
      "Baseline architecture-resistance relationship ",
      "differs between Italy and Manitoba"
    )

)


cat("\nAim 2 interaction CR2 test:\n")
print(aim2_omnibus)


write_csv(
  aim2_omnibus,
  file.path(
    OUT,
    "tables",
    "Aim2_reproducibility_CR2.csv"
  )
)


aim2_coefficients <- cr2_coefficients(

  fit =
    aim2_fit,

  cluster =
    analysis$cow_global

)


write_csv(
  aim2_coefficients,
  file.path(
    OUT,
    "tables",
    "Aim2_CR2_coefficients.csv"
  )
)


# ============================================================================
# STUDY-SPECIFIC EFFECTS
#
# These are descriptive estimates of the association within each cohort.
# Robust CR2 estimates are also calculated separately where estimable.
# ============================================================================

study_specific_results <- list()


for (s in PRIMARY_STUDIES) {

  d <- analysis %>%
    filter(
      study == s
    ) %>%
    droplevels()


  fit_s <- lm(

    ecological_resistance ~
      z_PC1 +
      z_PC2,

    data =
      d
  )


  saveRDS(
    fit_s,
    file.path(
      OUT,
      "models",
      paste0(
        "Aim2_",
        s,
        "_model.rds"
      )
    )
  )


  naive <- summary(
    fit_s
  )$coefficients


  robust <- tryCatch(

    cr2_coefficients(
      fit =
        fit_s,
      cluster =
        d$cow_global
    ),

    error =
      function(e) NULL
  )


  basic <- bind_rows(

    lapply(
      c(
        "z_PC1",
        "z_PC2"
      ),
      function(term) {

        tibble(

          study =
            s,

          term =
            term,

          estimate =
            naive[
              term,
              "Estimate"
            ],

          naive_SE =
            naive[
              term,
              "Std. Error"
            ],

          naive_p =
            naive[
              term,
              "Pr(>|t|)"
            ],

          trajectories =
            nrow(
              d
            ),

          cows =
            n_distinct(
              d$cow_global
            )

        )
      }
    )
  )


  if (!is.null(robust)) {

    robust_small <- robust %>%

      filter(
        term %in%
          c(
            "z_PC1",
            "z_PC2"
          )
      )


    basic <- basic %>%

      left_join(
        robust_small,
        by = "term",
        suffix = c(
          "",
          "_CR2"
        )
      )
  }


  study_specific_results[[s]] <-
    basic
}


study_specific_results <- bind_rows(
  study_specific_results
)


cat("\nStudy-specific slopes:\n")
print(study_specific_results)


write_csv(
  study_specific_results,
  file.path(
    OUT,
    "tables",
    "Aim2_study_specific_slopes.csv"
  )
)


# ============================================================================
# TAXONOMIC INTERPRETATION OF ECOLOGICAL AXES
# ============================================================================

TOP_N <- 25


top_PC1 <- loadings %>%

  arrange(
    desc(
      abs(
        PC1_loading
      )
    )
  ) %>%

  slice_head(
    n = TOP_N
  )


top_PC2 <- loadings %>%

  arrange(
    desc(
      abs(
        PC2_loading
      )
    )
  ) %>%

  slice_head(
    n = TOP_N
  )


write_csv(
  top_PC1,
  file.path(
    OUT,
    "tables",
    "top25_PC1_family_loadings.csv"
  )
)


write_csv(
  top_PC2,
  file.path(
    OUT,
    "tables",
    "top25_PC2_family_loadings.csv"
  )
)


# ============================================================================
# LINK PC LOADINGS TO DIRECTION OF RESISTANCE ASSOCIATION
#
# This gives an interpretable table:
#
# If beta_PC1 > 0:
#   positive PC1 loading -> associated direction toward greater resistance
#
# If beta_PC1 < 0:
#   positive PC1 loading -> associated direction toward greater susceptibility
#
# Same logic for PC2.
#
# This is DESCRIPTIVE INTERPRETATION, not taxon-level causal inference.
# ============================================================================

beta_PC1 <- coef(
  aim1_fit
)[
  "z_PC1"
]


beta_PC2 <- coef(
  aim1_fit
)[
  "z_PC2"
]


axis_interpretation <- loadings %>%

  mutate(

    PC1_resistance_direction_score =
      PC1_loading *
      beta_PC1,

    PC2_resistance_direction_score =
      PC2_loading *
      beta_PC2,

    combined_direction_score =
      PC1_resistance_direction_score +
      PC2_resistance_direction_score

  ) %>%

  arrange(
    desc(
      abs(
        combined_direction_score
      )
    )
  )


write_csv(
  axis_interpretation,
  file.path(
    OUT,
    "tables",
    "family_axis_resistance_interpretation.csv"
  )
)


write_csv(
  axis_interpretation %>%
    slice_head(
      n = 30
    ),
  file.path(
    OUT,
    "tables",
    "top30_family_axis_resistance_interpretation.csv"
  )
)


# ============================================================================
# FIGURE 1
# ECOLOGICAL RESISTANCE DISTRIBUTION
# ============================================================================

p1 <- ggplot(

  analysis,

  aes(
    x = study,
    y = ecological_resistance
  )

) +

  geom_boxplot(
    outlier.shape = NA
  ) +

  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.75,
    size = 2
  ) +

  labs(

    x =
      NULL,

    y =
      "Ecological resistance (1 - Bray-Curtis)",

    title =
      "Healthy mammary-quarter ecological resistance"

  ) +

  theme_classic(
    base_size = 12
  )


ggsave(

  file.path(
    OUT,
    "figures",
    "Figure1_ecological_resistance.pdf"
  ),

  p1,

  width = 6,
  height = 4.5

)


# ============================================================================
# FIGURE 2
# BASELINE ECOLOGICAL ARCHITECTURE
# ============================================================================

p2 <- ggplot(

  analysis,

  aes(
    x = PC1,
    y = PC2,
    shape = study
  )

) +

  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +

  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +

  geom_point(
    size = 2.6,
    alpha = 0.8
  ) +

  labs(

    x =
      paste0(
        "PC1 (",
        round(
          variance[1] *
            100,
          1
        ),
        "%)"
      ),

    y =
      paste0(
        "PC2 (",
        round(
          variance[2] *
            100,
          1
        ),
        "%)"
      ),

    title =
      "Baseline ecological architecture of healthy mammary quarters"

  ) +

  theme_classic(
    base_size = 12
  )


ggsave(

  file.path(
    OUT,
    "figures",
    "Figure2_baseline_ecological_architecture.pdf"
  ),

  p2,

  width = 6.5,
  height = 5

)


# ============================================================================
# FIGURE 3
# PC1 vs RESISTANCE
# ============================================================================

p3 <- ggplot(

  analysis,

  aes(
    x = z_PC1,
    y = ecological_resistance,
    shape = study
  )

) +

  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +

  geom_smooth(
    method = "lm",
    se = TRUE
  ) +

  labs(

    x =
      "Baseline ecological architecture PC1 (SD)",

    y =
      "Ecological resistance (1 - Bray-Curtis)",

    title =
      "Baseline PC1 and subsequent ecological resistance"

  ) +

  theme_classic(
    base_size = 12
  )


ggsave(

  file.path(
    OUT,
    "figures",
    "Figure3_PC1_resistance.pdf"
  ),

  p3,

  width = 6,
  height = 4.5

)


# ============================================================================
# FIGURE 4
# PC2 vs RESISTANCE
# ============================================================================

p4 <- ggplot(

  analysis,

  aes(
    x = z_PC2,
    y = ecological_resistance,
    shape = study
  )

) +

  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +

  geom_smooth(
    method = "lm",
    se = TRUE
  ) +

  labs(

    x =
      "Baseline ecological architecture PC2 (SD)",

    y =
      "Ecological resistance (1 - Bray-Curtis)",

    title =
      "Baseline PC2 and subsequent ecological resistance"

  ) +

  theme_classic(
    base_size = 12
  )


ggsave(

  file.path(
    OUT,
    "figures",
    "Figure4_PC2_resistance.pdf"
  ),

  p4,

  width = 6,
  height = 4.5

)


# ============================================================================
# MASTER HYPOTHESIS TABLE
# ============================================================================

master <- bind_rows(

  aim1_omnibus %>%

    transmute(

      aim =
        "Aim 1",

      biological_question =
        paste0(
          "Does baseline ecological architecture predict ",
          "subsequent ecological resistance?"
        ),

      hypothesis =
        hypothesis,

      Fstat =
        Fstat,

      df_num =
        df_num,

      df_denom =
        df_denom,

      p_value =
        p_value,

      interpretation =
        ifelse(
          p_value < 0.05,
          "Supported",
          "Not supported"
        )

    ),


  aim2_omnibus %>%

    transmute(

      aim =
        "Aim 2",

      biological_question =
        paste0(
          "Does the architecture-resistance relationship ",
          "differ between Italy and Manitoba?"
        ),

      hypothesis =
        hypothesis,

      Fstat =
        Fstat,

      df_num =
        df_num,

      df_denom =
        df_denom,

      p_value =
        p_value,

      interpretation =
        ifelse(
          p_value < 0.05,
          "Evidence of cohort-specific relationship",
          "Insufficient evidence of cohort heterogeneity"
        )

    )

)


write_csv(
  master,
  file.path(
    OUT,
    "tables",
    "PRIMARY_TWO_AIMS_RESULTS.csv"
  )
)


# ============================================================================
# HUMAN-READABLE SUMMARY
# ============================================================================

summary_file <- file.path(
  OUT,
  "PRIMARY_TWO_AIMS_SUMMARY.txt"
)


sink(
  summary_file
)


cat(
  "MMER PRIMARY HEALTHY TWO-AIM ANALYSIS\n"
)

cat(
  "======================================\n\n"
)


cat(
  "PRIMARY STUDIES\n"
)

cat(
  "Italy and Manitoba only.\n"
)

cat(
  "Wisconsin excluded from this analysis.\n\n"
)


cat(
  "ANALYSIS POPULATION\n"
)

print(
  sample_summary
)


cat(
  "\n\nECOLOGICAL RESISTANCE\n"
)

print(
  phenotype_summary
)


cat(
  "\n\nPCA VARIANCE\n"
)

print(
  variance_table %>%
    slice_head(
      n = 5
    )
)


cat(
  "\n\nAIM 1 — BASELINE ARCHITECTURE -> RESISTANCE\n"
)

print(
  aim1_omnibus
)


cat(
  "\nAim 1 robust coefficients:\n"
)

print(
  aim1_coefficients
)


cat(
  "\nAim 1 effect size:\n"
)

print(
  aim1_effect_size
)


cat(
  "\n\nAIM 2 — ITALY vs MANITOBA REPRODUCIBILITY\n"
)

print(
  aim2_omnibus
)


cat(
  "\nStudy-specific slopes:\n"
)

print(
  study_specific_results
)


cat(
  "\n\nTOP PC1 LOADINGS\n"
)

print(
  top_PC1,
  n = 15
)


cat(
  "\n\nTOP PC2 LOADINGS\n"
)

print(
  top_PC2,
  n = 15
)


cat(
  "\n\nTOP FAMILY-AXIS RESISTANCE INTERPRETATION\n"
)

print(
  axis_interpretation %>%
    slice_head(
      n = 20
    ),
  n = 20
)


sink()


# ============================================================================
# SESSION INFO
# ============================================================================

capture.output(

  sessionInfo(),

  file =
    file.path(
      OUT,
      "sessionInfo.txt"
    )
)


# ============================================================================
# FINAL CONSOLE REPORT
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PRIMARY TWO-AIM RESULTS\n")
cat("============================================================\n")

print(
  master
)


cat("\n")
cat("Aim 1 CR2 coefficients:\n")
print(
  aim1_coefficients
)


cat("\n")
cat("Study-specific slopes:\n")
print(
  study_specific_results
)


cat("\n")
cat("Top family-axis resistance interpretation:\n")
print(
  axis_interpretation %>%
    slice_head(
      n = 15
    ),
  n = 15
)


cat("\n")
cat("============================================================\n")
cat("[PASS] PRIMARY HEALTHY TWO-AIM ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "Output directory:",
  OUT,
  "\n"
)
