# MMER

**Mammary Metataxonomy Ecological Resistance**

MMER is a reproducible multi-cohort reanalysis framework for testing whether the **baseline ecological organization of the bovine mammary microbiome is associated with its subsequent ecological resistance**.

Ecological resistance is defined as:

`ecological_resistance = 1 - Bray-Curtis(baseline, follow-up)`

Higher values indicate greater compositional retention relative to the unit's own baseline.

## Current scientific question

The current MMER hypothesis is:

> **Does baseline multivariate mammary-community architecture contain information about how strongly that community resists subsequent ecological displacement, and does that phenomenon generalize across independent cohorts?**

The project now distinguishes clearly between **discovery/development cohorts** and **external generalization cohorts**.

## Current four-study evidence structure

| Role | Cohort | Platform | Primary unit / analysis | Current result |
|---|---|---|---|---|
| Discovery | Italy / Biscarini et al. / PRJEB38332 | 16S | healthy quarter trajectories | contributes to significant pooled baseline-architecture signal |
| Discovery | Manitoba / Derakhshani et al. | 16S | healthy quarter trajectories | contributes to significant pooled baseline-architecture signal |
| External generalization | Porcellato | 16S | 154 same-quarter Sampling1→Sampling2 trajectories | **strong positive external validation** |
| External generalization | Patangia | shotgun metagenomics / Kraken2 family profiles | 24 cows with M0, M2, M4, M6 | **primary taxonomic test negative; later-window signal exploratory only** |

Van Beeck and Wisconsin remain useful contextual/sensitivity cohorts but are **not part of the frozen four-study primary architecture evidence summarized here**.

---

## Discovery analysis — Italy + Manitoba

The frozen primary healthy analysis contains **49 quarter trajectories from 14 cows**:

- Italy: 20 trajectories / 5 cows
- Manitoba: 29 trajectories / 9 cows

Baseline family-level composition is represented using CLR-transformed community profiles with study-level centering before pooled PCA.

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Cow-clustered CR2 inference:

- joint PC1 + PC2 HTZ F = **10.0**
- df = **2, 7.81**
- p = **0.00699**
- PC1 beta = **-0.0286**, p = 0.124
- PC2 beta = **+0.0622**, p = **0.00302**
- descriptive R² = **0.334**
- adjusted R² = **0.289**

Interpretation:

> Baseline multivariate community architecture is associated with subsequent ecological resistance across the pooled Italy + Manitoba healthy-quarter analysis.

The study-interaction test was not conventionally significant (`p = 0.0766`), supporting a shared multivariate phenomenon while allowing cohort heterogeneity.

---

## External validation — Porcellato 16S cohort

Porcellato was independently reprocessed as a same-quarter longitudinal 16S cohort.

Trajectory-depth QC identified **189 complete Sampling1→Sampling2 trajectories**. The locked pre-model criterion of at least 1,000 DADA2 non-chimeric reads at both timepoints retained **154 trajectories**:

- Farm A: 93
- Farm K: 61

The DADA2/taxonomy pipeline produced:

- 378 samples
- 15,567 ASVs
- 213 bacterial families overall
- 308 samples × 211 families in the locked 154-trajectory resistance matrix
- 50 baseline families retained for PCA

Baseline PCA variance explained:

- PC1: **19.779%**
- PC2: **8.823%**

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + farm`

with cow-clustered CR2 inference because cows can contribute multiple quarter trajectories.

Primary result (`n = 154`):

- PC1 beta = **-0.0892**, p < 0.001
- PC2 beta = **+0.0636**, p = **0.0106**
- farm K beta = **+0.3061**, p < 0.001
- joint PC1 + PC2 HTZ F = **9.49**
- df = **2, 23.5**
- p < **0.001**
- descriptive R² = **0.410**

### Assigned-family-read robustness

The result remained stable under stricter post-taxonomy depth requirements:

- ≥1,000 assigned family reads at both timepoints: `n = 143`, joint F = **11.9**, p < 0.001
- ≥2,000 assigned family reads at both timepoints: `n = 121`, joint F = **7.87**, p = **0.00278**

Porcellato therefore provides strong **independent phenomenon-level replication** of baseline ecological architecture predicting subsequent same-quarter ecological resistance.

Because Porcellato's PCA is learned independently, this does **not** establish transport of identical discovery PC axes across cohorts.

---

## External generalization — Patangia shotgun-metagenomic cohort

Patangia contains **24 cows with complete M0, M2, M4, and M6 samples**. Family-level taxonomic profiles were derived from shotgun reads with Kraken2.

The validated family matrix contains:

- 96 samples
- 852 families overall
- finite relative-abundance profiles with row sums equal to 1
- 680 M0 families retained after the ≥10% baseline-prevalence filter

M0 baseline PCA:

- PC1 variance explained: **8.958%**
- PC2 variance explained: **6.580%**

### Prespecified M0→M2 result

Mean ecological resistance = **0.441**.

- PC1 beta = **+0.0334**
- PC2 beta = **-0.0252**
- R² = **0.0275**
- no significant joint baseline PC1 + PC2 association

**Conclusion:** the primary Patangia taxonomic architecture test did not validate the discovery phenomenon at M0→M2.

### Additional intervals

M2→M4 and M4→M6 were also negative under standard and HC3 inference.

For exploratory M2→M6:

- mean ecological resistance = **0.503**
- PC1 beta = +0.0430
- PC2 beta = +0.0723
- R² = 0.1448
- standard joint F-test p = 0.194
- HC3 PC2 p = 0.0618
- HC3 asymptotic joint Wald p = 0.0231

Because the conventional joint test is null and only 24 independent cows are available, M2→M6 is treated as **suggestive/exploratory**, not confirmatory.

The current evidence does **not** justify claiming a 16S-versus-shotgun platform effect. Platform, sample size, biological window, host/background DNA, and cohort context remain confounded.

---

## Current synthesis

The evidence currently supports the following statement:

> **Baseline mammary microbial community architecture is associated with subsequent ecological resistance in the pooled Italy + Manitoba discovery analysis and independently generalizes in the Porcellato 16S cohort, whereas taxonomic baseline architecture does not robustly predict resistance in the smaller Patangia shotgun-metagenomic cohort.**

The data do not establish:

- a universal taxonomic PC axis;
- universal effects of individual families;
- a causal mechanism;
- a sequencing-platform-specific effect;
- guaranteed transportability across all biological windows.

This pattern is scientifically informative because MMER tests **generalization**, not whether every cohort can be made statistically positive.

---

## Core aims

1. **Quantify mammary ecological resistance** at the longitudinal biological-unit level.
2. **Determine whether baseline ecological architecture predicts subsequent resistance**, while separating the phenomenon from treatment/time-average effects.
3. **Characterize the dimensions and generalizability of ecological resistance** across cohorts, biological contexts, and sequencing modalities.

---

## Key documents

- [`docs/four_study_ecological_resistance_freeze_2026-08-16.md`](docs/four_study_ecological_resistance_freeze_2026-08-16.md) — current four-study evidence freeze
- [`docs/external_generalization_status_2026-08-16.md`](docs/external_generalization_status_2026-08-16.md) — detailed Porcellato and Patangia external analyses
- [`docs/frozen_primary_healthy_two_aims_2026-08-16.md`](docs/frozen_primary_healthy_two_aims_2026-08-16.md) — frozen Italy + Manitoba primary healthy analysis
- [`docs/MMER_source_papers_critical_comparison.md`](docs/MMER_source_papers_critical_comparison.md) — conceptual distinction from the Italy and Manitoba source studies

---

## Reproducibility guardrails

- Preserve the Italy + Manitoba healthy-quarter analysis as the discovery/development framework.
- Preserve Porcellato `n = 154` as the primary external 16S analysis; `n = 143` and `n = 121` are depth robustness checks.
- Preserve Patangia M0→M2 as the primary shotgun taxonomic test.
- Keep Patangia M2→M6 explicitly exploratory.
- Use cow-clustered inference when cows contribute multiple quarter trajectories.
- Do not use cow-clustered CR2 when an endpoint-specific analysis has one independent observation per cow; use ordinary or heteroskedasticity-robust inference as appropriate.
- Do not infer a 16S-versus-shotgun effect from one shotgun cohort.
- Do not reinterpret taxonomic families defining PCA axes as individually causal biomarkers.
- Exact cross-cohort axis transport requires harmonized features and projection of frozen discovery loadings; independent PCA supports phenomenon-level replication only.
- Do not promote secondary time windows or sensitivity analyses to primary because they produce smaller p-values.

## Status

**Italy + Manitoba discovery analysis:** complete and frozen.

**Porcellato external 16S validation:** complete; strong positive replication.

**Patangia shotgun taxonomic external test:** complete; primary null with limited exploratory later-window evidence.

**Current priority:** consolidate the four-study result into publication-ready Aim 1–3 results, figures, and interpretation without changing the frozen primary hypotheses.
