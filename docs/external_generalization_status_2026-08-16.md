# MMER external generalization status — 2026-08-16

This note freezes the current cross-cohort evidence for the MMER hypothesis that **baseline mammary microbial community architecture is associated with subsequent ecological resistance**, where ecological resistance is defined as:

`ecological_resistance = 1 - Bray-Curtis(baseline, follow-up)`

Higher values indicate greater compositional retention/stability.

## 1. Primary healthy analysis: Italy + Manitoba

The frozen primary healthy analysis uses 49 quarter trajectories from 14 cows (Italy: 20 trajectories/5 cows; Manitoba: 29 trajectories/9 cows).

Baseline family-level composition was transformed using CLR, centered within study, and summarized by pooled PCA. The primary model was:

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Cow-clustered CR2 inference gave a significant joint PC1+PC2 architecture test:

- joint HTZ F = 10.0
- df = 2, 7.81
- p = 0.00699
- PC1 beta = -0.0286, p = 0.124
- PC2 beta = +0.0622, p = 0.00302
- descriptive R2 = 0.334
- adjusted R2 = 0.289

The study interaction test was not conventionally significant (joint interaction p = 0.0766), so the evidence supports a shared multivariate architecture signal while allowing cohort heterogeneity.

## 2. Porcellato external generalization

Porcellato was reprocessed independently as a same-quarter longitudinal 16S cohort.

There were 189 complete Sampling1 -> Sampling2 trajectories before sequencing-depth filtering. The pre-model DADA2 non-chimeric-read criterion retained 154 trajectories for the external analysis.

Within-farm CLR centering was performed before PCA. The model was:

`ecological_resistance ~ z_PC1 + z_PC2 + farm`

with cow-clustered CR2 inference because cows contributed multiple quarter trajectories.

### Primary external result (n = 154)

- PC1 beta = -0.0892, CR2 SE = 0.0227, p < 0.001
- PC2 beta = +0.0636, CR2 SE = 0.0228, p = 0.0106
- farm K beta = +0.3061, p < 0.001
- joint PC1+PC2 HTZ F = 9.49
- df = 2, 23.5
- p < 0.001
- descriptive R2 = 0.410

### Assigned-family-read robustness

Post-taxonomy assigned bacterial family-read thresholds gave nearly unchanged conclusions:

- >=1000 assigned family reads at both timepoints: n = 143; joint F = 11.9, p < 0.001
- >=2000 assigned family reads at both timepoints: n = 121; joint F = 7.87, p = 0.00278

These are robustness analyses. The n = 154 DADA2-filtered analysis remains the external primary analysis because it was the criterion used before viewing the model result.

### Interpretation

Porcellato independently supports the **general phenomenon** that baseline multivariate community architecture is associated with subsequent same-quarter ecological resistance.

Its PCA was learned independently. Therefore this is a phenomenon-level replication, not yet proof that the exact Italy/Manitoba PC axes transport unchanged across cohorts.

## 3. Patangia shotgun-metagenomic cohort

Patangia contains 24 cows with complete M0, M2, M4, and M6 samples. Family-level taxonomic profiles were generated from shotgun reads using Kraken2.

Because there is one trajectory per cow for each interval, ordinary independent-observation regression with HC3 robust covariance is used rather than cow-clustered CR2.

### M0 -> M2

- n = 24
- mean ecological resistance = 0.441
- PC1 beta = +0.0334
- PC2 beta = -0.0252
- R2 = 0.0275
- standard joint F-test p = 0.746
- HC3 joint Wald p = 0.784

Conclusion: no evidence that M0 baseline taxonomic architecture predicts M2 resistance.

### M2 -> M4

- n = 24
- mean ecological resistance = 0.576
- PC1 beta = +0.0555
- PC2 beta = +0.0102
- R2 = 0.0506
- standard joint F-test p = 0.580
- HC3 joint Wald p = 0.636

Conclusion: no evidence that M2 architecture predicts M4 resistance over this interval.

### M4 -> M6

- n = 24
- mean ecological resistance = 0.530
- PC1 beta = -0.0696
- PC2 beta = -0.0111
- R2 = 0.0819
- standard joint F-test p = 0.408
- HC3 joint Wald p = 0.392

Conclusion: no evidence that M4 architecture predicts M6 resistance.

### Exploratory M2 -> M6

- n = 24
- mean ecological resistance = 0.503
- PC1 beta = +0.0430
- PC2 beta = +0.0723
- R2 = 0.1448
- standard joint F-test p = 0.194
- HC3 PC2 p = 0.0618
- HC3 asymptotic joint Wald p = 0.0231

The ordinary joint test is not significant, whereas the asymptotic HC3 Wald test is. With only 24 independent cows, this discrepancy is treated as **suggestive/exploratory**, not confirmatory evidence.

### Patangia interpretation

Across the biologically sensible short intervals, Patangia does not show a robust taxonomic architecture -> ecological resistance association. The M2 -> M6 exploratory result is insufficient to overturn that conclusion.

This cohort differs from the positive 16S cohorts in both sequencing platform and biological context. Therefore the current evidence does **not** justify claiming that the phenomenon is 16S-specific. Platform, sample size, postpartum/lactation window, host/background DNA, and cohort biology remain confounded.

A future shotgun-specific extension can test whether baseline **functional metagenomic architecture** (e.g., pathway/gene-family composition) predicts subsequent resistance, but this should be predefined as a distinct hypothesis rather than used to rescue a null taxonomic result.

## 4. Current cross-cohort synthesis

Current evidence supports the following cautious statement:

> Baseline mammary community architecture is associated with subsequent ecological resistance in the primary Italy+Manitoba analysis and independently generalizes to Porcellato, while taxonomic architecture does not robustly predict resistance in the smaller Patangia shotgun-metagenomic cohort.

The current evidence does not establish a universal taxonomic axis, a universal direction for individual taxa, a causal mechanism, or a sequencing-platform-specific effect.

## 5. Next cohort

Van Beeck / PRJEB63336 is already processed locally with 60 complete 3k trajectories and family/genus-level outputs. It will be revisited next using the current baseline-architecture -> ecological-resistance hypothesis, without promoting prior secondary analyses to primary status after seeing results.

## Reproducibility / inference guardrails

- Do not claim causality from these retrospective observational associations.
- Keep Porcellato n=154 as the external primary analysis; n=143 and n=121 are post-taxonomy robustness checks.
- Keep Patangia M2 -> M6 explicitly exploratory.
- Do not infer a 16S-vs-shotgun platform effect from one shotgun cohort.
- Do not search additional time windows or taxa solely because a predefined test was null.
- Independent cohort PCA tests phenomenon-level generalization; exact axis transportability requires projecting frozen discovery loadings into a harmonized external feature space.
