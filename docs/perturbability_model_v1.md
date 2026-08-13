# MMER Study 1 Perturbability Model v1

## Purpose

This model converts the Q3 discovery result into a portable, falsifiable prediction rule for **future mammary microbiome displacement**. It is intentionally not a mastitis-risk model and is not intended for clinical decision-making.

## Discovery cohort

The model is trained on the Biscarini et al. PRJEB38332 production analysis. The modeling unit is the **cow-quarter trajectory**. T2 and T3 Bray-Curtis displacement values are averaged, yielding 20 independent trajectory-level targets from 5 cows.

### Target

`mean_bray_displacement = mean[ d(T1,T2), d(T1,T3) ]`

where each distance is calculated from the same quarter's own T1 microbiome.

### Baseline predictors

Only four transportable ecological features are used in v1:

- baseline richness
- baseline Shannon diversity
- baseline evenness
- baseline dominance

Taxonomic predictors and ASV-level PC axes are excluded from v1 because the planned Manitoba validation cohort uses a different 16S region, making direct ASV-space transport inappropriate.

## Why ridge regression?

The discovery cohort is small (20 trajectories) and the four ecological predictors are correlated. A flexible machine-learning model would overfit. Ridge regression shrinks correlated coefficients while retaining all four prespecified ecological variables.

No treatment variable is included in the portable model. The objective is not to recreate the original treatment-effect model but to estimate an intrinsic baseline **perturbability phenotype** that can be evaluated under a different perturbation regime.

## Internal validation

Performance is estimated using **nested leave-one-cow-out cross-validation**:

1. One entire cow (all four quarters) is held out as the outer test fold.
2. Ridge penalty strength is selected only within the remaining cows by inner leave-one-cow-out CV.
3. The locked inner model predicts all four quarters of the held-out cow.
4. This is repeated for all five cows.

This prevents quarter-level leakage between train and test sets.

Reported discovery performance includes RMSE, MAE, Pearson correlation, and Spearman rank correlation between predicted and observed trajectory displacement.

## Bootstrap uncertainty

Coefficient uncertainty is estimated by **cow-level bootstrap resampling**. A cow is resampled as a block, carrying all four quarter trajectories together. The default is 5,000 bootstrap replicates (`MMER_BOOTSTRAPS=5000`).

Bootstrap coefficient draws are exported so uncertainty can be propagated into external predictions without using external follow-up outcomes.

## Frozen external model

After internal performance is estimated by nested CV, the final model is fit to all 20 discovery trajectories. Ridge lambda is selected using cow-blocked leave-one-cow-out CV across the discovery cohort.

The frozen model exports:

- `MMER_Study1_Perturbability_Ridge_v1.rds`
- `frozen_model_coefficients.csv`
- `frozen_model_scaling.csv`
- `cow_bootstrap_raw_coefficients.csv`
- `model_manifest.csv`

These files should be archived before inspecting validation outcomes.

## Manitoba validation plan

The Manitoba cohort is intended as an **external validation of ecological perturbability**, not as a treatment-effect replication.

For each paired pre-dry-off milk sample:

1. Run the same taxonomic/QC pipeline as far as technically compatible.
2. Calculate the four baseline ecological predictors.
3. Apply the frozen MMER model using `scripts/07b_apply_perturbability_model.R`.
4. Save and freeze predicted perturbability before calculating/unblinding post-calving displacement.
5. Calculate observed pre-DCT-to-post-calving Bray-Curtis displacement.
6. Compare predicted and observed perturbability using Spearman rank correlation as the primary transportability endpoint, with Pearson correlation, RMSE, MAE, and calibration as secondary metrics.

Because sequencing region and laboratory protocol differ across cohorts, **rank transportability is the primary validation target**. Exact absolute calibration may require study-specific recalibration and must be reported separately from the locked validation.

## Interpretation

A successful validation would support the statement:

> Baseline mammary ecological architecture contains transportable information about the magnitude of subsequent community displacement in an independent cohort.

It would not by itself show that perturbability causes disease, that high displacement is harmful, or that the model predicts mastitis.

## Planned next-stage disease test

Only after perturbability is externally validated should the frozen score be tested in longitudinal mastitis cohorts. The proposed next question is whether predicted ecological perturbability at apparently healthy preclinical timepoints (for example, 6 or 4 weeks before S. aureus clinical mastitis) is associated with subsequent infection burden or disease emergence.

That disease analysis is a separate hypothesis and should not be used to modify the Study 1 model before the external perturbability validation is completed.
