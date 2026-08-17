# MMER five-study ecological-resistance evidence freeze — 2026-08-16

## Status

**Current evidence freeze.** This document supersedes the earlier four-study freeze created before completion of the Van Beeck external-validation analysis.

The MMER hypothesis is:

> **Baseline multivariate mammary microbial community architecture is associated with the subsequent ecological resistance of the same longitudinal biological unit.**

Ecological resistance is defined as:

`R = 1 - Bray-Curtis(baseline, follow-up)`

Higher values represent greater compositional retention relative to the unit's own baseline.

The evidence structure now contains **five source studies** with distinct scientific roles:

1. Italy — discovery/development;
2. Manitoba — discovery/development;
3. Porcellato — independent 16S external validation;
4. Van Beeck — independent 16S external validation;
5. Patangia — cross-platform shotgun taxonomic external generalization.

---

## 1. Study-role table

| Study | Role | Platform | Primary longitudinal unit | Frozen outcome/model | Result |
|---|---|---|---|---|---|
| Italy / Biscarini / PRJEB38332 | discovery/development | 16S | healthy mammary-quarter trajectory | pooled Italy+Manitoba architecture model | positive pooled discovery signal |
| Manitoba / Derakhshani | discovery/development | 16S | healthy mammary-quarter trajectory | pooled Italy+Manitoba architecture model | positive pooled discovery signal |
| Porcellato / PRJEB35792 | independent external validation | 16S | same-quarter Sampling1→Sampling2 trajectory | `R ~ PC1 + PC2 + farm` | **positive** |
| Van Beeck / PRJEB63336 | independent external validation | 16S | complete longitudinal milk trajectory | `overall_R ~ PC1 + PC2 + treatment + dairy` | **positive** |
| Patangia | external cross-platform generalization | shotgun / Kraken2 family profiles | cow trajectory | prespecified M0→M2 `R ~ PC1 + PC2` | **primary null** |

This structure deliberately preserves both successful validation and principled non-replication.

---

# 2. Frozen discovery: Italy + Manitoba

The primary healthy discovery analysis contains **49 mammary-quarter trajectories from 14 cows**:

- Italy: 20 trajectories / 5 cows;
- Manitoba: 29 trajectories / 9 cows.

Baseline family-level compositions are CLR transformed and centered within study before pooled PCA.

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Cow-clustered CR2 inference is used because cows contribute multiple quarter trajectories.

Frozen result:

- joint PC1 + PC2 HTZ F = **10.0**;
- df = **2, 7.81**;
- p = **0.00699**;
- PC1 beta = **-0.0286**, p = 0.124;
- PC2 beta = **+0.0622**, p = **0.00302**;
- descriptive R² = **0.334**;
- adjusted R² = **0.289**.

The study-interaction test was not conventionally significant (`p = 0.0766`).

### Discovery conclusion

> Baseline multivariate community architecture is associated with subsequent ecological resistance across the pooled Italy + Manitoba healthy-quarter analysis.

This is the phenomenon tested externally. It does not imply a universal fixed PC1 or PC2 axis.

---

# 3. Porcellato: independent 16S external validation

Porcellato provides a same-quarter independent longitudinal validation cohort.

Locked primary analysis:

- N = **154** Sampling1→Sampling2 quarter trajectories;
- Farm A = 93;
- Farm K = 61;
- family-level community representation;
- 50 baseline families retained for PCA;
- independently fitted Porcellato PCA.

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + farm`

Cow-clustered CR2 is appropriate because some cows contribute multiple quarter trajectories.

Frozen result:

- PC1 beta = **-0.0892**, p < 0.001;
- PC2 beta = **+0.0636**, p = **0.0106**;
- farm K beta = **+0.3061**, p < 0.001;
- joint PC1 + PC2 HTZ F = **9.49**;
- df = **2, 23.5**;
- p < **0.001**;
- descriptive R² = **0.410**.

Post-taxonomy assigned-read sensitivity remains positive at stricter thresholds:

- >=1,000 assigned family reads: N = 143, joint F = 11.9, p < 0.001;
- >=2,000 assigned family reads: N = 121, joint F = 7.87, p = 0.00278.

### Porcellato conclusion

> Porcellato provides strong independent phenomenon-level validation that baseline family-level community architecture is associated with subsequent same-quarter ecological resistance.

---

# 4. Van Beeck: independent 16S external validation

## 4.1 Source structure

Van Beeck et al. sampled milk and teat-skin microbiota at three California dairies at:

- Baseline / dry-off;
- 7 Days;
- 55–75 DIM in the subsequent lactation.

The MMER external-validation workflow uses only the locked longitudinal milk trajectories.

Locked external cohort:

- **60 trajectories**;
- **180 samples**;
- control low SCC = 14;
- control high SCC = 12;
- CB = 18;
- CH = 16.

## 4.2 Primary phenotype

For trajectory *i*:

`R_early,i = 1 - Bray(Baseline_i, 7Days_i)`

`R_late,i = 1 - Bray(Baseline_i, 55-75DIM_i)`

Frozen overall phenotype:

`R_overall,i = mean(R_early,i, R_late,i)`

The use of both post-baseline observations avoids selecting the more favorable time window after seeing the outcome.

## 4.3 Baseline architecture

- bacterial family relative abundance;
- baseline prevalence >=10%;
- zeros replaced with half the smallest positive baseline relative abundance;
- re-closure;
- CLR transformation;
- feature centering within dairy;
- independent Van Beeck PCA;
- standardized PC1 + PC2.

Treatment is **not** centered out of the PCA. It is retained as a design covariate in the regression.

## 4.4 Primary model

Reduced:

`overall_resistance ~ treatment + dairy`

Full:

`overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy`

Primary null:

`H0: beta_PC1 = beta_PC2 = 0`

There is one outcome observation per trajectory; singleton-cluster CR2 is therefore not used.

## 4.5 Frozen result

- N = **60**;
- full R² = **0.316**;
- adjusted R² = **0.224**;
- partial R² PC1+PC2 = **0.216**;
- classical F(2,52) = **7.18**, p = **0.00177**;
- HC3 joint Wald χ²(2) = **9.13**, p = **0.0104**;
- 9,999-stratum Freedman–Lane permutation p = **0.0183**.

Individual axes:

- PC1 beta = **-0.0172**, HC3 p = 0.283;
- PC2 beta = **-0.0665**, HC3 p = **0.0246**.

Bootstrap:

- 5,000 / 5,000 valid stratified complete-trajectory replicates;
- PC1 95% interval = **-0.0439 to +0.00665**;
- PC2 95% interval = **-0.104 to -0.0217**;
- 99.6% of PC2 bootstrap estimates negative;
- partial-R² bootstrap median = **0.223**;
- partial-R² percentile interval = **0.0499 to 0.566**.

### Van Beeck conclusion

> Baseline family-level ecological architecture significantly predicts subsequent ecological resistance in an independent multifactorial 16S cohort after adjustment for treatment and dairy.

This constitutes **positive external validation**.

The replicated object is the ecological-architecture phenomenon. Van Beeck PC directions are not treated as identical to discovery PCs because the PCA is independently fitted.

Canonical documentation:

`docs/vanbeeck_external_validation_workflow_2026-08-16.md`

Canonical analysis:

`scripts/26j_vanbeeck_external_validation_only.R`

Canonical figures:

`scripts/26k_vanbeeck_external_validation_plots.R`

---

# 5. Patangia: shotgun taxonomic external generalization

Patangia contains **24 cows** with complete M0, M2, M4, and M6 observations and family-level shotgun taxonomic profiles derived from Kraken2.

The prespecified primary taxonomic external test is M0→M2.

Primary result:

- mean resistance = **0.441**;
- PC1 beta = **+0.0334**;
- PC2 beta = **-0.0252**;
- R² = **0.0275**;
- standard joint p ≈ **0.746**;
- HC3 joint p ≈ **0.784**.

### Patangia conclusion

> The prespecified shotgun taxonomic architecture test does not validate the discovery phenomenon at M0→M2.

Later-window M2→M6 evidence remains exploratory because conventional and HC3 joint tests disagree and the cohort contains only 24 independent cows.

Patangia is therefore retained as a **primary null external generalization result**, not converted into a positive cohort by time-window searching.

---

# 6. Five-study synthesis

The strongest current synthesis is:

> **Baseline mammary microbial community architecture predicts subsequent ecological resistance in the pooled Italy + Manitoba discovery analysis and independently replicates in both Porcellato and Van Beeck 16S cohorts. The prespecified Patangia shotgun taxonomic test is null, indicating that the taxonomic architecture-resistance phenomenon is not universally detectable across all cohorts, windows, or feature-generation settings.**

This evidence structure is stronger than requiring every cohort to be positive because it distinguishes reproducible signal from context dependence.

A useful conceptual model is:

`Ecological resistance = f(baseline community architecture, host/background context, biological window, environment, perturbation context)`

with baseline architecture acting as a reproducible but not universally sufficient component.

---

# 7. What the current evidence supports

1. Ecological resistance varies substantially among longitudinal mammary biological units.
2. Baseline multivariate family-level architecture contains information about subsequent resistance in the Italy + Manitoba discovery analysis.
3. The phenomenon independently replicates in Porcellato.
4. The phenomenon independently replicates in Van Beeck after adjustment for dairy and treatment structure.
5. The signal is therefore not confined to the discovery cohorts.
6. Independent PCA can validate the general architecture-resistance phenomenon across heterogeneous studies.
7. The Patangia primary null demonstrates that the relationship is not automatically detectable in every taxonomic dataset.

---

# 8. What the current evidence does not support

The five-study evidence does **not** establish:

- a universal fixed PC1 or PC2 axis;
- universal directionality of independently fitted PCs;
- a single causal bacterial family;
- exact loading transport between cohorts;
- a sequencing-platform-specific biological effect;
- causality of baseline architecture;
- guaranteed prediction in every biological window;
- treatment-effect modification in Van Beeck as part of the external-validation claim.

---

# 9. Statistical guardrails

1. Italy + Manitoba remain the discovery/development framework.
2. Porcellato N = 154 remains its primary external-validation cohort.
3. Van Beeck N = 60 all-trajectory analysis remains its primary external-validation cohort.
4. Van Beeck primary outcome remains overall resistance averaged across the two baseline-relative follow-ups.
5. Van Beeck treatment and dairy remain adjustment covariates, not primary scientific targets.
6. Van Beeck subgroup and treatment-interaction analyses cannot replace the all-60 validation test.
7. Patangia M0→M2 remains the prespecified primary shotgun taxonomic test.
8. Patangia M2→M6 remains exploratory.
9. Use cow-clustered inference when multiple quarter trajectories come from a cow.
10. Do not manufacture cluster-robust inference when each model has one observation per cow/trajectory.
11. Independent PCA validates a phenomenon, not exact axis transport.
12. Cross-study PC signs must not be compared as if independently fitted axes are identical.
13. Exact transportability requires harmonized features and projection using frozen discovery loadings.
14. Do not promote individual PCA loading families to causal biomarkers.
15. Do not infer a 16S-versus-shotgun effect from the current cohorts because platform and biological context are confounded.
16. Do not redefine primary windows or cohorts according to the smallest observed p-value.

---

# 10. Publication-level interpretation

> **Across multiple independent longitudinal bovine mammary microbiome studies, baseline multivariate community organization is associated with subsequent ecological resistance. The relationship is detected in a pooled Italy + Manitoba discovery analysis and independently validated in Porcellato and Van Beeck, despite substantial differences in cohort structure and study context. A prespecified taxonomic analysis in the smaller Patangia shotgun-metagenomic cohort is null, demonstrating that the phenomenon is not universally detectable and motivating explicit study of the ecological and analytical contexts under which baseline architecture carries predictive information.**

---

# 11. Provenance

The earlier document:

`docs/four_study_ecological_resistance_freeze_2026-08-16.md`

is retained as a historical snapshot of the evidence state **before Van Beeck external validation was completed**. It should not be used as the current evidence summary.
