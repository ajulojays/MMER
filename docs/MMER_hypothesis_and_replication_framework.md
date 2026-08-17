# MMER hypothesis and replication framework

Updated: 2026-08-16

## Scientific objective

MMER tests whether the **ecological state of the bovine mammary microbiome before a major longitudinal transition contains information about its subsequent ecological resistance**.

Ecological resistance is defined as:

`ecological_resistance = 1 - Bray-Curtis(baseline, follow-up)`

Higher values indicate greater compositional retention relative to the same biological unit's own baseline.

The working biological model is:

`Ecological resistance = f(baseline community architecture, ecological context, biological window, background)`

rather than ecological response being determined solely by treatment exposure or by one universal diversity rule.

This document freezes the distinction between **development/discovery**, **independent external validation**, **cross-platform external generalization**, and **exploratory follow-up**.

---

## Current study roles

### Development/discovery — Italy + Manitoba

The frozen primary healthy analysis combines two independent 16S cohorts:

- **Italy / Biscarini et al. / PRJEB38332**
- **Manitoba / Derakhshani et al.**

Together they contribute **49 healthy mammary-quarter trajectories from 14 cows**:

- Italy: 20 trajectories / 5 cows
- Manitoba: 29 trajectories / 9 cows

These cohorts define the development-stage hypothesis that **baseline multivariate family-level community architecture predicts subsequent ecological resistance**.

### Independent external 16S validation — Porcellato

Porcellato is an independently reprocessed same-quarter longitudinal 16S cohort.

The locked primary analysis contains **154 Sampling1→Sampling2 quarter trajectories** after the fixed sequencing-depth criterion.

Porcellato tests whether the architecture-resistance phenomenon replicates outside the discovery cohorts.

### Independent external 16S validation — Van Beeck

Van Beeck / PRJEB63336 contains a multifactorial three-dairy longitudinal milk dataset sampled at Baseline, 7 Days and 55–75 DIM.

The canonical MMER external-validation cohort contains **60 complete trajectories / 180 samples**.

The frozen Van Beeck question is whether independently learned baseline family-level architecture predicts overall baseline-relative ecological resistance after adjustment for treatment and dairy.

Primary model:

`overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy`

Van Beeck is now part of the primary external-validation evidence structure; it is no longer treated only as a contextual/sensitivity cohort.

### External cross-platform generalization — Patangia

Patangia is a shotgun-metagenomic cohort with **24 cows sampled at M0, M2, M4 and M6**. Family-level taxonomic profiles were generated with Kraken2.

Patangia tests whether a taxonomic baseline architecture-resistance relationship is detectable under a different data modality and biological context.

### Contextual/sensitivity cohort

Wisconsin remains scientifically useful but is not currently part of the frozen five-study primary architecture evidence.

---

## Discovery analysis — Italy + Manitoba

Baseline family-level profiles are:

1. prevalence filtered;
2. compositionally transformed;
3. centered within study to reduce cohort-centroid effects;
4. summarized using pooled PCA.

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Cow-clustered CR2 inference is used because individual cows contribute multiple quarter trajectories.

Frozen discovery result:

- joint PC1 + PC2 HTZ F = **10.0**
- df = **2, 7.81**
- p = **0.00699**
- PC1 beta = **-0.0286**, p = 0.124
- PC2 beta = **+0.0622**, p = **0.00302**
- descriptive R² = **0.334**
- adjusted R² = **0.289**

Interpretation:

> Baseline multivariate community architecture is associated with subsequent ecological resistance across the pooled Italy + Manitoba healthy-quarter analysis.

This is the signal taken forward for external testing.

---

## External-validation hierarchy

MMER distinguishes three forms of replication.

### Level I — phenomenon-level replication

An independently processed cohort shows that its own baseline multivariate community architecture is jointly associated with subsequent ecological resistance.

### Level II — robustness within the external cohort

The association remains under reasonable inferential, sequencing-depth or resampling checks.

### Level III — exact axis transportability

Frozen discovery centering/loading vectors are projected into a harmonized external feature space and retain predictive value.

Current status:

- **Porcellato:** Level I + strong Level II evidence.
- **Van Beeck:** Level I + strong Level II evidence.
- **Level III exact-axis transportability:** not yet established.

Independent PCA validates the phenomenon, not identity of the numerical axes.

---

## Porcellato external 16S validation

Locked cohort:

- N = **154** trajectories
- Farm A = 93
- Farm K = 61
- 50 baseline families retained for PCA

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + farm`

Cow-clustered CR2 is appropriate because cows can contribute more than one quarter trajectory.

Primary result:

- PC1 beta = **-0.0892**, p < 0.001
- PC2 beta = **+0.0636**, p = **0.0106**
- farm K beta = **+0.3061**, p < 0.001
- joint PC1 + PC2 HTZ F = **9.49**
- df = **2, 23.5**
- p < **0.001**
- descriptive R² = **0.410**

Assigned-family-read robustness preserves the conclusion at stricter thresholds.

Interpretation:

> Porcellato provides strong independent 16S phenomenon-level validation of baseline architecture predicting ecological resistance.

---

## Van Beeck external 16S validation

### Locked cohort

- 60 complete longitudinal milk trajectories
- 180 samples
- low-SCC untreated control = 14
- high-SCC untreated control = 12
- CB = 18
- CH = 16

### Frozen phenotype

`overall_resistance = mean[1 - Bray(Baseline,7Days), 1 - Bray(Baseline,55-75DIM)]`

### Baseline architecture

- family-level relative abundance;
- >=10% baseline prevalence filter;
- data-derived zero replacement;
- CLR transformation;
- within-dairy centering;
- independent Van Beeck PCA;
- standardized PC1 + PC2.

### Primary model

`overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy`

The scientific target is the joint PC1 + PC2 contribution. Treatment and dairy are design covariates.

### Primary result

- N = **60**
- full R² = **0.316**
- adjusted R² = **0.224**
- partial R² for PC1 + PC2 = **0.216**
- classical F(2,52) = **7.18**, p = **0.00177**
- HC3 joint Wald χ²(2) = **9.13**, p = **0.0104**
- 9,999 Freedman–Lane permutation p = **0.0183**
- PC1 beta = **-0.0172**, HC3 p = 0.283
- PC2 beta = **-0.0665**, HC3 p = **0.0246**

Bootstrap:

- 5,000 / 5,000 valid stratified complete-trajectory replicates
- PC1 95% interval = **-0.0439 to +0.00665**
- PC2 95% interval = **-0.104 to -0.0217**
- PC2 negative in **99.6%** of replicates
- partial-R² bootstrap median = **0.223**
- percentile interval = **0.0499 to 0.566**

Interpretation:

> Van Beeck independently validates the baseline-architecture → ecological-resistance phenomenon after accounting for its treatment and dairy structure.

The canonical Van Beeck validation intentionally excludes subgroup and treatment-interaction inference from the primary claim.

Detailed workflow: [`vanbeeck_external_validation_workflow_2026-08-16.md`](vanbeeck_external_validation_workflow_2026-08-16.md).

---

## Patangia shotgun-metagenomic generalization

Patangia contains:

- 24 cows
- M0, M2, M4 and M6 for every cow
- 96 samples total
- 852 family-level Kraken2 taxa overall

Prespecified primary M0→M2 result:

- mean ecological resistance = **0.441**
- PC1 beta = **+0.0334**
- PC2 beta = **-0.0252**
- descriptive R² = **0.0275**
- standard joint p ≈ **0.746**
- HC3 joint p ≈ **0.784**

Interpretation:

> The prespecified Patangia taxonomic architecture test is negative.

M2→M4 and M4→M6 are also negative. M2→M6 remains exploratory because conventional and HC3 joint inference disagree and only 24 independent cows are available.

Patangia therefore defines an important boundary condition rather than being optimized into a positive result.

---

## Statistical inference rules

### Multiple quarter trajectories per cow

When a cohort contains multiple quarter trajectories from the same cow, inference must account for within-cow dependence. The framework uses cow-clustered CR2 where appropriate.

### One observation per cow/trajectory

When an endpoint model contains one independent observation per cow or trajectory, singleton-cluster CR2 is not used. Ordinary or heteroskedasticity-robust inference is appropriate.

Van Beeck therefore uses classical nested-model inference, HC3, stratified permutation and bootstrap rather than one-observation-per-core cluster correction.

### PCA interpretation

PCA axes describe multivariate baseline ecological states. Families with large positive or negative loadings define those states but are **not individually inferred to be causal, protective, harmful or independently significant**.

Independent PCA axes can flip sign, rotate or exchange order. Cross-cohort coefficient signs are therefore not biologically comparable unless the same frozen loading basis is projected into each cohort.

---

## Current five-study conclusion

The strongest defensible statement is:

> **Baseline mammary microbial community architecture is associated with subsequent ecological resistance in the pooled Italy + Manitoba development analysis and independently validates in both Porcellato and Van Beeck. The prespecified Patangia shotgun taxonomic analysis is null, indicating that the architecture-resistance relationship is not universally detectable across every cohort, biological window or data modality.**

This supports baseline-state dependence as a reproducible component of mammary ecological resistance while rejecting a simplistic universal-rule interpretation.

The current results do not establish:

- a universal taxonomic axis;
- universal individual-family effects;
- a causal baseline mechanism;
- a sequencing-platform-specific effect;
- exact PC-axis transportability;
- guaranteed generalization across every temporal window or mammary-health context.

---

## Frozen guardrails

1. Keep Italy + Manitoba as discovery/development.
2. Keep Porcellato N = 154 as the primary external 16S validation.
3. Keep Van Beeck N = 60 all-trajectory analysis as the primary external 16S validation.
4. Keep Van Beeck overall resistance as the frozen primary outcome.
5. Keep Van Beeck treatment and dairy as adjustment covariates rather than changing the validation question into a treatment-effect question.
6. Do not replace the Van Beeck all-60 test with subgroup or interaction analyses.
7. Keep Patangia M0→M2 as the primary shotgun taxonomic generalization test.
8. Keep Patangia M2→M6 exploratory.
9. Do not infer a 16S-versus-shotgun effect from the current cohorts.
10. Do not search time windows or taxa solely to rescue a null predefined test.
11. Do not claim exact PC-axis transportability from independent PCA.
12. Do not interpret PCA-loading taxa as causal biomarkers.
13. Preserve non-replication as part of the scientific result.

---

## Publication framing

The five-study design supports a replication-focused ecological story:

> **MMER identifies substantial heterogeneity in longitudinal mammary-community resistance, shows that baseline multivariate architecture predicts this phenotype in pooled Italy + Manitoba development cohorts, validates the phenomenon independently in Porcellato and Van Beeck, and defines an important boundary condition through the prespecified primary null result in Patangia.**

Current numerical freeze: [`five_study_ecological_resistance_freeze_2026-08-16.md`](five_study_ecological_resistance_freeze_2026-08-16.md).
