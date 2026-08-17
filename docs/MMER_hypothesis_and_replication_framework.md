# MMER hypothesis and replication framework

Updated: 2026-08-16

## Scientific objective

MMER tests whether the **ecological state of the bovine mammary microbiome before a major longitudinal transition contains information about its subsequent ecological resistance**.

Ecological resistance is defined as:

`ecological_resistance = 1 - Bray-Curtis(baseline, follow-up)`

Higher values indicate greater compositional retention relative to baseline.

The current biological model is:

`Ecological resistance = f(baseline community architecture, ecological context, biological window, background)`

rather than ecological response being determined solely by treatment exposure or by one universal diversity rule.

This document freezes the distinction between **development/discovery**, **external validation**, **cross-platform generalization**, and **exploratory follow-up**.

---

## Current study roles

### Development/discovery: Italy + Manitoba

The current frozen primary healthy analysis combines two independent 16S cohorts:

- **Italy / Biscarini et al. / PRJEB38332**
- **Manitoba / Derakhshani et al.**

Together they contribute **49 healthy mammary-quarter trajectories from 14 cows**:

- Italy: 20 trajectories / 5 cows
- Manitoba: 29 trajectories / 9 cows

These two cohorts define the development-stage hypothesis that **baseline multivariate family-level community architecture predicts subsequent ecological resistance**.

### Independent external 16S validation: Porcellato

Porcellato is an independently reprocessed same-quarter longitudinal 16S cohort.

The locked primary analysis contains **154 Sampling1→Sampling2 quarter trajectories** after requiring at least 1,000 DADA2 non-chimeric reads at both timepoints.

This dataset tests whether the general architecture-resistance phenomenon replicates outside the discovery cohorts.

### External cross-platform generalization: Patangia

Patangia is a shotgun-metagenomic cohort with **24 cows sampled at M0, M2, M4, and M6**. Family-level taxonomic profiles were generated with Kraken2.

Patangia tests whether a baseline taxonomic architecture-resistance relationship is detectable under a different data modality and biological context.

### Contextual/sensitivity cohorts

Van Beeck and Wisconsin remain scientifically useful, but they are **not part of the frozen four-study primary architecture evidence** summarized in this framework.

---

## Discovery analysis: Italy + Manitoba

### Baseline representation

Baseline family-level profiles are:

1. prevalence filtered;
2. transformed using compositional log-ratio methodology;
3. centered within study to reduce cohort-centroid effects;
4. summarized using pooled PCA.

The primary predictor is the **joint multivariate baseline architecture represented by PC1 and PC2**, not an isolated taxon or a single diversity index.

### Primary model

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Cow-clustered CR2 inference is used because individual cows contribute multiple quarter trajectories.

### Frozen discovery result

- joint PC1 + PC2 HTZ F = **10.0**
- df = **2, 7.81**
- p = **0.00699**
- PC1 beta = **-0.0286**, p = 0.124
- PC2 beta = **+0.0622**, p = **0.00302**
- descriptive R² = **0.334**
- adjusted R² = **0.289**

The study interaction test was not conventionally significant (`p = 0.0766`).

### Interpretation

The pooled Italy + Manitoba analysis supports a shared association between baseline multivariate community architecture and subsequent ecological resistance while allowing cohort-specific heterogeneity.

This is the signal taken forward for external testing.

---

## External validation hierarchy

MMER distinguishes several forms of replication.

### Level I — phenomenon-level replication

An independently processed cohort shows that its own baseline multivariate community architecture is jointly associated with subsequent ecological resistance.

### Level II — robustness within the external cohort

The association remains under reasonable prespecified or biologically justified sequencing-depth and quality thresholds.

### Level III — exact axis transportability

Frozen discovery loadings are projected into a harmonized external feature space and retain predictive/inferential value.

The current Porcellato result establishes **Level I and strong Level II evidence**. It does not yet establish Level III exact-axis transportability because Porcellato PCA was learned independently.

---

## Porcellato external 16S validation

### Locked cohort

Trajectory reconstruction produced **189 complete same-quarter Sampling1→Sampling2 trajectories**.

The primary pre-model depth rule retained **154 trajectories**:

- Farm A: 93
- Farm K: 61

### Taxonomic processing

DADA2/SILVA processing produced:

- 378 samples
- 15,567 ASVs
- 213 bacterial families overall

The locked primary resistance analysis used:

- 308 samples
- 211 families
- 50 baseline families retained for PCA

PCA variance explained:

- PC1 = **19.779%**
- PC2 = **8.823%**

### Primary model

`ecological_resistance ~ z_PC1 + z_PC2 + farm`

Cow-clustered CR2 inference is appropriate because cows can contribute more than one quarter trajectory.

### Primary result

- PC1 beta = **-0.0892**, SE = 0.0227, p < 0.001
- PC2 beta = **+0.0636**, SE = 0.0228, p = **0.0106**
- farm K beta = **+0.3061**, p < 0.001
- joint PC1 + PC2 HTZ F = **9.49**
- df = **2, 23.5**
- p < **0.001**
- descriptive R² = **0.410**

### Robustness

Assigned bacterial family-read thresholds preserved the conclusion:

- >=1000 reads at both timepoints: `n = 143`, joint F = **11.9**, p < 0.001
- >=2000 reads at both timepoints: `n = 121`, joint F = **7.87**, p = **0.00278**

### Interpretation

Porcellato provides strong independent 16S validation of the **general phenomenon** that baseline multivariate community architecture is associated with subsequent same-quarter ecological resistance.

The finding must not be described as proof that the exact Italy + Manitoba PC axes are universally conserved.

---

## Patangia shotgun-metagenomic generalization

### Cohort and matrix

Patangia contains:

- 24 cows
- M0, M2, M4, and M6 for every cow
- 96 samples total
- 852 family-level Kraken2 taxa overall

Validated matrix QC showed no non-finite relative abundances and sample row sums equal to 1.

At M0:

- 780 families were observed
- 680 remained after the >=10% prevalence filter
- PC1 variance explained = **8.958%**
- PC2 variance explained = **6.580%**

### Primary M0→M2 test

Mean ecological resistance = **0.441**.

- PC1 beta = **+0.0334**
- PC2 beta = **-0.0252**
- descriptive R² = **0.0275**
- no significant joint PC1 + PC2 association

### Interpretation

The prespecified Patangia taxonomic architecture test is **negative**. This result should remain a genuine external null rather than being reframed after inspecting later time windows.

### Additional intervals

M2→M4:

- mean resistance = 0.576
- R² = 0.0506
- standard joint p = 0.5797
- HC3 joint p = 0.6360

M4→M6:

- mean resistance = 0.530
- R² = 0.0819
- standard joint p = 0.4079
- HC3 joint p = 0.3925

Exploratory M2→M6:

- mean resistance = 0.503
- R² = 0.1448
- standard joint p = 0.1936
- HC3 PC2 p = 0.0618
- HC3 asymptotic joint Wald p = 0.0231

Because the ordinary joint test is null, the HC3 result is asymptotic, and only 24 independent cows are available, M2→M6 is **suggestive/exploratory only**.

---

## Statistical inference rules

### Multiple quarter trajectories per cow

When a cohort contains multiple quarter trajectories from the same cow, inference must account for within-cow dependence. The current framework uses cow-clustered CR2 where appropriate.

### One trajectory per cow

For endpoint-specific Patangia models there is one observation per cow. In that setting, cow-clustered CR2 with one observation per cluster is unnecessary and can generate misleadingly small effective degrees of freedom.

Use ordinary independent-observation inference or heteroskedasticity-robust covariance such as HC3 as appropriate.

### PCA interpretation

PCA axes describe multivariate baseline ecological states. Families with large positive or negative loadings define those states but are **not individually inferred to be causal, protective, harmful, or independently significant**.

---

## Current cross-study conclusion

The strongest defensible statement is:

> **Baseline mammary microbial community architecture is associated with subsequent ecological resistance in the pooled Italy + Manitoba development analysis and independently generalizes to Porcellato, while taxonomic architecture does not robustly predict resistance in the smaller Patangia shotgun-metagenomic cohort.**

This evidence supports baseline-state dependence as a biologically meaningful component of mammary ecological resistance while rejecting a simplistic universal-rule interpretation.

The current results do not establish:

- a universal taxonomic axis;
- universal individual-family effects;
- a causal baseline mechanism;
- a sequencing-platform-specific phenomenon;
- generalizability across every temporal window or mammary-health context.

---

## Frozen guardrails

1. Keep Italy + Manitoba as the discovery/development analysis.
2. Keep Porcellato `n = 154` as the primary external 16S validation.
3. Keep Porcellato `n = 143` and `n = 121` as post-taxonomy robustness checks.
4. Keep Patangia M0→M2 as the primary shotgun taxonomic generalization test.
5. Keep Patangia M2→M6 exploratory.
6. Do not infer a 16S-versus-shotgun effect from one shotgun cohort.
7. Do not search additional time windows or taxa solely to rescue a null predefined test.
8. Do not claim exact PC-axis transportability from independent PCA analyses.
9. Do not interpret PCA-loading taxa as causal biomarkers.
10. Preserve non-replication as part of the scientific result.

---

## Publication framing

The four-study design supports a replication-focused ecological story rather than a sequence of independent positive analyses:

> **MMER identifies substantial heterogeneity in longitudinal mammary-community resistance, shows that baseline multivariate architecture predicts this phenotype in two pooled development cohorts, validates the phenomenon independently in Porcellato, and defines an important boundary condition through the primary null result in Patangia.**

Detailed numerical freeze: [`four_study_ecological_resistance_freeze_2026-08-16.md`](four_study_ecological_resistance_freeze_2026-08-16.md).
