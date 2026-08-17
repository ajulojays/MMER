# MMER four-study ecological resistance evidence freeze — 2026-08-16

## Purpose

This document freezes the current four-study evidence structure for the MMER hypothesis that **baseline mammary microbial community architecture is associated with subsequent ecological resistance**.

Ecological resistance is defined as:

`R = 1 - Bray-Curtis(baseline, follow-up)`

Higher values indicate greater compositional retention relative to baseline.

This freeze separates:

1. **development/discovery evidence** from Italy + Manitoba;
2. **independent 16S external validation** in Porcellato;
3. **cross-platform shotgun taxonomic generalization** in Patangia.

Van Beeck and Wisconsin are not counted as part of this frozen four-study architecture evidence.

---

## 1. Study roles

| Study | Role | Data type | Longitudinal unit | Status |
|---|---|---|---|---|
| Italy / Biscarini / PRJEB38332 | discovery/development | 16S | mammary-quarter trajectory | positive contribution to pooled primary architecture signal |
| Manitoba / Derakhshani | discovery/development | 16S | mammary-quarter trajectory | positive contribution to pooled primary architecture signal |
| Porcellato | independent external validation | 16S | same-quarter Sampling1→Sampling2 trajectory | strong positive validation |
| Patangia | external cross-platform generalization | shotgun metagenomics, Kraken2 family profiles | cow trajectory | primary null; later interval exploratory only |

The goal is not to require statistical significance in every cohort. The purpose of external generalization is to test where the ecological-resistance phenomenon does and does not transport.

---

## 2. Frozen discovery analysis: Italy + Manitoba

### Cohort

The primary healthy analysis contains **49 quarter trajectories from 14 cows**:

- Italy: 20 trajectories / 5 cows
- Manitoba: 29 trajectories / 9 cows

### Baseline architecture

Family-level baseline compositions were CLR transformed and centered within study before pooled PCA, limiting the extent to which the principal axes simply encode cohort identity.

### Primary model

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Cow-clustered CR2 inference was used because cows contribute multiple quarter trajectories.

### Result

- joint PC1 + PC2 HTZ F = **10.0**
- numerator df = **2**
- denominator df = **7.81**
- p = **0.00699**
- PC1 beta = **-0.0286**, p = 0.124
- PC2 beta = **+0.0622**, p = **0.00302**
- descriptive R² = **0.334**
- adjusted R² = **0.289**

The study interaction test was not conventionally significant (`p = 0.0766`).

### Frozen interpretation

The pooled Italy + Manitoba analysis supports an association between **baseline multivariate community architecture** and subsequent ecological resistance while permitting cohort heterogeneity.

This is the development/discovery signal to be tested externally.

---

## 3. Porcellato: independent 16S external validation

### Trajectory reconstruction and locked depth criterion

Correct parsing of the sample identifiers produced **189 unique Sampling1→Sampling2 quarter trajectories**.

The primary criterion was fixed before model fitting at:

`>= 1000 DADA2 non-chimeric reads at both timepoints`

This retained **154 trajectories**:

- Farm A: 93
- Farm K: 61

Sensitivity sets were:

- >=500 reads: 167 trajectories
- >=2000 reads: 127 trajectories

The locked primary set remains the 154-trajectory cohort.

### Taxonomic processing

DADA2 and SILVA taxonomy produced:

- 378 samples
- 15,567 ASVs
- 213 bacterial families overall

For the locked resistance analysis:

- 308 samples
- 211 families
- 50 families retained for baseline PCA

PCA variance explained:

- PC1 = **0.19779**
- PC2 = **0.08823**

### Resistance by farm

- Farm A: `n = 93`, mean resistance = **0.297**, median = 0.213
- Farm K: `n = 61`, mean resistance = **0.603**, median = 0.668

### Primary external model

`ecological_resistance ~ z_PC1 + z_PC2 + farm`

Cow-clustered CR2 inference was used because cows may contribute multiple quarter trajectories.

### Primary result

- intercept = 0.2968, p < 0.001
- PC1 beta = **-0.0892**, CR2 SE = 0.0227, p < 0.001
- PC2 beta = **+0.0636**, CR2 SE = 0.0228, p = **0.0106**
- farm K beta = **+0.3061**, p < 0.001
- joint PC1 + PC2 HTZ F = **9.49**
- df = **2, 23.5**
- p < **0.001**
- descriptive R² = **0.4100375**

### Post-taxonomy depth robustness

Using assigned bacterial family reads at both timepoints:

- >=500 reads: `n = 154`; joint F = 9.49; p < 0.001
- >=1000 reads: `n = 143`; joint F = 11.9; p < 0.001
- >=2000 reads: `n = 121`; joint F = 7.87; p = **0.00278**

PC1 and PC2 retained the same directions across all three thresholds.

### Frozen interpretation

Porcellato provides **strong independent phenomenon-level replication** that baseline multivariate community architecture is associated with subsequent same-quarter ecological resistance.

Because Porcellato PCA was learned independently, the result does not demonstrate that the exact Italy + Manitoba PC loading vectors transport unchanged. Exact axis transportability requires projection of frozen discovery loadings into a harmonized feature space.

---

## 4. Patangia: shotgun-metagenomic external generalization

### Data structure

Patangia contains **24 cows**, each with M0, M2, M4, and M6 samples:

- 96 samples total
- 24 samples at each timepoint

Shotgun family-level taxonomic profiles were generated from Kraken2 reports.

Validated matrix QC showed:

- 852 families overall
- no NA or non-finite relative abundances
- every sample relative-abundance row sum = 1
- 780 families observed at M0
- 680 families retained at >=10% baseline prevalence

For the M0 baseline PCA:

- PC1 variance explained = **0.08958**
- PC2 variance explained = **0.06580**

### Baseline-relative ecological resistance

| Endpoint | n | Mean | SD | Median | Min | Max |
|---|---:|---:|---:|---:|---:|---:|
| M0→M2 | 24 | 0.441 | 0.252 | 0.357 | 0.150 | 0.831 |
| M0→M4 | 24 | 0.412 | 0.283 | 0.296 | 0.0736 | 0.910 |
| M0→M6 | 24 | 0.381 | 0.206 | 0.298 | 0.121 | 0.801 |

### Prespecified primary M0→M2 test

- PC1 beta = **+0.0334**
- PC2 beta = **-0.0252**
- descriptive R² = **0.02748**
- no significant joint PC1 + PC2 association

**Frozen conclusion:** Patangia does not validate the baseline taxonomic-architecture hypothesis at the prespecified M0→M2 endpoint.

### M2→M4

- mean resistance = **0.576**
- PC1 beta = +0.0555
- PC2 beta = +0.0102
- R² = 0.0506
- standard joint p = 0.5797
- HC3 joint Wald p = 0.6360

Conclusion: negative.

### M4→M6

- mean resistance = **0.530**
- PC1 beta = -0.0696
- PC2 beta = -0.0111
- R² = 0.08185
- standard joint p = 0.4079
- HC3 joint Wald p = 0.3925

Conclusion: negative.

### Exploratory M2→M6

- mean resistance = **0.503**
- PC1 beta = +0.04305
- PC2 beta = +0.07232
- R² = **0.14476**
- standard joint F-test p = **0.1936**
- HC3 PC2 p = **0.0618**
- HC3 asymptotic joint Wald p = **0.0231**

The discrepancy between the ordinary joint test and asymptotic HC3 Wald test, together with only 24 independent cows, means this result remains **suggestive/exploratory**.

### Frozen interpretation

Patangia demonstrates that the taxonomic architecture-resistance relationship is **not uniformly detectable across cohorts, time windows, and data modalities**.

This null external result is informative and must not be converted into a positive result by promoting a later exploratory interval.

It also does not establish that the phenomenon is 16S-specific. Sequencing platform, cohort biology, lactation window, sample size, host/background DNA, and taxonomic resolution are confounded.

A future shotgun-specific hypothesis may test functional metagenomic architecture, but that would be a distinct predefined analysis.

---

## 5. Four-study synthesis

### Supported

The current evidence supports:

> **Baseline mammary microbial community architecture is associated with subsequent ecological resistance in the pooled Italy + Manitoba discovery analysis and independently generalizes in Porcellato.**

### Not universally supported

The current evidence also shows:

> **Taxonomic baseline architecture does not robustly predict ecological resistance in the smaller Patangia shotgun-metagenomic cohort under the prespecified primary test.**

### Therefore

The strongest current model is not a universal fixed taxonomic rule. Instead:

`Ecological resistance = f(baseline community architecture, ecological context, biological window, background)`

with baseline architecture acting as a meaningful but context-dependent component.

---

## 6. Statistical guardrails

1. Italy + Manitoba remain the frozen discovery/development analysis.
2. Porcellato `n = 154` remains the primary external 16S validation set.
3. Porcellato `n = 143` and `n = 121` are post-taxonomy depth robustness analyses.
4. Patangia M0→M2 remains the primary shotgun taxonomic test.
5. Patangia M2→M6 remains explicitly exploratory.
6. When multiple quarter trajectories arise from the same cow, use cow-clustered inference.
7. When an endpoint-specific Patangia model contains one observation per cow, do not manufacture clustering with one-observation clusters; ordinary or heteroskedasticity-robust inference is appropriate.
8. Do not infer a sequencing-platform effect from one shotgun cohort.
9. Do not interpret PCA loading families as individually causal taxa.
10. Independent PCA validates a phenomenon, not exact PC-axis transportability.
11. Exact cross-study transport requires harmonized features and frozen loading projection.
12. Do not search time windows or redefine primary outcomes based on which analysis gives the smallest p-value.

---

## 7. Publication-level interpretation

The scientifically defensible publication story is:

> **Across independent longitudinal bovine mammary microbiome datasets, ecological resistance varies substantially among biological units. Baseline multivariate community architecture predicts this resistance in a pooled two-cohort 16S discovery analysis and replicates strongly in an independent 16S cohort. A smaller shotgun-metagenomic cohort does not reproduce the primary taxonomic association, indicating that the signal is not universally detectable and may depend on ecological context, temporal window, feature representation, or study design.**

This combination of successful replication and principled non-replication is stronger evidence than post hoc optimization for significance in every dataset.
