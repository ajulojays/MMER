# Van Beeck Study 2 analysis status — 2026-08-14

## Purpose

Study 2 uses Van Beeck et al. / PRJEB63336 as an independent replication/generalization cohort for the MMER Q1–Q3 framework developed in Study 1.

The main replication hypothesis is whether pre-treatment baseline mammary ecological state predicts subsequent within-cow microbiome perturbability.

Study 2 is treated as a replication cohort, not an unrestricted feature-discovery dataset.

---

## Longitudinal design

Three milk timepoints are used:

- Baseline: dry-off, before intramammary treatment
- 7 Days: 7 days after dry-off/treatment
- 55–75 DIM: 55–75 days in milk in the subsequent lactation

This is not temporally equivalent to Study 1, where T2 is calving and T3 is 5 DIM.

A further design caveat remains open: the source study selected a high-SCC quarter for microbiome sampling, but the deposited metadata reconstructed so far do not explicitly encode anatomical quarter identity across all three samples. Exact same-quarter continuity should not be claimed until verified directly.

---

## Metadata and sequencing inventory

Reconstructed project inventory:

- 530 ENA sequencing runs
- 514 BioSamples
- 208 longitudinal core IDs
- 114 complete 3-timepoint trajectories across sample types
- 67 complete milk trajectories
- 47 complete teat-skin trajectories

The 67 complete milk trajectories correspond to 201 biological samples.

Technical replicate structure:

- 208 run-level FASTQs attached to the 201 milk BioSamples
- 194 BioSamples with one run
- 7 BioSamples with two runs

Technical replicate runs were collapsed before DADA2.

---

## DADA2 production

Ion Torrent single-end V4 processing used:

- trimLeft = 21
- truncLen = 290
- maxEE = 2
- Ion Torrent homopolymer-aware DADA2 settings
- consensus chimera removal

Production output:

- input reads: 8,584,765
- filtered reads: 4,819,927
- denoised reads: 4,478,375
- non-chimeric reads: 4,068,186
- overall non-chimeric retention: 47.39%
- final ASVs: 14,102
- median ASV length: 269 bp
- minimum non-chimeric depth: 1,394
- median non-chimeric depth: 16,627
- maximum non-chimeric depth: 89,607

---

## SILVA 138.2 taxonomy

Fast genus-level SILVA 138.2 NR99 assignment completed in 1m43s.

Taxonomic assignment:

- Kingdom: 14,084 / 14,102 ASVs (99.9%)
- Phylum: 12,705 (90.1%)
- Class: 12,661 (89.8%)
- Order: 12,497 (88.6%)
- Family: 11,766 (83.4%)
- Genus: 8,980 (63.7%)

Bacterial-only table:

- bacterial ASVs: 12,561 / 14,102
- bacterial reads: 3,856,793
- median bacterial depth: 16,255
- minimum bacterial depth: 1,007
- maximum bacterial depth: 73,947

---

## Sequencing-depth decision

Complete trajectory survival by bacterial depth:

| Minimum reads/sample | Samples retained | Complete trajectories | Trajectory retention |
|---:|---:|---:|---:|
| 1,000 | 201 | 67 | 100.0% |
| 1,500 | 196 | 62 | 92.5% |
| 2,000 | 195 | 61 | 91.0% |
| 2,500 | 195 | 61 | 91.0% |
| 3,000 | 194 | 60 | 89.6% |
| 3,500 | 191 | 57 | 85.1% |
| 4,000 | 188 | 54 | 80.6% |
| 5,000 | 181 | 47 | 70.1% |

The primary Study 2 cohort was fixed at 3,000 bacterial reads/sample, retaining 60 complete trajectories / 180 milk samples.

Primary 3k cohort:

- Dairy 1: 14
- Dairy 2: 16
- Dairy 3: 30

Treatment/control counts:

- CB: 19
- CH: 15
- control high SCC: 12
- control low SCC: 14

Dairy × treatment:

| Dairy | CB | CH | control high SCC | control low SCC | Total |
|---|---:|---:|---:|---:|---:|
| Dairy 1 | 5 | 3 | 2 | 4 | 14 |
| Dairy 2 | 6 | 4 | 3 | 3 | 16 |
| Dairy 3 | 8 | 8 | 7 | 7 | 30 |
| Total | 19 | 15 | 12 | 14 | 60 |

At 4k, 54 trajectories remain. At 5k, only 47 remain.

---

## Q1 — ASV-level longitudinal displacement

Mean Bray–Curtis displacement:

- Baseline → 7 Days: 0.9013
- Baseline → 55–75 DIM: 0.9122
- 7 Days → 55–75 DIM: 0.9173
- mean baseline displacement: 0.9068

Treatment means for Baseline → 7 Days ranged from approximately 0.885 to 0.913. Dairy means ranged from approximately 0.882 to 0.938.

Individual adjusted models found no significant overall treatment effect for baseline-relative displacement.

Repeated mixed model:

`displacement ~ interval + treatment + dairy + interval:treatment + (1|core)`

Fixed-effect tests:

- interval: P = 0.9810
- treatment: P = 0.7816
- dairy: P = 0.3035
- interval × treatment: P = 0.8832

Likelihood-ratio tests were similarly null.

Baseline-to-7d versus baseline-to-late displacement showed Pearson r = 0.433, P = 0.00055, but Spearman rho = 0.192, P = 0.143. This is interpreted cautiously as moderate linear persistence without strong rank-stable perturbability memory.

Q1 conclusion: large ecological movement occurs, but treatment does not explain the magnitude or temporal pattern of individual baseline-relative displacement in this cohort after accounting for dairy.

---

## Q3 — primary ASV-level replication

The frozen Study 1-style baseline predictors were:

- richness
- Shannon diversity
- evenness
- dominance

Primary Study 2 target:

`mean baseline displacement = mean[D(Baseline,7 Days), D(Baseline,55–75 DIM)]`

Adjusted inference model:

`outcome ~ baseline ecological feature + dairy + treatment`

### ASV-level result

For the primary mean-displacement target:

- richness: beta negative, P = 0.533
- Shannon: beta negative, P = 0.394
- evenness: beta = -0.379, P = 0.0499, BH q = 0.168
- dominance: beta = +0.369, P = 0.0840, BH q = 0.168

The Study 2 coefficient pattern is the opposite of Study 1:

Study 1:

`+ richness, + Shannon, + evenness, - dominance`

Study 2:

`- richness, - Shannon, - evenness, + dominance`

Therefore the strong universal-direction replication hypothesis fails at the primary ASV endpoint.

### ASV alpha-ecology leave-one-dairy-out

For mean baseline displacement:

- null RMSE: 0.07009
- treatment RMSE: 0.07041
- ecology RMSE: 0.07151
- ecology + treatment RMSE: 0.07283

Ecology was 2.02% worse than null by RMSE and pooled CV R² remained negative.

### ASV CLR/PCA leave-one-dairy-out

Baseline ASV composition also failed to beat null reliably.

For mean baseline displacement:

- null RMSE: 0.07009
- CLR/PCA RMSE: 0.07046
- RMSE change versus null: -0.53%
- CV R²: -0.064

No further model escalation is justified on the same 60 trajectories.

---

## Secondary genus-level resolution analysis

The primary ASV target showed strong ceiling compression, so ASVs were aggregated to genus/higher-taxonomy bins as a taxonomic-resolution sensitivity analysis.

Genus table:

- 1,032 genus/higher-taxonomy bins before 3k rarefaction
- 916 bins retained after rarefaction

Genus displacement means:

- Baseline → 7 Days: 0.6804
- Baseline → 55–75 DIM: 0.7297
- 7 Days → late: 0.7190
- mean baseline displacement: 0.7051

ASV and genus displacement remained strongly correlated:

| Target | Pearson r | Spearman rho |
|---|---:|---:|
| Baseline → 7 Days | 0.811 | 0.744 |
| Baseline → late | 0.767 | 0.719 |
| Mean baseline displacement | 0.821 | 0.766 |

Thus ASV resolution magnified the magnitude of turnover but did not invent a separate ecological phenomenon.

---

## Genus-level Q3

At genus level, baseline ecological state strongly associated with future displacement after adjustment for dairy and treatment.

For genus mean baseline displacement:

| Feature | Beta direction | P | BH q | Incremental R² |
|---|---|---:|---:|---:|
| richness | negative | 0.03497 | 0.03497 | 0.0693 |
| Shannon | negative | 0.00488 | 0.00651 | 0.1196 |
| evenness | negative | 0.000507 | 0.00101 | 0.1756 |
| dominance | positive | 0.000367 | 0.00101 | 0.1833 |

This is statistically strong but directionally opposite to Study 1.

Interpretation: Study 2 supports a weaker ecological-susceptibility hypothesis — baseline state matters — but not a universal direction of susceptibility.

### Genus CLR/PCA leave-one-dairy-out

High-dimensional genus composition did not generalize:

- null RMSE: 0.14572
- genus CLR/PCA RMSE: 0.15394
- RMSE change versus null: -5.65%
- pooled CV R²: -0.219

### Genus alpha-ecology leave-one-dairy-out

Simple alpha ecology showed modest, inconsistent portability:

For mean genus displacement:

- null RMSE: 0.14572
- alpha-ecology RMSE: 0.13881
- RMSE improvement: 4.74%
- pooled CV R²: 0.009
- Pearson r: 0.265

For Baseline → 7 Days, RMSE improved by 6.97% versus null but pooled CV R² remained slightly negative. For Baseline → late, alpha ecology was 3.95% worse than null.

Conclusion: broad ecological organization may be somewhat more transportable than exact taxonomic composition, but predictive generalization remains weak and dairy-dependent.

---

## Low-SCC untreated control sensitivity

A targeted sensitivity analysis restricted Q3 to untreated low-SCC cows to address the possibility that the Study 2 reversal was driven by antibiotics or high baseline SCC.

At 3k:

- n = 14 complete trajectories
- Dairy 1 = 4
- Dairy 2 = 3
- Dairy 3 = 7

The reverse Study 2 direction remained.

ASV mean displacement:

- evenness: beta = -1.157, P = 0.0333, incremental R² = 0.325
- dominance: beta = +0.967, P = 0.0570, incremental R² = 0.272

Genus mean displacement retained the same directions for all four features:

`- richness, - Shannon, - evenness, + dominance`

For genus Baseline → 7 Days:

- richness: beta = -0.00106, P = 0.0255
- Shannon: beta = -0.140, P = 0.0447

Because n = 14, this is interpreted by direction and effect magnitude rather than treated as an independently powered subgroup study.

The key conclusion is that the Study 2 reversal does not disappear in low-SCC untreated cows. Antibiotic treatment and high-SCC status alone therefore do not explain the difference from Study 1.

---

## Cross-study interpretation

Study 1 and Study 2 produce different Q3 directions despite sharing a dry-off pre-treatment baseline.

Study 1 discovery pattern:

`higher richness/Shannon/evenness -> greater displacement; higher dominance -> lower displacement`

Study 2 pattern:

`higher richness/Shannon/evenness -> lower displacement; higher dominance -> greater displacement`

The studies also differ in the biological transition measured:

- Study 1: dry-off → calving → 5 DIM
- Study 2: dry-off → 7 days → 55–75 DIM

High SCC is an important Study 2 design factor, but targeted low-SCC untreated analysis shows that SCC/treatment alone are insufficient to explain the reversal.

A major remaining mechanistic possibility is that the organisms/community regimes present at baseline differ fundamentally between Study 1 and Study 2. Alpha-diversity summaries can be numerically similar while representing biologically different communities.

---

## Next hypothesis-driven analysis

Before moving to another predictive model, compare the baseline taxonomic ecology of Study 1 and Study 2 at a harmonized genus level.

The next question is:

> Are Study 1 and Study 2 starting from fundamentally different baseline mammary community regimes, and could this explain the opposite Q3 directional relationships?

Planned tests:

1. shared versus study-specific genera;
2. abundance-weighted shared-genus coverage;
3. baseline genus-composition separation and effect size by study;
4. dominant/core genus overlap;
5. support overlap for richness, Shannon, evenness, and dominance;
6. if appropriate, Q3 comparison in a harmonized shared-genus representation.

This is a mechanistic follow-up to failed directional replication, not unrestricted exploratory taxa mining.
