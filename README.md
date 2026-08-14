# MMER

**Mammary Metataxonomy Ecological Resistance**

MMER is a reproducible reanalysis framework for quantifying **ecological resistance and perturbability of the bovine mammary microbiome** to dry-cow interventions using public longitudinal microbiome data.

## Core scientific questions

1. **Q1 — How far does each mammary community move from its own baseline?**
2. **Q2 — What ecological dimensions constitute resistance/perturbation?**
3. **Q3 — Does baseline ecological state predict later perturbability, and does that relationship generalize across biological contexts?**

For a longitudinal unit `i` at post-baseline time `t`:

`D(i,t) = d(M(i,t), M(i,T1))`

Lower displacement indicates greater ecological resistance.

The project is now explicitly structured as a multi-cohort hypothesis-testing program:

`Study 1 discovery -> Study 2 independent replication/generalization -> refined cross-study hypothesis`

See [`docs/MMER_hypothesis_and_replication_framework.md`](docs/MMER_hypothesis_and_replication_framework.md) for the frozen scientific framework.

---

## Study 1 — Biscarini et al. 2020 / PRJEB38332

Five healthy Holstein-Friesian cows contributed four mammary quarters each and three longitudinal sampling points, giving 60 total samples and 20 quarter trajectories.

### Experimental structure

| Quarter | Anatomical code | Intervention |
|---|---|---|
| Front left | FL / AS | Untreated control |
| Front right | FR / AD | Teat sealant |
| Rear right | RR / PD | Cephalonium |
| Rear left | RL / PS | Cloxacillin |

Timepoints:

- **T1:** dry-off / pre-intervention baseline
- **T2:** calving
- **T3:** 5 days in milk

Treatment is fixed to anatomical quarter, so treatment and quarter anatomy cannot be separated statistically. Results are therefore interpreted as treatment/quarter-condition associations rather than pure causal drug effects.

### Study 1 Q1

Aitchison displacement showed the clearest treatment/quarter-condition signal. Cloxacillin-assigned quarters showed approximately **+20 Aitchison units** of excess displacement relative to untreated controls; the 20-quarter aggregated contrast remained significant after BH correction (`q = 0.019`). Bray–Curtis showed the same broad direction but no significant overall treatment effect.

### Study 1 Q2

Cloxacillin-associated trajectories showed a coherent multidimensional perturbation profile involving reduced core retention, increased turnover, and larger richness/dominance restructuring. Individual Q2 components did not survive full-family BH-FDR correction and are treated as ecological decomposition rather than independent confirmatory endpoints.

### Study 1 Q3 — discovery signal

The portable baseline predictor set is:

- richness
- Shannon diversity
- evenness
- dominance

Study 1 directional signature:

`+ richness, + Shannon, + evenness, - dominance`

The trajectory-level target is:

`mean Bray perturbability = mean[d(T1,T2), d(T1,T3)]`

Nested leave-one-cow-out validation of baseline ecology alone achieved:

- RMSE: **0.10655**
- MAE: **0.08796**
- Pearson `r`: **0.5574**
- Spearman `rho`: **0.5895**
- CV `R²`: **0.2941**

Treatment-only prediction was substantially worse (`CV R² = -0.239`), and adding treatment to baseline ecology worsened RMSE by 8.6%.

Study 1 interpretation:

> **Treatment helps define the perturbation context; baseline ecology helps define susceptibility to perturbation.**

Because Study 1 contains only five cows, this is discovery evidence rather than proof of a universal relationship.

Detailed record: [`docs/production_analysis_status_2026-08-12.md`](docs/production_analysis_status_2026-08-12.md).

---

## Study 2 — Van Beeck et al. / PRJEB63336

Study 2 is the independent replication/generalization cohort.

### Longitudinal design

The three milk timepoints are:

- **Baseline:** dry-off, before intramammary treatment
- **7 Days:** 7 days after dry-off/treatment
- **55–75 DIM:** 55–75 days in milk in the subsequent lactation

This temporal structure is **not equivalent** to Study 1. In particular, Study 1 T2 is calving whereas Study 2 T2 is only 7 days after dry-off.

The source study selected a high-SCC quarter for microbiome sampling, but the currently reconstructed deposited metadata do not explicitly encode anatomical quarter identity across all three longitudinal samples. Exact same-quarter continuity will not be claimed until directly verified.

### Metadata reconstruction

Recovered project inventory:

- 530 sequencing runs
- 514 unique BioSamples
- 208 longitudinal core IDs
- 114 complete 3-timepoint trajectories across sample types
- 67 complete milk trajectories
- 47 complete teat-skin trajectories

The complete milk set contains 201 biological samples. Seven BioSamples had two technical sequencing runs; technical runs were collapsed before DADA2.

### DADA2 and taxonomy

Study 2 single-end Ion Torrent V4 production processing produced:

- input reads: 8,584,765
- non-chimeric reads: 4,068,186
- final ASVs: 14,102
- bacterial ASVs: 12,561
- bacterial reads: 3,856,793
- median bacterial depth: 16,255
- genus-assigned ASVs: 8,980 / 14,102 (63.7%)

SILVA 138.2 NR99 genus-level taxonomy was used.

### Primary sequencing-depth decision

Primary Study 2 analyses use **3,000 bacterial reads/sample**, retaining:

- **60 of 67 complete trajectories (89.6%)**
- **180 milk samples**

Sensitivity thresholds are 2,000 and 4,000 reads/sample.

At 3k:

| Dairy | CB | CH | control high SCC | control low SCC | Total |
|---|---:|---:|---:|---:|---:|
| Dairy 1 | 5 | 3 | 2 | 4 | 14 |
| Dairy 2 | 6 | 4 | 3 | 3 | 16 |
| Dairy 3 | 8 | 8 | 7 | 7 | 30 |
| **Total** | **19** | **15** | **12** | **14** | **60** |

At 4k, 54 trajectories remain; at 5k, only 47 remain.

---

## Study 2 Q1 — displacement magnitude

ASV-level mean Bray–Curtis displacement:

- Baseline → 7 Days: **0.9013**
- Baseline → 55–75 DIM: **0.9122**
- 7 Days → 55–75 DIM: **0.9173**

Repeated mixed model:

`displacement ~ interval + treatment + dairy + interval:treatment + (1|core)`

Fixed-effect tests:

- interval: `P = 0.981`
- treatment: `P = 0.782`
- dairy: `P = 0.304`
- interval × treatment: `P = 0.883`

Thus large longitudinal ecological displacement occurred, but treatment did not explain the magnitude or temporal pattern of individual baseline-relative displacement after accounting for dairy.

---

## Study 2 Q3 — primary replication result

The primary Q3 target preserves the Study 1 structure:

`mean baseline displacement = mean[D(Baseline,7 Days), D(Baseline,55–75 DIM)]`

The same frozen baseline alpha-ecology predictors were used: richness, Shannon diversity, evenness, and dominance.

### Primary ASV endpoint

The **Study 1 directional rule did not replicate**.

Study 1:

`+ richness, + Shannon, + evenness, - dominance`

Study 2 ASV target:

`- richness, - Shannon, - evenness, + dominance`

For the primary ASV mean-displacement target, none of the four features survived within-outcome BH-FDR correction.

Leave-one-dairy-out prediction also failed to show portable ASV-level alpha ecology:

| Model | RMSE for mean displacement |
|---|---:|
| Null mean | **0.07009** |
| Treatment | 0.07041 |
| Alpha ecology | 0.07151 |
| Ecology + treatment | 0.07283 |

ASV CLR/PCA baseline-composition models likewise failed to beat the null across dairies.

**Primary conclusion:** the strong universal-direction version of the Study 1 Q3 hypothesis is not replicated in Study 2 at the ASV-level endpoint.

---

## Study 2 genus-level sensitivity analysis

Because ASV-level Bray was strongly compressed near the upper bound, a secondary taxonomic-resolution sensitivity analysis aggregated bacterial ASVs to genus/higher-taxonomy bins before calculating Bray displacement.

Genus-level mean displacement was substantially lower:

- Baseline → 7 Days: **0.6804**
- Baseline → 55–75 DIM: **0.7297**
- mean baseline displacement: **0.7051**

Yet ASV and genus displacement remained strongly correlated. For mean baseline displacement:

- Pearson `r = 0.821`
- Spearman `rho = 0.766`

Thus ASV resolution amplified fine-scale turnover but tracked the same underlying ecological movement.

### Genus-level Q3 association

At genus resolution, baseline ecology strongly associated with future displacement after dairy and treatment adjustment, **but in the opposite direction from Study 1**.

For genus mean baseline displacement:

| Feature | Direction | BH q | Incremental R² |
|---|---|---:|---:|
| richness | negative | 0.03497 | 0.0693 |
| Shannon | negative | 0.00651 | 0.1196 |
| evenness | negative | 0.00101 | 0.1756 |
| dominance | positive | 0.00101 | 0.1833 |

This supports the weaker statement that baseline ecological state relates to perturbability, but not a universal direction of susceptibility.

### Genus-level transportability

High-dimensional genus CLR/PCA composition did not generalize across dairies.

Simple baseline alpha ecology showed modest, inconsistent portability for genus mean displacement:

- null RMSE: **0.14572**
- alpha-ecology RMSE: **0.13881**
- improvement versus null: **4.74%**
- pooled CV `R²`: **0.009**

This is not yet a strong predictive model.

---

## Low-SCC untreated control sensitivity

A targeted sensitivity analysis restricted Q3 to the **14 untreated low-SCC complete trajectories**:

- Dairy 1: 4
- Dairy 2: 3
- Dairy 3: 7

The reverse Study 2 direction persisted.

For ASV mean displacement:

- evenness: beta `-1.157`, `P = 0.033`, incremental `R² = 0.325`
- dominance: beta `+0.967`, `P = 0.057`, incremental `R² = 0.272`

For genus mean displacement, all four directions remained:

`- richness, - Shannon, - evenness, + dominance`

Therefore antibiotic exposure and high-SCC status alone are insufficient explanations for the Study 1 versus Study 2 Q3 difference.

---

## Current cross-study Q3 synthesis

The two studies now give different Q3 directional results despite both using dry-off pre-treatment milk as baseline.

### Strong universal hypothesis

`One universal baseline alpha-diversity rule predicts mammary perturbability with the same direction across cohorts.`

**Current result: not supported.**

### Weaker ecological-susceptibility hypothesis

`Baseline ecological state contains information about subsequent perturbability.`

**Current result: supported within both studies, but the direction is cohort-dependent and Study 2 cross-dairy portability is weak.**

The combined evidence therefore favors a refined model:

`Perturbability = f(baseline ecological structure, community identity, perturbation regime, ecological background)`

This refined formulation is a hypothesis to test, not a justification for unrestricted model searching within Study 2.

---

## Next hypothesis-driven question

A major mechanistic possibility is that **the organisms present at baseline differ fundamentally between Study 1 and Study 2**.

Alpha diversity measures community shape, not organism identity. Similar richness/Shannon/evenness can describe biologically different taxonomic communities.

The next cross-study question is therefore:

> **Are the baseline mammary microbiomes in Study 1 and Study 2 compositionally different ecological regimes, and could that domain shift explain the opposite Q3 directional relationships?**

The next analysis will compare baseline communities at a harmonized genus level using:

1. shared versus study-specific genera;
2. abundance-weighted shared-genus coverage;
3. baseline genus-composition separation and effect size by study;
4. dominant/core genus overlap;
5. overlap of baseline alpha-ecology support;
6. if justified, a shared-genus harmonized Q3 comparison.

This is a mechanistic follow-up to the failed directional replication, not unrestricted taxa hunting.

Detailed Study 2 record: [`docs/vanbeeck_study2_analysis_status_2026-08-14.md`](docs/vanbeeck_study2_analysis_status_2026-08-14.md).

---

## Reproducibility principles

- Preserve ENA accessions and original aliases.
- Treat repeated samples from the same cow/core ID as longitudinally linked observations.
- Separate ecological resistance from antimicrobial resistance terminology.
- Treat treatment and anatomical quarter as confounded in Study 1.
- Treat Study 1 Q3 as discovery/hypothesis-generating because Study 1 contains only five cows.
- Never evaluate repeated-quarter models with quarter-level random train/test splitting; hold out cows as blocks.
- Learn treatment effects inside training folds when treatment is used predictively.
- For Study 2, explicitly evaluate dairy/background transportability.
- Do not claim exact same-quarter continuity in Study 2 until anatomical-quarter tracking is directly verified.
- Preserve the primary ASV-level replication endpoint; genus analyses remain taxonomic-resolution sensitivities.
- Do not promote secondary results to primary because they are statistically stronger.
- Do not escalate model complexity solely to improve performance on the same replication cohort.
- Do not reinterpret ecological perturbability as disease risk without a separate disease-outcome test.

## Status

**Study 1:** discovery-complete for Q1–Q3.

**Study 2 metadata and sequence processing:** complete for the primary milk cohort.

**Study 2 Q1:** complete at the primary 3k ASV endpoint.

**Study 2 primary Q3 replication:** complete; universal Study 1 direction not replicated.

**Study 2 genus-level sensitivity:** complete; strong ecological association with opposite direction and modest cross-dairy portability.

**Low-SCC untreated sensitivity:** complete; Study 2 reverse direction persists.

**Next:** harmonized Study 1 versus Study 2 baseline genus-community comparison to test whether baseline taxonomic regime/domain shift can explain heterogeneous Q3 behavior.
