# MMER production analysis status — updated 2026-08-13

## Production pipeline

The complete 60-sample MMER workflow has now been executed end-to-end on BioProject `PRJEB38332`.

Workflow stages:

1. primer removal / QC
2. DADA2 ASV inference
3. SILVA 138.2 taxonomy assignment
4. bacterial-community ecology layer
5. Q1–Q3 resistance/perturbability analyses
6. Q1/Q2/Q3 inferential hardening and sensitivity analyses
7. nested leave-one-cow-out perturbability modeling and treatment-adjusted predictive sensitivity analysis

### Final DADA2 production output

- Samples: 60
- Input read pairs: 8,816,802
- Filtered: 6,492,116 (73.6%)
- Merged: 5,633,741 (86.8% of filtered)
- Non-chimeric: 5,384,222 (~61.1% of input)
- Non-chimeric ASVs: 9,859
- Median ASV length: 403 bp

### Final taxonomy output

- Bacterial ASVs: 7,736 / 9,859 (78.5%)
- Bacterial reads: 5,074,031 / 5,384,222 (94.2%)
- Genus-assigned reads: 4,358,406 / 5,384,222 (~81.0%)
- Species-assigned reads: 3,025,328 / 5,384,222 (~56.2%)

The primary ecological analyses remain ASV-level, with taxonomy used for biological interpretation.

---

# Scientific questions

## Q1 — How far does each mammary community move from its own baseline?

Q1 quantifies baseline-relative displacement from T1 to T2/T3 using Bray–Curtis and Aitchison distances, with the untreated quarter in the same cow used as the physiological-transition reference.

### Main result

Bray–Curtis displacement showed directional treatment/quarter differences but no significant overall treatment effect in the repeated-measures mixed model (`P = 0.167`).

Aitchison displacement showed a significant overall treatment/quarter-condition effect (`P = 0.018`). The cloxacillin-assigned rear-left quarters produced the most consistent excess displacement relative to untreated front-left quarters.

Key cloxacillin-associated estimates:

- T2 excess Aitchison displacement: +19.8
- T3 excess Aitchison displacement: +20.6
- 20-quarter aggregated estimate: +20.19, `P = 0.0032`, BH `q = 0.0192`
- T3-only estimate: +20.58, `P = 0.0105`, BH `q = 0.0628`
- paired within-cow T3 mean excess: +20.58; 95% CI 10.76–30.41; `P = 0.00435`; BH `q ≈ 0.052`

Interpretation: the clearest treatment/quarter-associated signal occurs in compositional log-ratio space rather than abundance-weighted Bray–Curtis space. Treatment therefore helps define the perturbation context, but it does not by itself explain the substantial between-quarter heterogeneity in ecological displacement.

### Design caveat

Treatment is fixed to anatomical quarter in the source experiment (FL untreated, FR teat sealant, RR cephalonium, RL cloxacillin). Treatment and quarter anatomy therefore cannot be separated statistically. Results should be described as treatment/quarter-condition associations rather than pure causal treatment effects.

---

## Q2 — What ecological dimensions constitute resistance/perturbation?

Q2 decomposes the response phenotype into ecological dimensions rather than treating resistance as a single beta-diversity statistic.

Dimensions include:

- core retention
- membership instability / turnover
- absolute Shannon change
- absolute richness change
- absolute evenness change
- absolute dominance change
- exploratory multidimensional resistance score

### Main result

Cloxacillin-associated trajectories show the most coherent perturbable profile, but individual dimensions do not survive BH-FDR correction across the full Q2 test family.

In the 20-quarter analysis, cloxacillin relative to untreated showed:

- core retention: −0.131 (`P = 0.068`)
- membership instability: +0.131 (`P = 0.068`)
- absolute Shannon change: +0.522 (`P = 0.213`)
- absolute richness change: +102.5 ASVs (`P = 0.093`)
- absolute evenness change: +0.0487 (`P = 0.288`)
- absolute dominance change: +0.0628 (`P = 0.0289`)
- multidimensional resistance: −0.850 (`P = 0.0517`)

After family-wise BH correction, none of these individual dimensions remained significant.

Interpretation: Q2 is best treated as an ecological decomposition of the perturbation phenotype rather than evidence that every component differs independently by treatment/quarter condition.

Important: membership instability is mathematically `1 - core_retention`; these are not independent outcomes. The composite multidimensional score is derived from component metrics and remains exploratory.

---

## Q3 — Does baseline ecological state predict later perturbability?

Q3 tests whether the pre-intervention T1 community predicts later baseline-relative displacement.

### Association-level results

Primary baseline predictors include richness, Shannon diversity, evenness, dominance, CLR-PC1, and CLR-PC2. Per 1-SD increase in the baseline predictor, effects on later Bray–Curtis displacement were:

| Predictor | beta per SD | BH q |
|---|---:|---:|
| richness | +0.133 | 0.0136 |
| Shannon | +0.155 | 0.00166 |
| evenness | +0.157 | 0.00124 |
| dominance | −0.132 | 0.00375 |
| PC1 | −0.154 | 0.0206 |
| PC2 | −0.102 | 0.0388 |

All six Bray–Curtis associations retained the same direction and BH-FDR significance in the 20-quarter aggregated sensitivity analysis and in the T3-only sensitivity analysis.

Baseline richness was not strongly driven by sequencing depth, and direct genus-level validation did not identify a single genus that survived FDR across the major sensitivity analyses. This supports an emergent whole-community interpretation rather than a single-taxon biomarker.

### Portable perturbability phenotype

For predictive modeling, Q3 is summarized at the independent cow-quarter trajectory level as:

`mean Bray perturbability = mean[d(T1,T2), d(T1,T3)]`

This gives 20 quarter-level trajectories from 5 cows. The portable v1 predictor set is deliberately restricted to baseline richness, Shannon diversity, evenness, and dominance.

### Nested leave-one-cow-out prediction

The original ecology-only ridge model was evaluated with nested leave-one-cow-out (LOCO) validation, holding out entire cows rather than individual quarters.

Performance:

- RMSE: **0.10655**
- MAE: **0.08796**
- Pearson `r`: **0.5574**
- Spearman `rho`: **0.5895**
- cross-validated `R²`: **0.2941**

Thus baseline ecological state carries moderate out-of-cow predictive information about subsequent overall Bray–Curtis perturbability in this discovery cohort.

### Treatment-adjusted Q3 sensitivity analysis

A direct nested-LOCO comparison was run using the exact same 20-trajectory target and baseline ecological predictors, with treatment identity added only inside the training folds.

| Model | RMSE | MAE | Pearson r | Spearman rho | CV R² |
|---|---:|---:|---:|---:|---:|
| baseline ecology only | **0.10655** | **0.08796** | **0.5574** | **0.5895** | **0.2941** |
| baseline ecology + treatment | 0.11576 | 0.09824 | 0.4526 | 0.4226 | 0.1668 |
| treatment only | 0.14115 | 0.12554 | -0.2884 | -0.2977 | -0.2387 |

Relative to treatment alone, ecology alone reduced LOCO RMSE by **24.5%**. Adding treatment to baseline ecology worsened RMSE by **8.6%** relative to ecology alone.

Interpretation: in this small discovery cohort, the predictive Q3 signal is not explained by treatment identity. Baseline ecological state predicts overall perturbability better than treatment identity, and treatment does not add out-of-cow predictive value to the ecology-only model.

This does **not** imply that treatment has no ecological effect. Q1 and Q3 answer different questions: Q1 characterizes the perturbation/treatment-quarter context, whereas Q3 asks which starting communities are more susceptible to subsequent displacement.

---

# Current synthesis

The Study 1 results now support the following working model:

`baseline ecological state -> susceptibility to displacement`

within a broader perturbation process:

`baseline ecological state + treatment/physiological context -> observed longitudinal ecological response`

The major findings are:

1. Natural longitudinal turnover is substantial.
2. Cloxacillin-assigned quarters show the strongest reproducible excess displacement in Aitchison space, with the important caveat that treatment is confounded with quarter anatomy.
3. The ecological response is multidimensional, although individual Q2 components do not survive full-family FDR correction.
4. Baseline ecological architecture predicts later Bray–Curtis displacement across association-level sensitivity analyses.
5. Higher baseline richness, Shannon diversity, and evenness predict greater subsequent displacement, whereas greater baseline dominance predicts lower displacement.
6. The portable four-feature baseline ecological model retains moderate predictive performance under nested leave-one-cow-out validation (`RMSE = 0.10655`, `r = 0.557`, `rho = 0.589`, `CV R² = 0.294`).
7. Treatment identity alone does not generalize as a predictor of the Q3 perturbability phenotype (`CV R² = -0.239`).
8. Adding treatment to baseline ecology does not improve Q3 and worsens LOCO RMSE by 8.6%.
9. Therefore the Q3 signal is not simply a proxy for treatment assignment: in Study 1, baseline ecology carries more transferable information about who is perturbable than treatment identity.
10. No single genus adequately explains the baseline-state effect after stringent sensitivity analysis.

A concise interpretation is:

> **Treatment helps define the perturbation context; baseline ecology helps define susceptibility to perturbation.**

## Scope and caution

Study 1 contains only **20 quarter trajectories from 5 cows**. The LOCO results are therefore discovery evidence for transportable ecological susceptibility, not proof of a universally generalizable predictor. External cohorts are required to test whether the relationship survives different farms, populations, sequencing protocols, treatment regimes, and temporal structures.

## Relationship to the original study

Biscarini et al. primarily asked whether dry-cow treatment groups differed in microbiome composition/diversity and reported weak treatment separation relative to strong temporal effects. MMER asks a different set of questions: how far each quarter moves from its own baseline, what ecological dimensions underlie that displacement, and whether the pre-intervention community predicts future perturbability.

These conclusions are complementary rather than contradictory: weak average treatment-group separation can coexist with substantial quarter-specific ecological displacement and baseline-dependent susceptibility to perturbation.

---

# Study 1 decision

Study 1 is now considered **discovery-complete for the core Q1–Q3 hypothesis**. Further model complexity should not be used to extract stronger claims from five cows.

The next priority is cross-study replication using publicly available longitudinal bovine milk microbiome datasets, preserving the central questions while adapting them to each study's design:

- quantify ecological displacement relative to an appropriate baseline;
- model exposure/treatment context where available;
- test whether starting ecological state predicts subsequent displacement;
- use cow-blocked validation whenever repeated quarters/samples permit it;
- evaluate transportability across farms, geography, sequencing protocols, and perturbation types.
