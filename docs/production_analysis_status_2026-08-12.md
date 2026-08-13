# MMER production analysis status — 2026-08-12

## Production pipeline

The complete 60-sample MMER workflow has now been executed end-to-end on BioProject `PRJEB38332`.

Workflow stages:

1. primer removal / QC
2. DADA2 ASV inference
3. SILVA 138.2 taxonomy assignment
4. bacterial-community ecology layer
5. Q1–Q3 resistance/perturbability analyses
6. Q1/Q2/Q3 inferential hardening and sensitivity analyses

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

Interpretation: the clearest treatment/quarter-associated signal occurs in compositional log-ratio space rather than abundance-weighted Bray–Curtis space.

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

Cloxacillin-assigned trajectories show the most coherent perturbable profile, but individual dimensions do not survive BH-FDR correction across the full Q2 test family.

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

Primary predictors:

- baseline richness
- Shannon diversity
- evenness
- dominance
- CLR-PC1
- CLR-PC2

Primary outcomes:

- Bray–Curtis displacement
- Aitchison displacement (supportive compositional sensitivity outcome)

### Repeated-measures results

Per 1-SD increase in the baseline predictor, effects on later Bray–Curtis displacement were:

| Predictor | beta per SD | BH q |
|---|---:|---:|
| richness | +0.133 | 0.0136 |
| Shannon | +0.155 | 0.00166 |
| evenness | +0.157 | 0.00124 |
| dominance | −0.132 | 0.00375 |
| PC1 | −0.154 | 0.0206 |
| PC2 | −0.102 | 0.0388 |

All six Bray–Curtis associations retained the same direction and BH-FDR significance in the 20-quarter aggregated sensitivity analysis and in the T3-only sensitivity analysis.

### Sequencing-depth diagnostic

Baseline richness was not strongly driven by sequencing depth:

- richness vs DADA2 non-chimeric reads: Spearman rho = 0.297, `P = 0.203`
- richness vs bacterial analysis reads: rho = 0.260, `P = 0.268`
- adjusted depth models were also non-significant

Thus the richness–perturbability relationship is not obviously a sequencing-depth artifact.

### Community configuration

CLR-PC1 and CLR-PC2 both remained predictors of later Bray displacement, indicating that baseline community composition contributes information beyond alpha diversity alone.

Signed PCA loadings identified coherent compositional states. Taxa contributing to the more perturbable multivariate state included UCG-005, *Bacteroides*, *Acinetobacter*, *Flavobacterium*, *Solibacillus*, Rikenellaceae RC9 gut group, *Corynebacterium*, *Alistipes*, *Xylanibacter*, and *Romboutsia*.

Candidate resistance-associated states included *Alcanivorax*, *Kocuria*, *Empedobacter*, *Reyranella*, *Novosphingobium*, *Sphingobium*, *Methylobacterium*, *Paenalcaligenes*, *Klebsiella*, and *Anaerococcus*.

Direct genus-level validation did **not** identify a single genus that survived FDR across repeated-measures, 20-quarter, and T3-only analyses. This supports interpretation of ecological resistance as an emergent whole-community property rather than a single-taxon biomarker.

---

# Current synthesis

The results support the following conceptual model:

`baseline ecological state -> intervention / physiological transition -> multidimensional ecological response`

The major findings are:

1. Natural longitudinal turnover is substantial, consistent with the original Biscarini et al. analysis.
2. Cloxacillin-assigned quarters show the strongest reproducible excess displacement in Aitchison space.
3. The ecological response is multidimensional and includes community membership loss, richness/diversity restructuring, and dominance change, although individual Q2 components do not survive full-family FDR correction.
4. Baseline ecological architecture robustly predicts later Bray–Curtis displacement.
5. Higher baseline richness, Shannon diversity, and evenness predict greater subsequent displacement, whereas greater baseline dominance predicts lower displacement.
6. Baseline multivariate community configuration (PC1/PC2) also predicts later perturbability.
7. No single genus adequately explains the baseline-state effect after stringent sensitivity analysis.

The emerging interpretation is that mammary microbiome resistance is an ecological phenotype of the starting community, not simply a treatment effect or a single-taxon feature.

## Relationship to the original study

Biscarini et al. primarily asked whether dry-cow treatment groups differed in microbiome composition/diversity and reported weak treatment separation relative to strong temporal effects. MMER asks a different set of questions: how far each quarter moves from its own baseline, what ecological dimensions underlie that displacement, and whether the pre-intervention community predicts future perturbability.

These conclusions are complementary rather than contradictory: weak average treatment-group separation can coexist with substantial quarter-specific ecological displacement and baseline-dependent susceptibility to perturbation.
