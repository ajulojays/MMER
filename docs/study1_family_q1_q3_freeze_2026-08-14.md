# Study 1 family-level Q1/Q3 freeze

Date: 2026-08-14
Status: **FROZEN**

This document records the current canonical Study 1 family-level Q1/Q3 reconstruction and the prespecified nested leave-one-cow-out baseline-composition extension. Earlier Study 1 ridge/alpha-only results are retained in the repository as historical analyses and should not be silently overwritten.

## Cohort and processing

Study 1: Biscarini et al. 2020 / PRJEB38332, Italy.

- 5 cows
- 4 quarter trajectories per cow
- 20 quarter trajectories
- 60 samples total
- T1 = dry-off / pre-intervention baseline
- T2 = calving
- T3 = 5 DIM
- 60/60 samples passed the frozen 3,000 bacterial-read threshold
- family-level bacterial abundance was generated from the existing SILVA 138.2 taxonomy
- each retained sample was rarefied to 3,000 bacterial reads
- trajectory target: `mean[d(T1,T2), d(T1,T3)]` using family-level Bray-Curtis distance

Treatment is structurally confounded with anatomical quarter in Study 1 and is therefore interpreted as treatment/quarter context rather than a separable causal drug effect.

## Q1 — treatment/quarter context and perturbation magnitude

Mean family-level Bray displacement:

| Treatment/quarter context | T1→T2 | T1→T3 | Mean displacement |
|---|---:|---:|---:|
| Untreated | 0.555 | 0.426 | 0.490 |
| Teat sealant | 0.355 | 0.366 | 0.361 |
| Cephalonium | 0.461 | 0.361 | 0.411 |
| Cloxacillin | 0.609 | 0.555 | 0.582 |

Blocked Friedman tests:

| Target | chi-square | P | BH q |
|---|---:|---:|---:|
| T1→T2 | 7.32 | 0.0624 | 0.1158 |
| T1→T3 | 6.84 | 0.0772 | 0.1158 |
| Mean displacement | 5.40 | 0.1447 | 0.1447 |

Cow-adjusted treatment tests were also nonsignificant:

- T1→T2: P = 0.1767
- T1→T3: P = 0.1483
- mean displacement: P = 0.1306

**Frozen Q1 interpretation:** substantial differences in observed perturbation magnitude exist across treatment/quarter contexts, with cloxacillin-assigned quarters showing the largest mean movement, but Study 1 does not provide statistically supported evidence for an overall treatment effect at family-level Bray. Treatment alone is not a useful predictor of held-out-cow perturbability.

## Q3 — baseline family-level ecological structure

The frozen alpha-ecology direction is:

`+ richness, + Shannon, + evenness, - dominance`

For mean family-level perturbability:

| Baseline feature | beta | P | BH q | R² |
|---|---:|---:|---:|---:|
| Richness | +0.111 | 1.86e-4 | 2.83e-4 | 0.549 |
| Shannon | +0.112 | 1.49e-4 | 2.83e-4 | 0.560 |
| Evenness | +0.110 | 2.15e-4 | 2.83e-4 | 0.542 |
| Dominance | -0.109 | 2.83e-4 | 2.83e-4 | 0.528 |

The signal is much stronger for T1→T3 than for T1→T2. For T1→T3, individual baseline alpha features explain approximately 80–85% of trajectory-level variation in displacement. The four alpha metrics are highly correlated descriptors of ecological organization and should not be interpreted as four independent mechanisms.

The four-feature ecology-only multivariable model had:

- R² = 0.5703
- adjusted R² = 0.4557
- model P = 0.00936

Individual coefficients in the joint model were not significant because of strong collinearity among richness, Shannon, evenness, and dominance.

## Frozen alpha/treatment leave-one-cow-out comparison

| Model | RMSE | Pearson | Spearman | CV R² |
|---|---:|---:|---:|---:|
| Null | 0.14914 | -0.299 | -0.307 | -0.050 |
| Treatment only | 0.14874 | 0.269 | 0.313 | -0.045 |
| Alpha4 | 0.13802 | 0.508 | 0.612 | 0.101 |
| Alpha4 + treatment | 0.13301 | 0.550 | 0.605 | 0.165 |

**Interpretation:** baseline ecological structure carries predictive information. Treatment alone does not generalize, but in the alpha-only formulation treatment context adds complementary predictive information.

## Nested baseline-composition extension: PC1 + PC2

A follow-up composition test asked whether the first two family-level baseline compositional axes improve held-out-cow prediction beyond alpha ecology and treatment.

To prevent leakage, PCA was re-fit inside every leave-one-cow-out training fold using only baseline family compositions from the training cows. Held-out cow baselines were projected into the corresponding training-fold PCA space. PCA used CLR-transformed family counts after the frozen 3k rarefaction.

Nested LOCO results:

| Model | RMSE | Pearson | Spearman | CV R² |
|---|---:|---:|---:|---:|
| **PC1 + PC2** | **0.11301** | **0.631** | **0.714** | **0.397** |
| PC1 + PC2 + treatment | 0.11582 | 0.617 | 0.597 | 0.367 |
| Alpha4 + PC1 + PC2 + treatment | 0.13299 | 0.533 | 0.544 | 0.165 |
| Alpha4 + treatment | 0.13301 | 0.550 | 0.605 | 0.165 |
| Alpha4 | 0.13802 | 0.508 | 0.612 | 0.101 |
| Alpha4 + PC1 + PC2 | 0.14058 | 0.477 | 0.546 | 0.067 |
| Treatment only | 0.14874 | 0.269 | 0.313 | -0.045 |
| Null | 0.14914 | -0.299 | -0.307 | -0.050 |

Relative to the frozen Alpha4 + treatment model, PC1 + PC2 alone improved RMSE by **15.0%** and increased CV R² by **0.232**.

Adding treatment to PC1 + PC2 slightly worsened prediction, and stacking alpha metrics, PCs, and treatment did not improve generalization. Given n = 20 trajectories from only 5 cows, the larger combined models are considered over-parameterized exploratory models.

## Frozen Study 1 interpretation

The current Study 1 hierarchy is:

`baseline family-level community state > alpha-ecological summaries > treatment alone`

The strongest held-out-cow signal is carried by the first two baseline family-composition axes. This refines the biological hypothesis from a simple diversity rule toward a community-state formulation:

`Perturbability = f(baseline community identity/state, ecological structure, perturbation context)`

For Study 1 specifically, baseline community state appears to dominate predictive performance, while treatment adds little once PC1/PC2 are known.

## Claim boundary

This remains a discovery result from only five cows. The nested PCA procedure prevents direct test-fold leakage but does not solve the small-sample problem. PC1/PC2 should therefore be treated as a **frozen hypothesis-generating extension** that must be tested in independent cohorts rather than optimized further in Study 1.

Do not add additional PCs, taxa-selected predictors, interactions, or nonlinear models to Study 1 solely to improve performance.

## Next confirmatory step

1. Inspect the taxonomic loadings of the frozen PC1/PC2 axes for biological interpretation.
2. Reproduce the same family-level Q1/Q3 framework in Study 2 without changing the Study 1 specification.
3. Test whether family-level baseline composition/community state is portable across ecological domains.
