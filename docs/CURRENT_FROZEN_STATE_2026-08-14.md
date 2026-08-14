# MMER current frozen state — 2026-08-14

This file is the short pointer to the current analysis state. Detailed historical analyses remain in their dated documents.

## Study identities

- Study 1 — Biscarini et al. 2020 / PRJEB38332 — Italy
- Study 2 — Van Beeck et al. / PRJEB63336 — California, USA
- Study 3 — Manitoba / PRJNA448966 — Manitoba, Canada
- Study 4 — Porcellato et al. / PRJEB35792 — Norway

## Cross-study taxonomic level

Family is the primary cross-study taxonomic level because it is more robust across heterogeneous 16S protocols and matches the ecological resolution used for the current Study 1 reconstruction. Genus remains a secondary biological-interpretation layer.

Top family signatures show clear domain heterogeneity. Study 1 is enriched for Pseudomonadaceae (17.34%), Streptococcaceae (12.17%), Alcaligenaceae (10.17%), and Propionibacteriaceae (6.36%). Study 2 is led by Staphylococcaceae (11.42%), Peptostreptococcaceae (11.16%), and Lachnospiraceae (7.25%). Study 3 is diffuse and led by Lachnospiraceae (10.96%), Oscillospiraceae (6.24%), Peptostreptococcaceae (5.99%), and Comamonadaceae (5.69%). Study 4 is concentrated around Staphylococcaceae (20.47%), Aerococcaceae (18.25%), and Corynebacteriaceae (16.34%).

Families required to reach 90% cumulative mean abundance:

- Study 1: 36
- Study 2: 35
- Study 3: 45
- Study 4: 19

## Study 1 family-level Q1

Current family-level Study 1 Q1 does not support a statistically significant overall treatment/quarter-context effect on Bray perturbation magnitude.

Friedman P values:

- T1→T2: 0.0624
- T1→T3: 0.0772
- mean displacement: 0.1447

Cow-adjusted treatment P values:

- T1→T2: 0.1767
- T1→T3: 0.1483
- mean displacement: 0.1306

Cloxacillin-assigned quarters show the largest observed mean displacement (0.582), but treatment is confounded with anatomical quarter and is not interpreted as a separable causal effect.

## Study 1 family-level Q3

Frozen alpha-ecology direction:

`+ richness, + Shannon, + evenness, - dominance`

Mean-perturbability univariate results:

- richness: beta +0.111, q = 2.83e-4, R² = 0.549
- Shannon: beta +0.112, q = 2.83e-4, R² = 0.560
- evenness: beta +0.110, q = 2.83e-4, R² = 0.542
- dominance: beta -0.109, q = 2.83e-4, R² = 0.528

The association is substantially stronger for T1→T3 than for T1→T2.

Alpha/treatment nested leave-one-cow-out:

| Model | RMSE | Pearson | Spearman | CV R² |
|---|---:|---:|---:|---:|
| Null | 0.14914 | -0.299 | -0.307 | -0.050 |
| Treatment only | 0.14874 | 0.269 | 0.313 | -0.045 |
| Alpha4 | 0.13802 | 0.508 | 0.612 | 0.101 |
| Alpha4 + treatment | 0.13301 | 0.550 | 0.605 | 0.165 |

## Frozen nested PCA extension

PCA is fitted inside each LOCO training fold on CLR-transformed baseline family composition, with the held-out cow projected into the training-fold coordinate system.

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

PC1+PC2 improves RMSE by 15.0% relative to Alpha4+treatment and raises CV R² by 0.232.

## Frozen interpretation

For Study 1:

`baseline family-level community state > alpha summaries > treatment alone`

The current refined hypothesis is:

`Perturbability = f(baseline community identity/state, ecological structure, perturbation context)`

The PC1+PC2 result is frozen as hypothesis-generating discovery evidence. No additional PCs, taxa-selected predictors, interactions, or nonlinear models should be added in Study 1 solely to improve performance.

## Canonical detailed record

See `docs/study1_family_q1_q3_freeze_2026-08-14.md`.
