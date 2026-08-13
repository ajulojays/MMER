# Study 1 perturbability model findings

Date: 2026-08-13

## Objective

Formalize MMER Q3 as a predictive model of future mammary microbiome perturbability using only pre-perturbation ecological features.

## Modeling unit and target

The biological unit is one cow-quarter trajectory. The discovery cohort contributes 20 trajectories from 5 cows. The target is the mean Bray-Curtis displacement of each quarter from its own T1 baseline across T2 and T3:

`mean_bray = mean[d(T1,T2), d(T1,T3)]`

Baseline predictors are richness, Shannon diversity, evenness, and dominance. Treatment and future measurements are intentionally excluded from the portable perturbability model.

## Prespecified ridge model

The primary CowPAL/MMER v1 architecture uses ridge regression with nested leave-one-cow-out validation and cow-level bootstrap uncertainty.

Observed nested LOCO performance:

| Metric | Value |
|---|---:|
| n trajectories | 20 |
| n cows | 5 |
| RMSE | 0.106 |
| MAE | 0.0875 |
| Pearson r | 0.559 |
| Spearman rho | 0.592 |
| Final ridge lambda | 20.72315 |
| Cow-level bootstrap fits retained | 5000 / 5000 |

Interpretation: using only the baseline ecological state, the model shows moderate out-of-cow prediction of which mammary communities subsequently undergo greater ecological displacement.

## Machine-learning benchmark

A second benchmark compared regularized linear and nonlinear models using the same 20 trajectories, the same four baseline predictors, and nested leave-one-cow-out validation. Hyperparameter tuning occurred only within outer training cows.

| Model | RMSE | MAE | Pearson r | Spearman rho | CV R2 |
|---|---:|---:|---:|---:|---:|
| Elastic Net | 0.10617 | 0.08692 | 0.56632 | 0.60000 | 0.29918 |
| Ridge | 0.10679 | 0.08754 | 0.55646 | 0.58947 | 0.29093 |
| Random Forest | 0.11150 | 0.08740 | 0.54459 | 0.58496 | 0.22712 |
| Gradient Boosting | 0.10899 | 0.08289 | 0.54659 | 0.54032 | 0.26142 |

Elastic Net was the numerical winner by rank prediction, but its improvement over Ridge was very small (Spearman 0.600 vs 0.589). The nonlinear models did not improve held-out-cow rank performance. The current interpretation is therefore that the available predictive signal is largely low-dimensional and can be captured by a parsimonious regularized linear architecture.

## Model-selection decision

Ridge remains the prespecified CowPAL perturbability v1 model until external validation is completed. Elastic Net is retained as a secondary benchmark. A more complex model will only replace Ridge if it demonstrates reproducible improvement in both out-of-cow discovery validation and the independent Manitoba cohort.

## Remaining discovery validation

A cow-blocked permutation null should be run for Ridge and Elastic Net. The four-quarter outcome profile from each cow must be permuted as a block so that within-cow outcome structure is preserved.

## Claim boundary

This model predicts ecological perturbability, not mastitis risk, pathogen burden, treatment efficacy, or clinical outcome. Disease-risk claims require a separate external disease-outcome study.
