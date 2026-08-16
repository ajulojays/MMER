# Frozen primary healthy-quarter analysis — 2026-08-16

## Scope

The primary inferential analysis is now restricted to two phenotype-comparable healthy-quarter cohorts:

- **Italy:** 20 quarter trajectories from 5 cows
- **Manitoba:** 29 quarter trajectories from 9 cows
- **Total:** 49 quarter trajectories from 14 cows

**Wisconsin is excluded from the primary analysis** because individual quarter-level Healthy/CHRON/NEWINF/POS labels could not be recovered from the available public metadata. It should not be described as a healthy validation cohort.

## Frozen phenotype

For quarter trajectory `iq`:

`ecological resistance = 1 - BrayCurtis(T1_iq, T2_iq)`

Higher values indicate greater compositional retention across the dry-off to early-postpartum transition; lower values indicate greater ecological displacement/susceptibility.

Mean resistance:

- Italy: **0.555**
- Manitoba: **0.448**

## Frozen baseline ecological architecture

Baseline (`T1`) bacterial family composition was:

1. filtered at 10% prevalence;
2. CLR transformed with pseudocount 0.5;
3. centered within study to remove cohort centroids;
4. decomposed by pooled PCA across Italy + Manitoba.

A total of **161 families** entered the PCA.

- PC1 variance: **10.03%**
- PC2 variance: **9.26%**
- PC1 + PC2: **19.29%**

## Aim 1 — baseline architecture predicts ecological resistance

Primary model:

`ecological_resistance ~ z_PC1 + z_PC2 + study`

Inference uses cow-clustered CR2 small-sample robust variance estimation.

Joint PC1 + PC2 test:

- **F(2, 7.81) = 10.0**
- **P = 0.00699**

Robust coefficients:

| Term | beta | CR2 SE | Satterthwaite P |
|---|---:|---:|---:|
| PC1 (per SD) | -0.0286 | 0.0164 | 0.124 |
| PC2 (per SD) | +0.0622 | 0.0152 | **0.00302** |
| Manitoba vs Italy | -0.107 | 0.0299 | **0.00625** |

Descriptive model fit:

- R² = **0.334**
- adjusted R² = **0.289**

The pooled architecture effect is therefore driven predominantly by PC2.

## Aim 2 — reproducibility across Italy and Manitoba

Interaction model:

`ecological_resistance ~ (z_PC1 + z_PC2) * study`

Joint PC × study interaction:

- **F(2, 5.55) = 4.23**
- **P = 0.0766**

Interpretation: there is insufficient evidence at the conventional 0.05 threshold that the architecture-resistance relationship differs between Italy and Manitoba, but the magnitude and dominant axis are not identical across cohorts and context dependence remains plausible.

Study-specific descriptive pattern:

- Italy PC2 association is strong and positive; PC2-resistance correlation = **0.735**.
- Manitoba PC2 association is weaker but directionally positive; correlation = **0.180**.

## PC2 ecological-state gradient

Resistance increases monotonically from low to high PC2 in both cohorts.

| Study | Low PC2 | Intermediate PC2 | High PC2 |
|---|---:|---:|---:|
| Italy | 0.422 | 0.600 | 0.656 |
| Manitoba | 0.402 | 0.449 | 0.492 |

The high-PC2 state is therefore interpreted as a **baseline community-organization gradient associated with greater ecological resistance**, not as a single protective taxon.

### Positive PC2 axis-defining families

Prominent positive loadings include:

- Lactobacillaceae
- Pseudomonadaceae
- Propionibacteriaceae
- Alcanivoracaceae1
- Beijerinckiaceae
- Mycobacteriaceae
- Streptococcaceae
- Enterobacteriaceae

### Negative PC2 axis-defining families

Prominent negative loadings include:

- Paracoccaceae
- Saccharimonadaceae
- Bacteroidaceae
- Caulobacteraceae
- Anaerovoracaceae
- Succinivibrionaceae
- Spirochaetaceae
- Ruminococcaceae
- Christensenellaceae
- Oscillospiraceae

These are **axis-associated families**. They should not be described as causal, protective, or risk taxa without separate testing.

## Frozen primary interpretation

> **Baseline ecological architecture predicts the magnitude of microbiome resistance across the dry-off–early postpartum transition in healthy bovine mammary quarters.**

A more mechanistic framing is:

> **Ecological resistance is associated with multivariate baseline community organization rather than simple diversity or preservation of exact quarter identity.**

No further primary outcome/model searching should be performed. Additional work should be limited to descriptive ecological interpretation, visualization, manuscript development, and prospectively motivated validation.

## Publication figure

`scripts/24_publication_figure_primary_healthy.R` generates the frozen main multi-panel figure using only existing primary-analysis outputs. It produces:

- vector PDF
- vector SVG
- 600-dpi PNG

Panels:

A. quarter-level ecological resistance by cohort;
B. baseline PC2 versus ecological resistance with frozen CR2 statistics;
C. low/intermediate/high PC2 resistance gradient in both cohorts;
D. baseline family composition defining the PC2 ecological-state gradient.

The figure script performs **no new hypothesis tests**.
