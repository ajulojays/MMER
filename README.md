# MMER

**Mammary Metataxonomy Ecological Resistance**

MMER is a reproducible reanalysis framework for quantifying **ecological resistance and perturbability of the bovine mammary microbiome** to dry-cow interventions using public longitudinal microbiome data.

## Discovery cohort

The first MMER analysis uses **Biscarini et al. (2020)**, BioProject **PRJEB38332**: five healthy Holstein-Friesian cows, four mammary quarters per cow, and three longitudinal sampling points (60 samples total). Each cow contributes all four intervention states, making the design a within-cow, quarter-level experiment.

> ENA labels the material `milk metagenome`; the underlying assay is **16S rRNA-gene amplicon sequencing**, not shotgun metagenomics.

### Experimental structure

| Quarter | Anatomical code | Intervention |
|---|---|---|
| Front left | FL / AS | Untreated control |
| Front right | FR / AD | Teat sealant |
| Rear right | RR / PD | Cephalonium |
| Rear left | RL / PS | Cloxacillin |

Sampling blocks reconstructed from the deposited sample order:

- **S1–S20:** T1, dry-off / baseline
- **S21–S40:** T2, calving
- **S41–S60:** T3, 5 days in milk

## Core scientific questions

1. **Q1 — How far does each mammary community move from its own baseline?**
2. **Q2 — What ecological dimensions constitute resistance/perturbation?**
3. **Q3 — Does baseline ecological state predict later perturbability?**

For quarter `q` of cow `i` at post-baseline time `t`:

`D(i,q,t) = d(M(i,q,t), M(i,q,T1))`

Lower displacement indicates greater ecological resistance.

## Study 1 — current findings

The complete 60-sample production workflow has been executed end-to-end.

### Q1 — perturbation context

Aitchison displacement shows the clearest treatment/quarter-condition signal. Cloxacillin-assigned quarters showed approximately **+20 Aitchison units** of excess displacement relative to untreated controls; the 20-quarter aggregated contrast remained significant after BH correction (`q = 0.019`). Bray–Curtis showed the same broad direction but no significant overall treatment effect.

Because treatment is fixed to anatomical quarter in the source experiment, treatment and quarter anatomy cannot be separated. These are therefore **treatment/quarter-condition associations**, not pure causal drug effects.

### Q2 — ecological decomposition

Cloxacillin-associated trajectories show a coherent multidimensional perturbation profile involving reduced core retention, increased turnover, and larger richness/dominance restructuring. Individual Q2 components do not survive full-family BH-FDR correction, so Q2 is interpreted primarily as ecological decomposition of the response phenotype.

### Q3 — baseline ecological susceptibility

Baseline richness, Shannon diversity, evenness, dominance, and multivariate community configuration predict later Bray–Curtis displacement in the discovery analyses. The portable Q3 model uses four baseline alpha-ecology features:

- richness
- Shannon diversity
- evenness
- dominance

The trajectory-level prediction target is:

`mean Bray perturbability = mean[d(T1,T2), d(T1,T3)]`

This gives **20 cow-quarter trajectories from 5 cows**.

Under nested **leave-one-cow-out** validation, the ecology-only ridge model achieved:

- RMSE: **0.10655**
- MAE: **0.08796**
- Pearson `r`: **0.5574**
- Spearman `rho`: **0.5895**
- CV `R²`: **0.2941**

### Treatment-adjusted Q3

We directly tested whether the Q3 predictive signal was simply capturing treatment identity. Treatment was learned only inside each nested LOCO training fold.

| Model | RMSE | MAE | Pearson r | Spearman rho | CV R² |
|---|---:|---:|---:|---:|---:|
| **baseline ecology only** | **0.10655** | **0.08796** | **0.5574** | **0.5895** | **0.2941** |
| baseline ecology + treatment | 0.11576 | 0.09824 | 0.4526 | 0.4226 | 0.1668 |
| treatment only | 0.14115 | 0.12554 | -0.2884 | -0.2977 | -0.2387 |

Ecology alone reduced LOCO RMSE by **24.5%** relative to treatment alone. Adding treatment to baseline ecology worsened RMSE by **8.6%**.

The Study 1 interpretation is therefore:

> **Treatment helps define the perturbation context; baseline ecology helps define susceptibility to perturbation.**

The Q3 signal is not simply a proxy for treatment assignment in this cohort. However, Study 1 contains only five cows, so this is discovery evidence rather than proof of universal transportability.

See [`docs/production_analysis_status_2026-08-12.md`](docs/production_analysis_status_2026-08-12.md) for the detailed analysis record.

## Study 1 perturbability model

Q3 is formalized as a portable ecological perturbability model in `scripts/07a_build_perturbability_model.R`.

Model design:

- ridge regression for correlated ecological predictors
- nested **leave-one-cow-out** cross-validation
- entire cows held out during validation
- cow-level bootstrap for coefficient uncertainty
- frozen model workflow for outcome-blinded external application

The companion script `scripts/07b_apply_perturbability_model.R` applies the frozen model to an external baseline table without using external follow-up outcomes.

## Cross-study replication phase

**Study 1 is now discovery-complete for the core Q1–Q3 hypothesis.** The next phase is replication across public longitudinal bovine milk microbiome datasets.

The cross-study framework will preserve the same conceptual questions while respecting each study's design:

`current/baseline ecological state + perturbation context + longitudinal history -> subsequent ecological displacement`

Priority datasets identified for replication include public bovine milk studies spanning dry-cow antimicrobial exposure, longitudinal lactation, farm/geographic variation, and shotgun/amplicon sequencing. The aim is to determine which Study 1 signals replicate across cows, farms, geography, sequencing protocols, and perturbation types rather than further optimizing the five-cow discovery cohort.

## Repository layout

```text
MMER/
├── config/
├── data/
├── docs/
├── manuscript/
├── results/
├── scripts/
├── tests/
├── workflow/
├── environment.yml
└── Makefile
```

## Quick start

```bash
conda env create -f environment.yml
conda activate mmer
python scripts/00_build_metadata.py
python tests/test_metadata.py
bash scripts/01_download_fastq.sh
```

## Data provenance

- Study: Biscarini F. et al. (2020), *A Randomized Controlled Trial of Teat-Sealant and Antibiotic Dry-Cow Treatments for Mastitis Prevention Shows Similar Effect on the Healthy Milk Microbiome.* Frontiers in Veterinary Science 7:581.
- DOI: `10.3389/fvets.2020.00581`
- ENA BioProject: `PRJEB38332`

## Reproducibility principles

- Preserve ENA accessions and original aliases.
- Treat quarter nested within cow as the biological unit for longitudinal displacement.
- Separate ecological resistance from antimicrobial resistance terminology.
- Treat treatment and anatomical quarter as confounded in Study 1.
- Treat baseline-predictor analyses as discovery/hypothesis-generating because Study 1 contains only five cows.
- Never evaluate model performance with quarter-level random train/test splitting; hold out cows as blocks.
- Learn treatment effects inside training folds when treatment is used predictively.
- Freeze external predictions before inspecting follow-up outcomes when performing prospective-style validation.
- Do not reinterpret ecological perturbability as disease risk without a separate disease-outcome test.

## Status

**Study 1 production pipeline:** complete.

**Q1:** complete.

**Q2:** complete.

**Q3:** complete, including nested LOCO and treatment-adjusted sensitivity analysis.

**Study 1 decision:** discovery-complete; avoid further model escalation on five cows.

**Next:** cross-study replication on public bovine milk microbiome cohorts.
