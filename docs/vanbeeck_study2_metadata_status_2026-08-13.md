# Van Beeck Study 2 metadata status — 2026-08-13

> **Historical metadata checkpoint.** Sequence processing and Q1–Q3 analyses have since been completed for the primary Study 2 milk cohort. See [`vanbeeck_study2_analysis_status_2026-08-14.md`](vanbeeck_study2_analysis_status_2026-08-14.md) for the current scientific status and [`MMER_hypothesis_and_replication_framework.md`](MMER_hypothesis_and_replication_framework.md) for the frozen discovery/replication framework.

## Dataset

Study 2 uses the public longitudinal bovine milk/teat-skin microbiome dataset deposited under ENA BioProject `PRJEB63336`.

The purpose of this cohort within MMER is cross-study replication of the Study 1 perturbability hypothesis using a larger multi-dairy longitudinal design.

## Run and BioSample inventory

Metadata were reconstructed from ENA run reports and EMBL-EBI BioSamples JSON records.

Current inventory at this checkpoint:

- 530 sequencing runs
- 514 unique BioSamples
- 208 unique longitudinal sample/core identifiers

The sample naming structure is `core.suffix`, for example `16922.1`, `16922.2`, `16922.3`.

## Longitudinal suffix mapping

Direct BioSamples metadata verification establishes:

- `.1` = `Baseline`
- `.2` = `7 Days`
- `.3` = `55-75DIM`

For example, core `16922` is deposited as:

- `16922.1`: Baseline
- `16922.2`: 7 Days
- `16922.3`: 55-75DIM

The suffix-derived timepoint mapping matched deposited timepoint metadata for all 530 run-level rows in the reconstructed manifest.

## Complete trajectories

Across the 208 core identifiers, timepoint coverage is:

- complete `1,2,3`: 114 cores
- `1,2` only: 36
- `1,3` only: 22
- `2,3` only: 20
- `1` only: 6
- `2` only: 4
- `3` only: 6

Thus there are **114 complete three-timepoint trajectories**.

## Sample type

Among complete trajectories:

- milk: **67**
- teat skin: **47**

The primary Study 2 MMER replication cohort begins from the **67 complete milk trajectories**.

## Dairy distribution

Complete milk trajectories:

- Dairy 1: 14
- Dairy 2: 21
- Dairy 3: 32

This multi-dairy structure provides the first MMER cohort suitable for direct testing of background/geographic transportability.

## Treatment/control structure

Complete milk trajectories are distributed as:

- `CB`: 22
- `CH`: 15
- `control high SCC`: 14
- `control low SCC`: 16

Dairy-by-treatment distribution:

| Dairy | CB | CH | control high SCC | control low SCC | Total |
|---|---:|---:|---:|---:|---:|
| Dairy 1 | 5 | 3 | 2 | 4 | 14 |
| Dairy 2 | 8 | 4 | 5 | 4 | 21 |
| Dairy 3 | 9 | 8 | 7 | 8 | 32 |
| **Total** | **22** | **15** | **14** | **16** | **67** |

## Additional deposited covariates

The BioSamples records expose additional longitudinal/contextual metadata including:

- dairy
- treatment
- sample type
- timepoint
- SCC
- SCC group
- MUN
- collection date

These fields are retained in the production manifest and can be used for descriptive or adjusted analyses where scientifically appropriate.

## Original planned MMER questions at this checkpoint

### Q1 — displacement and perturbation context

For each complete milk trajectory:

- `D12 = d(Baseline, 7 Days)`
- `D13 = d(Baseline, 55-75DIM)`
- `D23 = d(7 Days, 55-75DIM)` as an interval-specific secondary analysis

### Q3 — baseline ecological susceptibility

The main replication question was whether baseline ecological state predicts later displacement in this larger cohort.

Primary formulation:

`future displacement ~ baseline ecological state + treatment + dairy`

Prediction was planned using cow/core-blocked validation and train-on-two-dairies/test-on-the-third transportability analysis.

These analyses have now been executed; see the 2026-08-14 Study 2 analysis-status document for results.
