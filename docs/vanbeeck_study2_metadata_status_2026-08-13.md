# Van Beeck Study 2 metadata status — 2026-08-13

## Dataset

Study 2 uses the public longitudinal bovine milk/teat-skin microbiome dataset deposited under ENA BioProject `PRJEB63336`.

The purpose of this cohort within MMER is cross-study replication of the Study 1 perturbability hypothesis using a larger multi-dairy longitudinal design.

## Run and BioSample inventory

Metadata were reconstructed from ENA run reports and EMBL-EBI BioSamples JSON records.

Current inventory:

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

The primary Study 2 MMER replication cohort is the **67 complete milk trajectories**.

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

These fields will be retained in the production manifest and used for descriptive or adjusted analyses where scientifically appropriate.

## Planned MMER Study 2 questions

### Q1 — displacement and perturbation context

For each complete milk trajectory:

- `D12 = d(Baseline, 7 Days)`
- `D13 = d(Baseline, 55-75DIM)`
- `D23 = d(7 Days, 55-75DIM)` as an interval-specific secondary analysis

Primary Q1 tests will evaluate whether displacement differs by treatment/control stratum and dairy.

### Q3 — baseline ecological susceptibility

The main replication question is whether baseline ecological state predicts later displacement in this larger cohort.

Primary formulation:

`future displacement ~ baseline ecological state + treatment + dairy`

Prediction will use cow/core-blocked validation. Because Study 2 includes three dairies, an additional transportability analysis is planned:

`train on two dairies -> test on the third`

This directly tests the working hypothesis that baseline microbiome background varies by location and can influence ecological trajectory or perturbability.

## Cross-study objective

Study 1 suggested that treatment helps define perturbation context while baseline ecology helps define susceptibility. Study 2 will test whether that relationship survives:

- a larger number of longitudinal trajectories;
- multiple dairies;
- different treatment/control strata;
- a different cohort and sequencing experiment.

The emphasis is replication and transportability rather than further optimization of the Study 1 model.

## Immediate next step

Construct the milk-only FASTQ run manifest for the 67 complete milk trajectories, corresponding to 201 biological timepoint samples plus any attached technical replicate runs, then begin production sequence processing.
