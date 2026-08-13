# Manitoba external validation protocol

## Goal

Evaluate whether the Study-1 MMER perturbability model generalizes to an independent longitudinal bovine milk microbiome cohort.

## Cohort

Derakhshani et al. (2018), Journal of Dairy Science, DOI 10.3168/jds.2018-14858.

The study includes 29 paired mammary secretion trajectories: milk collected at dry-off and corresponding colostrum collected immediately after calving from clinically healthy Holstein udder quarters. The quarters received penicillin G plus novobiocin and an internal teat sealant. Sequencing targeted the V1-V2 region of the bacterial 16S rRNA gene. Teat-canal samples were also collected, but the first MMER validation will use the paired mammary secretion samples.

## Primary question

Does the perturbability score learned in Study 1 from baseline mammary ecology correlate with the magnitude of later microbiome displacement in the independent Manitoba cohort?

## Analysis sequence

1. Build a complete sample/accession manifest and verify all 29 baseline/follow-up pairs.
2. Process Manitoba reads through a cohort-specific DADA2 workflow and perform QC without reference to external-validation outcomes.
3. Calculate the four prespecified baseline features from each pre-DCT milk sample: richness, Shannon diversity, evenness, and dominance.
4. Apply the frozen Study-1 model without refitting and generate one predicted perturbability score per baseline sample.
5. Freeze and archive the prediction table before calculating follow-up displacement.
6. Calculate observed dry-off-to-post-calving Bray-Curtis and Aitchison displacement for each paired trajectory.
7. Compare frozen predicted perturbability with observed displacement.

## Primary endpoint

Spearman correlation between frozen Study-1 perturbability prediction and observed Manitoba Bray-Curtis displacement.

Rank transportability is primary because the cohorts differ in 16S region, laboratory protocol, intervention, and sampling context.

## Secondary endpoints

- Pearson correlation for Bray-Curtis displacement.
- MAE and RMSE as descriptive calibration measures.
- Association between predicted perturbability and observed Aitchison displacement.
- Comparison of the prespecified Ridge model with the secondary Elastic Net benchmark.
- Sensitivity analyses based on predefined sequencing and sample QC.

## Rules to prevent leakage

- Do not retrain or tune model coefficients on Manitoba follow-up outcomes.
- Do not use post-calving information to construct baseline features.
- Finalize the baseline feature table before generating predictions.
- Freeze and archive predictions before calculating or merging observed follow-up displacement.
- Do not transfer Study-1 ASV-space PC1/PC2 because the cohorts use different amplicon regions.
- Do not use individual taxa as primary v1 predictors.

## Interpretation

Manitoba is not an external test of the Study-1 cloxacillin comparison because it does not contain the same untreated within-cow design. It is an external test of the Q3 concept: baseline ecological state predicts future ecological displacement.

The physiological dry-period transition is part of the observed Manitoba change. Therefore the validation concerns general mammary microbiome perturbability under a distinct longitudinal context, not an isolated antibiotic effect.

## Planned outputs

- manitoba_sample_manifest.csv
- manitoba_read_tracking.csv
- manitoba_baseline_features.csv
- manitoba_predictions_FROZEN.csv
- manitoba_observed_displacement.csv
- manitoba_external_validation_metrics.csv
- prediction-versus-observed figure

## Decision rule

External validation will be judged primarily from effect size, direction, uncertainty, and consistency with the Study-1 prediction rather than from a single significance threshold. A reproducible positive rank association would support transportability of the perturbability phenotype; weak or absent association would indicate limited transportability and will be reported accordingly.
