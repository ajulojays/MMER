# Van Beeck external validation of mammary ecological resistance

**Frozen MMER external-validation workflow — 2026-08-16**

## Executive conclusion

The Van Beeck cohort provides a second independent 16S external validation of the MMER hypothesis that **baseline multivariate mammary microbial architecture is associated with subsequent ecological resistance**.

The frozen Van Beeck external-validation analysis uses all **60 complete longitudinal milk trajectories** and defines ecological resistance relative to each trajectory's own baseline. Baseline family-level ecological architecture is learned independently within Van Beeck and tested using the prespecified joint PC1 + PC2 hypothesis while adjusting for treatment group and dairy.

Primary model:

`overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy`

Primary result:

- N = **60** complete trajectories
- full-model R² = **0.316**
- adjusted R² = **0.224**
- partial R² for PC1 + PC2 = **0.216**
- classical partial F = **7.18**, df = 2,52, p = **0.00177**
- HC3 joint Wald χ² = **9.13**, df = 2, p = **0.0104**
- 9,999-stratum Freedman–Lane permutation p = **0.0183**
- PC1 beta = **-0.0172**, HC3 p = 0.283
- PC2 beta = **-0.0665**, HC3 p = **0.0246**
- PC1 95% stratified-bootstrap interval = **-0.0439 to +0.00665**
- PC2 95% stratified-bootstrap interval = **-0.104 to -0.0217**
- PC2 was negative in **99.6%** of 5,000 bootstrap replicates
- partial-R² bootstrap median = **0.223**; percentile interval = **0.0499 to 0.566**

The external-validation conclusion is based on the **joint PC1 + PC2 hypothesis**, not on PC2 alone.

> **Frozen interpretation:** baseline family-level community architecture significantly predicted subsequent ecological resistance in the independent Van Beeck cohort after accounting for treatment and dairy, providing independent phenomenon-level validation of the ecological-resistance architecture identified in Italy + Manitoba.

---

# 1. Background and rationale

## 1.1 Ecological resistance as an individual longitudinal phenotype

Microbiome studies commonly ask whether average community composition differs between treatments, farms, physiological stages, or sampling times. Those analyses are important, but they do not directly answer why biological units exposed to broadly similar transitions can differ substantially in the magnitude of community change.

MMER reframes longitudinal bovine mammary microbiome data around **ecological resistance**. For biological unit *i*:

`R_i = 1 - Bray-Curtis(Baseline_i, Follow-up_i)`

Higher resistance means that the follow-up community retains more of the unit's own baseline compositional structure. Lower resistance means greater ecological displacement.

The core MMER question is therefore not simply whether a perturbation changes the microbiome on average. It is whether **the initial ecological organization of the community contains information about the magnitude of its subsequent displacement**.

This makes between-unit heterogeneity the phenotype rather than treating it only as residual variation around a treatment or time effect.

## 1.2 Discovery in Italy + Manitoba

The frozen discovery/development analysis combines healthy mammary-quarter trajectories from Italy and Manitoba. Family-level baseline compositions are CLR transformed, centered within study, and represented by pooled principal components. The joint PC1 + PC2 term predicts subsequent ecological resistance:

- 49 quarter trajectories from 14 cows
- joint PC1 + PC2 HTZ F = **10.0**
- df = 2,7.81
- p = **0.00699**
- descriptive R² = **0.334**

That analysis establishes the discovery hypothesis:

> **baseline multivariate mammary-community architecture is associated with subsequent ecological resistance.**

The next scientific requirement is not to repeatedly optimize the discovery datasets. It is to challenge the hypothesis in independent longitudinal cohorts with different farms, animals, treatment structures, sampling schedules, and laboratory histories.

## 1.3 Why Van Beeck is an informative external cohort

Van Beeck et al. studied milk and teat-skin microbiota across the dry-off/lactation transition at **three California dairies**, with sampling at baseline dry-off, 7 days later, and 55–75 days in milk in the next lactation. The source design included low-SCC untreated controls, high-SCC untreated controls, and high-SCC groups receiving intramammary cephapirin benzathine (CB) or ceftiofur hydrochloride (CH). The raw sequencing study is available under **PRJEB63336**; the source publication is Van Beeck et al., *Applied and Environmental Microbiology* (2026), DOI **10.1128/aem.02312-25**.

The source paper primarily evaluated microbiota variation associated with dairy, sampling time, SCC context, and intramammary cephalosporin use. MMER uses the longitudinal structure for a different purpose: each retained trajectory becomes an ecological system with a baseline state and a subsequently observed degree of compositional retention.

Van Beeck is particularly valuable for external validation because it differs from the discovery setting in several ways:

- three dairies rather than the Italy/Manitoba study structure;
- a cow-level selected-quarter longitudinal sampling framework;
- four SCC/treatment contexts;
- two post-baseline observations rather than a single follow-up endpoint;
- an independent sequencing experiment and laboratory history;
- a different geographic and production context.

These differences make Van Beeck a genuine test of **generalization**, not a technical replicate.

---

# 2. How the MMER question differs from the original Italy and Manitoba questions

The following distinctions are central to the novelty of MMER and should remain explicit in manuscripts, presentations, and reviewer responses.

## 2.1 From average treatment/time effects to response heterogeneity

The Italy and Manitoba source studies primarily characterized how mammary microbial communities differed across interventions, the dry period, and sampling times.

Their broad direction of inference was:

`Treatment / time -> microbiome change`

MMER instead asks:

`Baseline ecological state -> magnitude of later ecological displacement`

The object of interest is therefore **heterogeneity in response among biological units**, not only the average response of a treatment or time group.

## 2.2 A new phenotype is created from the longitudinal data

The original studies analyzed diversity, composition, temporal change, persistence, or taxonomic shifts as outcomes.

MMER constructs a new unit-level phenotype:

`ecological_resistance = 1 - Bray-Curtis(own baseline, own follow-up)`

This converts a pair of longitudinal microbiome samples into a quantitative measure of ecological retention for each trajectory.

## 2.3 Baseline community architecture becomes the predictor

The original studies did not primarily ask whether the *pre-transition multivariate community state* predicted the magnitude of subsequent change.

MMER explicitly makes baseline community organization the explanatory variable. Baseline family composition is represented using CLR transformation and multivariate ordination, and the joint baseline architecture is tested against later resistance.

## 2.4 MMER does not reduce baseline ecology to alpha diversity

Richness, Shannon diversity, and other alpha-diversity summaries cannot capture the full identity and abundance configuration of the baseline community.

MMER therefore uses **multivariate family-level architecture**. Two samples can have comparable diversity yet occupy very different compositional states; the ecological-resistance hypothesis is about this broader state.

## 2.5 The longitudinal biological unit is preserved

The key linkage is always the same biological unit through time:

`Baseline_i -> Follow-up_i`

In Italy/Manitoba this is a quarter trajectory; in Van Beeck it is the retained longitudinal milk trajectory associated with the selected quarter/cow sampling unit. MMER does not replace that linkage with group-average distances.

## 2.6 Between-unit variability is signal, not nuisance

Conventional longitudinal analyses often summarize average changes and treat remaining individual variability as unexplained error.

MMER asks whether part of that variability is **structured by initial ecological state**. Thus, heterogeneity around the average trajectory is the biological phenomenon under study.

## 2.7 MMER explicitly separates discovery from external validation

Italy + Manitoba are retained as discovery/development cohorts. Van Beeck is not used to redefine the discovery model, select more favorable PCs, or search for a new primary endpoint.

Instead, Van Beeck receives an independently constructed baseline ecological representation and the same phenomenon-level hypothesis is tested externally.

This discovery/validation separation was not a requirement of the original source studies because they were answering their own within-study biological questions.

## 2.8 MMER solves a cross-study harmonization problem

The original studies did not need to create a common ecological framework across independent cohorts.

MMER harmonizes the inference at the level of:

- bacterial family composition;
- CLR-transformed baseline states;
- a baseline-relative Bray-Curtis resistance phenotype;
- a joint two-axis multivariate architecture test.

The cohorts can differ in exact experimental details while still testing the same ecological phenomenon.

## 2.9 MMER tests initial-state dependence rather than only temporal difference

The temporal logic changes from:

`How different are communities at T1 and T2?`

into:

`Does M(T1) predict d[M(T1), M(T2)]?`

This is an initial-state-dependence question closely related to ecological stability, resistance, susceptibility, and context dependence.

## 2.10 MMER is complementary to, not a correction of, the source studies

The original Italy and Manitoba analyses were appropriate for their stated objectives. MMER does not claim those analyses were incomplete or wrong.

A concise synthesis is:

> **The source studies describe average longitudinal and intervention-associated trajectories; MMER asks why individual mammary ecological systems differ around those average trajectories.**

---

# 3. Frozen Van Beeck external-validation question

## Primary question

> **Does independently learned baseline family-level mammary microbial architecture predict subsequent ecological resistance in the Van Beeck cohort after accounting for the cohort's treatment and dairy structure?**

## Primary null hypothesis

`H0: beta_PC1 = beta_PC2 = 0`

## Primary alternative

At least one component of the two-dimensional baseline ecological architecture contributes information about overall ecological resistance beyond treatment and dairy.

## What this analysis is not testing

The canonical external-validation analysis does **not** test:

- whether CB differs from CH;
- whether treated animals differ from untreated animals in the architecture slope;
- whether high-SCC animals differ from low-SCC animals in the architecture slope;
- whether 7-day resistance is more predictable than 55–75-DIM resistance;
- whether a particular family is a causal biomarker;
- whether the sign of Van Beeck PC1 or PC2 matches the sign of discovery PCs;
- whether the discovery loading vectors transport exactly into Van Beeck.

Those are separate hypotheses. They must not be allowed to redefine the external-validation claim.

---

# 4. Methods

## 4.1 Source cohort and longitudinal design

The Van Beeck source study sampled dairy cows at three California dairies at:

1. **Baseline** at dry-off;
2. **7 Days** after baseline;
3. **55–75 DIM** in the subsequent lactation.

The source design contained four biological/treatment groups:

- low-SCC untreated control;
- high-SCC untreated control;
- cephapirin benzathine (CB);
- ceftiofur hydrochloride (CH).

The MMER locked sequencing-depth workflow retains **60 complete milk trajectories**, giving **180 samples** across the three timepoints.

Locked trajectory counts are:

| Group | Trajectories |
|---|---:|
| control low SCC | 14 |
| control high SCC | 12 |
| CB | 18 |
| CH | 16 |
| **Total** | **60** |

These groups are retained in the primary model as design covariates rather than converted into subgroup validation analyses.

## 4.2 Sequencing matrix

The canonical workflow starts from the locked rarefied bacterial ASV table:

`results/vanbeeck/analysis_3k/seqtab_bacterial_rarefied_3000.rds`

and the corresponding bacterial taxonomy:

`results/vanbeeck/taxonomy/taxonomy_bacterial.rds`

The metadata source is explicitly the locked file:

`results/vanbeeck/analysis_3k/metadata_60_complete_trajectories_locked.csv`

The script requires:

- 180 samples;
- 60 unique trajectories;
- exactly three samples per trajectory;
- one Baseline, one 7 Days, and one 55-75DIM observation per trajectory;
- the locked treatment totals above.

This prevents the previously observed unlocked-metadata duplication from entering the external-validation analysis.

## 4.3 Removal of non-bacterial organellar features

ASVs with taxonomy strings matching mitochondrial, chloroplast, plastid, or organelle annotations are removed before family aggregation.

## 4.4 Family aggregation

ASVs are aggregated to bacterial family. When Family is missing, conservative fallback labels are constructed using the deepest available higher rank:

1. Order;
2. Class;
3. Phylum;
4. `Unclassified_Bacteria`.

The family count matrix is then converted to relative abundance by sample closure.

## 4.5 Ecological-resistance phenotype

For trajectory *i*:

`R_early,i = 1 - Bray-Curtis(Baseline_i, 7Days_i)`

`R_late,i = 1 - Bray-Curtis(Baseline_i, 55-75DIM_i)`

The frozen primary Van Beeck phenotype is:

`R_overall,i = mean(R_early,i, R_late,i)`

Equivalently:

`R_overall,i = 1 - mean[d(Baseline_i,7Days_i), d(Baseline_i,55-75DIM_i)]`

This choice uses both post-baseline observations to summarize the trajectory's overall retention relative to baseline without selecting a post hoc time window according to significance.

The component distances are retained as transparent intermediate outputs but are **not separately hypothesis-tested in the external-validation-only script**.

## 4.6 Baseline prevalence filter

Only baseline samples are used to define the predictor feature space.

A family is retained when detected in at least **10% of the 60 baseline trajectories**.

This filter is applied before CLR transformation and PCA.

## 4.7 Zero replacement and CLR transformation

For the retained baseline relative-abundance matrix, the pseudocount is data-derived:

`pseudo = 0.5 × minimum positive baseline relative abundance`

Only zero entries are replaced. Samples are re-closed after zero replacement and transformed using the centered log ratio:

`CLR(x_ij) = log(x_ij) - mean_j[log(x_ij)]`

A large arbitrary pseudocount such as 0.5 on the relative-abundance scale is not used.

## 4.8 Within-dairy centering

Because Van Beeck includes three dairies and the source paper itself found substantial location-associated microbiome variation, CLR features are centered within dairy before PCA.

For each family feature, the dairy-specific CLR centroid is subtracted from samples in that dairy.

This reduces the risk that the baseline principal axes merely rediscover dairy identity.

**Treatment is not centered out before PCA.** Treatment/SCC context remains part of the observed ecological state, and treatment is handled as a design covariate in the outcome model.

## 4.9 Independent Van Beeck PCA

PCA is fitted independently on the within-dairy-centered Van Beeck baseline CLR matrix.

PC1 and PC2 scores are standardized before regression.

This analysis therefore tests **phenomenon-level external validation**:

> does multivariate baseline architecture predict resistance in another cohort?

It does not test exact transportability of the Italy/Manitoba PC loading vectors.

Because independent PCA has arbitrary sign and can rotate when feature covariance differs, the signs of PC coefficients must **not** be compared biologically across cohorts as if PC1 or PC2 were the same fixed axis.

Exact axis transport would require a different analysis with harmonized features, frozen discovery centering/loading vectors, and projection of external samples into that fixed discovery coordinate system.

## 4.10 Primary regression model

The reduced design-only model is:

`overall_resistance ~ treatment + dairy`

The full external-validation model is:

`overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy`

The inferential target is the incremental contribution of PC1 + PC2 beyond treatment and dairy.

There is one outcome observation per retained Van Beeck trajectory. Therefore, cow/core-clustered CR2 with one-observation clusters is not used.

## 4.11 Classical joint test and effect size

The nested models are compared with a 2-df partial F-test.

The architecture effect size is reported as partial R²:

`partial_R2 = (RSS_reduced - RSS_full) / RSS_reduced`

This quantifies the fraction of residual variation from the treatment+dairy model that is accounted for by adding the two baseline ecological axes.

## 4.12 HC3 robust inference

HC3 heteroskedasticity-consistent covariance estimates are calculated for the full model. The joint PC1 + PC2 hypothesis is evaluated using a 2-df robust Wald statistic.

Individual PC coefficient p-values are reported descriptively, but the prespecified validation decision is based on the **joint two-axis architecture test**.

## 4.13 Freedman–Lane permutation

A design-aware residual permutation is used as an additional robustness test.

1. Fit the reduced model containing treatment and dairy.
2. Preserve the reduced-model fitted values.
3. Permute reduced-model residuals within **dairy × treatment** strata.
4. Reconstruct the permuted outcome.
5. Refit reduced and full models.
6. Record the partial F statistic.
7. Repeat **9,999** times.

The Monte Carlo permutation p-value is:

`(1 + number[F_perm >= F_observed]) / (9999 + 1)`

This preserves the major Van Beeck design structure while evaluating whether the baseline architecture adds more explanatory information than expected under the null.

## 4.14 Stratified bootstrap

The workflow performs **5,000 complete-trajectory bootstrap replicates** within dairy × treatment strata.

The Van Beeck PCA coordinate system is held fixed during this bootstrap. Thus, the bootstrap quantifies uncertainty in the regression coefficients and partial R² **conditional on the independently learned Van Beeck ecological coordinate system**.

Reported bootstrap quantities include:

- PC1 coefficient percentile interval;
- PC2 coefficient percentile interval;
- coefficient sign stability;
- partial-R² bootstrap distribution.

The partial-R² interval is an effect-size uncertainty summary and is **not** interpreted as a separate null-hypothesis test, because adding predictors to nested OLS models constrains in-sample partial R² to be non-negative.

## 4.15 Reproducibility

Canonical analysis script:

`scripts/26j_vanbeeck_external_validation_only.R`

Canonical figure script:

`scripts/26k_vanbeeck_external_validation_plots.R`

Primary output directory:

`results/vanbeeck/external_validation_only/`

Figure output directory:

`figures/vanbeeck_external_validation/`

---

# 5. Results

## 5.1 Cohort lock

The external-validation cohort contains **60 complete longitudinal trajectories / 180 samples** distributed across the four locked groups:

- low-SCC untreated control: 14;
- high-SCC untreated control: 12;
- CB: 18;
- CH: 16.

All trajectories contribute one Baseline, one 7 Days, and one 55–75 DIM sample.

## 5.2 Primary joint architecture test

The full model explained **31.6%** of observed variation in overall ecological resistance (`R² = 0.316`; adjusted `R² = 0.224`).

Adding baseline PC1 + PC2 to the treatment+dairy model accounted for **21.6% of the variation remaining after the design covariates** (`partial R² = 0.216`).

The joint architecture hypothesis was supported by all three inferential routes:

| Test | Statistic | p |
|---|---:|---:|
| Classical partial F | F(2,52) = 7.18 | **0.00177** |
| HC3 joint Wald | χ²(2) = 9.13 | **0.0104** |
| Freedman–Lane permutation | 9,999 permutations | **0.0183** |

The concordance of the classical, heteroskedasticity-robust, and design-aware permutation analyses is the primary evidence for external validation.

## 5.3 Individual baseline axes

PC1:

- beta = **-0.0172**;
- HC3 SE = 0.0158;
- HC3 p = 0.283.

PC2:

- beta = **-0.0665**;
- HC3 SE = 0.0287;
- HC3 p = **0.0246**.

Although PC2 contributes more visibly to the observed association, the primary hypothesis remains joint PC1 + PC2. PC2 should not be retroactively promoted to a single-axis primary endpoint.

## 5.4 Bootstrap stability

All **5,000 / 5,000** requested bootstrap replicates produced valid estimates.

PC1:

- observed beta = -0.0172;
- bootstrap mean = -0.0192;
- bootstrap median = -0.0192;
- 95% percentile interval = **-0.0439 to +0.00665**;
- 93.2% of bootstrap coefficients were negative.

PC2:

- observed beta = -0.0665;
- bootstrap mean = -0.0663;
- bootstrap median = -0.0669;
- 95% percentile interval = **-0.104 to -0.0217**;
- 99.6% of bootstrap coefficients were negative.

Joint architecture effect size:

- observed partial R² = **0.216**;
- bootstrap mean = 0.251;
- bootstrap median = **0.223**;
- 95% percentile interval = **0.0499 to 0.566**.

The bootstrap therefore supports a stable negative Van Beeck PC2 association and a non-trivial joint baseline-architecture effect size, while also showing substantial uncertainty in the exact magnitude of the explained variation.

---

# 6. Discussion

## 6.1 Van Beeck independently validates the ecological-resistance phenomenon

The most important result is not the sign of PC2. It is that a baseline multivariate ecological representation learned independently in Van Beeck significantly predicts later ecological resistance after accounting for the study's dairy and treatment structure.

This is precisely the level of evidence required for **phenomenon-level replication** of the Italy + Manitoba discovery result.

The external cohort differs in geography, dairies, animals, experimental structure, treatment context, and longitudinal schedule. A significant joint architecture association under those differences argues that the discovery signal is not simply a numerical peculiarity of the pooled Italy/Manitoba matrix.

## 6.2 Treatment and dairy are adjustment variables, not the new scientific question

Van Beeck has a multifactorial source design, and the original publication found important microbiota variation associated with sampling location and time. MMER therefore retains dairy and treatment as design covariates.

The external-validation claim is:

> baseline architecture predicts resistance **beyond** those known structural factors.

The canonical analysis intentionally excludes treatment-interaction and SCC-subgroup hypothesis testing. Those questions may be scientifically interesting, but they are not needed to establish external validation and could distract from the predefined replication objective.

## 6.3 The negative PC2 coefficient is cohort-specific

Van Beeck PC2 is strongly and stably negatively associated with resistance. However, an independently fitted PCA does not create the same coordinate system as the Italy/Manitoba PCA or the Porcellato PCA.

PCA signs are arbitrary, and axes can rotate or exchange according to cohort covariance structure. Therefore:

- Van Beeck PC2 negative does not contradict a positive discovery PC2 coefficient;
- signs should not be compared as if they represented the same ecological gradient;
- families defining Van Beeck PC2 are useful for describing the local ecological state, not for claiming a universal PC2 biomarker signature.

The replicated object is the **relationship between baseline multivariate architecture and later resistance**, not a fixed numerical axis.

## 6.4 Relation to Porcellato and Patangia

The emerging MMER evidence structure now includes two independent positive 16S validations:

- Porcellato: strong positive phenomenon-level validation;
- Van Beeck: positive phenomenon-level validation with robust and permutation support.

Patangia's prespecified shotgun taxonomic test remains null.

That pattern is informative. It supports generalizability without implying universal detectability in every cohort or platform. The current data do not justify claiming that 16S inherently captures the phenomenon while shotgun metagenomics does not; platform, sample size, biological window, feature construction, host/background DNA, and cohort context are confounded.

## 6.5 Effect size is biologically meaningful but uncertain

A partial R² of 0.216 indicates that the two baseline ecological axes explain a substantial fraction of the variation left after treatment and dairy adjustment. At the same time, the bootstrap interval is broad, reflecting finite-sample uncertainty.

The correct interpretation is therefore neither “the baseline microbiome determines resistance” nor “21.6% is a universal constant.” Instead, Van Beeck provides evidence that baseline architecture contains appreciable information about future ecological retention in this independent cohort.

## 6.6 Association is not causation

The source data were not prospectively collected to manipulate baseline community architecture. MMER therefore identifies **initial-state dependence consistent with ecological resistance**, not a causal mechanism.

Potential explanations may involve taxonomic interactions, host state, milk environment, immune context, prior exposures, or latent variables correlated with baseline composition. These cannot be separated by the present retrospective reanalysis.

## 6.7 Why the null hypothesis matters

External validation is strongest when the analysis is capable of failing. Van Beeck was not allowed to redefine the hypothesis through subgroup selection, post hoc time-window selection, or promotion of whichever PC produced the smallest p-value.

The locked all-60 analysis tests a single clear statement. The positive result is therefore more informative than a collection of optimized secondary analyses.

---

# 7. Publication-ready result language

## Results paragraph

> In the independent Van Beeck cohort, baseline family-level mammary microbial architecture significantly predicted subsequent ecological resistance after adjustment for treatment and dairy (N = 60). Addition of the first two baseline ecological principal components explained 21.6% of the residual variation beyond these design covariates (partial R² = 0.216; F2,52 = 7.18, P = 0.00177). The joint association remained supported using heteroskedasticity-robust inference (HC3 Wald χ²2 = 9.13, P = 0.0104) and 9,999 Freedman–Lane permutations within dairy-by-treatment strata (P = 0.0183). PC2 showed the stronger individual association (β = -0.0665, HC3 P = 0.0246), with a 95% stratified-bootstrap interval of -0.104 to -0.0217 and a negative coefficient in 99.6% of 5,000 bootstrap replicates. Because the Van Beeck PCA was fitted independently, these axis directions are cohort-specific; the result constitutes phenomenon-level validation rather than transport of identical discovery PC axes.

## Discussion paragraph

> Van Beeck provides a second independent 16S validation of the MMER ecological-resistance hypothesis. Despite differences in geography, dairy, treatment structure, and longitudinal sampling context, baseline multivariate community organization retained predictive information about subsequent compositional retention. The convergence of classical, HC3-robust, permutation, and bootstrap analyses argues against the association being driven solely by a particular variance assumption. Importantly, the replicated object is not the sign or taxonomic loading pattern of a specific principal component, because PCA axes were learned independently in each cohort. Rather, the repeated finding is that individual mammary microbial systems occupy baseline compositional states that are associated with how strongly they retain their baseline organization through a subsequent biological transition.

---

# 8. Figures

The canonical plotting workflow is:

`scripts/26k_vanbeeck_external_validation_plots.R`

It produces the following publication-ready PNG and PDF figures from the frozen analysis outputs.

## Figure VB1 — distribution of overall ecological resistance

`Fig_VB1_overall_resistance_distribution`

Purpose: show the individual-level phenotype and heterogeneity across all 60 trajectories without subgroup inference.

## Figure VB2A — adjusted PC1 association

`Fig_VB2A_PC1_added_variable`

Purpose: added-variable visualization of the PC1 coefficient after adjustment for PC2, treatment, and dairy.

## Figure VB2B — adjusted PC2 association

`Fig_VB2B_PC2_added_variable`

Purpose: visualize the stronger Van Beeck PC2 contribution without redefining PC2 as the primary hypothesis.

## Figure VB3 — bootstrap coefficient intervals

`Fig_VB3_bootstrap_PC_coefficients`

Purpose: show observed PC1/PC2 coefficients with 95% stratified-bootstrap intervals.

## Figure VB4 — bootstrap partial-R² distribution

`Fig_VB4_bootstrap_partial_R2`

Purpose: show uncertainty in the magnitude of the joint architecture effect. This is an effect-size plot, not a bootstrap significance test.

## Figure VB5 — joint-test concordance

`Fig_VB5_joint_test_concordance`

Purpose: display the agreement among classical partial F, HC3 joint Wald, and Freedman–Lane permutation inference.

## Figure VB6 — baseline PCA space

`Fig_VB6_baseline_PCA`

Purpose: display the independently learned Van Beeck ecological coordinate system as ordination/QC context. The caption explicitly warns against cross-cohort sign equivalence.

Three lightweight vector summary figures are also committed under `figures/vanbeeck_external_validation/` for immediate GitHub rendering.

---

# 9. Output contract

The canonical analysis writes:

- `locked_60_trajectory_design.csv`
- `vanbeeck_family_ecological_resistance_60.csv`
- `vanbeeck_primary_resistance_summary.csv`
- `vanbeeck_baseline_family_prevalence.csv`
- `vanbeeck_PCA_variance_explained.csv`
- `vanbeeck_baseline_PCA_scores.csv`
- `vanbeeck_baseline_PCA_loadings.csv`
- `vanbeeck_PC1_top25_family_loadings.csv`
- `vanbeeck_PC2_top25_family_loadings.csv`
- `vanbeeck_external_validation_analysis_table.csv`
- `PRIMARY_HC3_coefficients.csv`
- `PRIMARY_external_validation_joint_PC1_PC2.csv`
- `PERMUTATION_9999_partial_F.csv`
- `BOOTSTRAP_5000_replicates.csv`
- `BOOTSTRAP_5000_summary.csv`
- `VANBEECK_EXTERNAL_VALIDATION_FREEZE.csv`
- `PRIMARY_OLS_model_summary.txt`
- `sessionInfo.txt`

---

# 10. Guardrails

1. Van Beeck is an **external-validation cohort**, not a discovery cohort.
2. The primary cohort is all **60 locked complete trajectories**.
3. The primary taxonomic level is **family**.
4. The primary outcome is **overall resistance**, defined before inferential testing as the mean of baseline-relative resistance at the two post-baseline timepoints.
5. The primary model is `overall_resistance ~ z_PC1 + z_PC2 + treatment + dairy`.
6. The primary hypothesis is the **joint PC1 + PC2** term.
7. Treatment and dairy are design covariates; they are not the external-validation scientific target.
8. No subgroup or treatment-interaction result can replace the all-60 primary test.
9. HC3 and stratified permutation are the principal robustness checks because there is one observation per trajectory in the endpoint model.
10. The bootstrap quantifies uncertainty conditional on frozen Van Beeck PCA coordinates.
11. Independent PCA supports phenomenon-level replication only; exact PC signs or loading identities are not transported across studies.
12. Individual loading families are axis-defining features, not proven causal biomarkers.
13. Patangia's primary null remains part of the evidence base and must not be hidden because Van Beeck is positive.
14. A positive external validation is evidence of reproducibility across contexts, not proof of universal ecological law or causality.

---

# 11. Source reference

Van Beeck W, Lemos MLP, Niesen AM, Finnegan P, Shih TM, Ho A, Rossow HA, Marco ML. **Variations in cow milk and teat skin microbiota across the lactation cycle with intramammary cephalosporin use at dry-off.** *Applied and Environmental Microbiology*. 2026;92(6):e02312-25. DOI: **10.1128/aem.02312-25**. ENA project: **PRJEB63336**.
