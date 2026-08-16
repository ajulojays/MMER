# How MMER Fundamentally Differs from the Two Source Studies

**Critical conceptual and methodological comparison**

This document summarizes the key distinctions between the original Italy and Manitoba mammary microbiome studies and the frozen MMER two-aim reanalysis. The original studies focused primarily on treatment- and time-associated microbiome change. MMER instead asks whether the baseline ecological state of an individual healthy mammary quarter is associated with that same quarter's subsequent ecological resistance across the dry-off to early-postpartum transition.

## 1. The Italy study asked an intervention-comparison question

The Italy study was primarily designed to compare dry-cow interventions in healthy mammary quarters. Its central inferential structure was approximately:

`Treatment / time -> microbiome composition`

The study compared diversity and community composition across untreated, teat-sealant, cephalonium, and cloxacillin-assigned quarters and across dry-off, calving, and early-lactation timepoints.

MMER does not use treatment response as its primary question. Instead, it treats between-quarter variation in longitudinal response as the phenotype of interest.

A further design consideration is that intervention was tied to anatomical quarter position in the deposited design. Consequently, treatment and quarter position cannot be cleanly separated in our reanalysis, which is one reason treatment modification is not part of the frozen primary MMER framework.

## 2. The Manitoba study asked a dry-period community-dynamics question

The Manitoba study focused on how teat-canal and intramammary microbial communities changed through the dry period and after antimicrobial dry-cow therapy plus internal teat sealant.

Its broad inferential structure was:

`Dry period + DCT -> community restructuring / persistence`

The paper emphasized temporal change, taxonomic shifts, richness changes, and persistence of organisms across the dry period.

MMER uses the same longitudinal information for a different purpose: to quantify each quarter's ecological response and ask whether that response is related to the quarter's baseline community state.

## 3. MMER changes the biological question

The original studies ask, in essence:

> What happens to the mammary microbiome after treatment or across time?

MMER asks:

> Why do individual healthy mammary quarters differ in how much their microbiomes change?

The central conceptual shift is therefore:

`Original studies: perturbation -> community change`

versus

`MMER: baseline community state -> susceptibility / resistance to perturbation`

The arrow in the MMER formulation represents temporal association, not demonstrated causality.

## 4. MMER elevates between-quarter heterogeneity into the phenotype

Group-level averages can conceal large differences among individual mammary quarters exposed to the same broad physiological transition.

MMER defines quarter-level ecological resistance as:

`R(iq) = 1 - Bray-Curtis[T1(iq), T2(iq)]`

Higher values indicate greater compositional retention; lower values indicate greater displacement.

Rather than treating between-quarter variability as residual noise around an average treatment or time effect, MMER makes that variability the central biological object of study.

## 5. MMER treats the mammary quarter as a longitudinal ecological system

The biological unit is the quarter trajectory, nested within cow.

Each quarter is treated as a small ecological system with:

- an initial community configuration at dry-off;
- exposure to a major dry-period / postpartum transition;
- a subsequent degree of ecological displacement.

The key linkage is same-quarter baseline to same-quarter outcome:

`T1_quarter -> T2_same quarter`

This transforms longitudinal sampling into an ecological susceptibility/resistance framework rather than a conventional cross-sectional differential-abundance comparison.

## 6. MMER does not reduce baseline ecology to diversity alone

The original studies appropriately examined familiar alpha-diversity and beta-diversity measures.

The frozen MMER analysis instead represents baseline ecological state using multivariate family-level community architecture:

`T1 family composition -> prevalence filtering -> CLR transformation -> within-study centering -> pooled PCA`

This matters because two communities can have similar richness or Shannon diversity while containing very different organisms and abundance structures.

In the frozen analysis, the joint PC1 + PC2 baseline architecture term is associated with subsequent ecological resistance:

`F(2, 7.81) = 10.0, P = 0.00699`

PC2 is the strongest pooled axis:

`beta = +0.0622, P = 0.00302`

## 7. MMER harmonizes independent studies into one ecological framework

The Italy and Manitoba datasets were not designed as a combined experiment. They differ in geography, exact sampling schedule, intervention context, and 16S variable region.

MMER constructs a common ecological representation by harmonizing:

- bacterial family-level profiles;
- a comparable dry-off / pre-intervention baseline;
- an early-postpartum ecological state;
- the same quarter-level Bray-Curtis resistance phenotype.

This cross-study harmonization is itself a major methodological contribution because it allows the same ecological hypothesis to be evaluated across independent cohorts despite differences in the original study designs.

## 8. MMER removes study centroids before learning the ecological axes

A naive pooled PCA would largely rediscover cohort identity because geography, farm, sequencing region, laboratory procedures, and biology all contribute to between-study compositional separation.

MMER therefore subtracts the study-specific CLR centroid before pooled PCA.

The resulting PCs are intended to represent within-study ecological variation rather than a simple Italy-versus-Manitoba axis.

This is a methodological problem the original papers did not need to solve because neither study was designed for cross-cohort ecological-state inference.

## 9. MMER interprets taxa as components of a multivariate state, not isolated biomarkers

The original papers appropriately reported taxa that increased, decreased, or persisted through time.

MMER asks a different question: which taxa collectively define a baseline multivariate state associated with later resistance?

PC2 is therefore interpreted as a weighted ecological configuration rather than as a list of individually causal organisms.

Positive PC2-loading families include Lactobacillaceae, Pseudomonadaceae, Propionibacteriaceae, Streptococcaceae, and Beijerinckiaceae.

Negative PC2-loading families include Paracoccaceae, Bacteroidaceae, Anaerovoracaceae, Ruminococcaceae, Saccharimonadaceae, and Caulobacteraceae.

These are axis-defining families. They are not claimed to be individually protective, harmful, causal, or independently significant.

## 10. MMER changes the temporal direction of inference

The source studies mainly characterize differences between timepoints.

MMER instead asks whether the initial ecological condition is associated with the magnitude of subsequent change:

`M(T1) -> d[M(T1), M(T2)]`

This is an initial-state dependence question and places the work closer to ecological theories of resistance, stability, response heterogeneity, and context dependence than to conventional microbiome differential-abundance analysis.

## 11. MMER reframes, rather than contradicts, the original findings

The MMER interpretation should not be positioned as evidence that the original studies were wrong.

The original studies characterize average temporal and intervention-associated community behavior.

MMER interrogates heterogeneity around those average trajectories.

A useful synthesis is:

> The original studies describe the average trajectory; MMER asks why individual quarters differ around that trajectory.

This framing makes MMER complementary to the source studies and allows their datasets to support a new ecological question that was not the primary objective of either original paper.

## 12. MMER inherits important limitations from the source datasets

The data were not prospectively collected to test the final MMER ecological-resistance hypothesis.

Important inherited limitations include:

- small sample size, particularly the Italy cohort;
- different farms and geographic settings;
- different 16S variable regions;
- different intervention regimes;
- different exact postpartum sampling schedules;
- possible differences in extraction and library preparation;
- the low-biomass nature of healthy milk;
- reliance on relative-abundance 16S profiles;
- no experimental manipulation of the baseline ecological state.

Accordingly, MMER should be framed as a hypothesis-driven retrospective reanalysis identifying evidence consistent with baseline-dependent ecological resistance, not as proof of a causal mechanism.

## Core novelty statement

> **Whereas the original studies asked how the healthy mammary microbiome changes across dry-cow interventions and the dry period, MMER asks whether inter-quarter variation in baseline community organization is associated with the magnitude of each quarter's subsequent ecological displacement.**

**Conceptual transformation:** from average perturbation effects to baseline-dependent ecological susceptibility/resistance.

## Source studies

- Biscarini F. et al. (2020). *A Randomized Controlled Trial of Teat-Sealant and Antibiotic Dry-Cow Treatments for Mastitis Prevention Shows Similar Effect on the Healthy Milk Microbiome.* Frontiers in Veterinary Science.
- Derakhshani H. et al. (2018). *Composition of the teat canal and intramammary microbiota of dairy cows subjected to antimicrobial dry cow therapy and internal teat sealant.* Journal of Dairy Science.

*This comparison is intended to define conceptual and analytical novelty. It does not imply that the original studies were deficient; their analyses were appropriate to their stated objectives.*