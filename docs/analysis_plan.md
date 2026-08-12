# MMER analysis plan — Paper 1: ecological resistance

## Primary question
How much does quarter-level mammary microbial community structure change from its own pre-intervention state, and how does that displacement differ across interventions within the same cow?

## Biological unit
A **mammary quarter nested within cow**, observed longitudinally at T1, T2, and T3.

## Aim 1 — Quantify ecological resistance
Primary outcome: baseline-relative Bray-Curtis displacement from T1 to T2 and T1 to T3.

Secondary outcomes:
- Aitchison distance (compositionally appropriate sensitivity analysis)
- taxon retention
- absolute change in Shannon/Simpson/richness
- dominance stability

Interpretation is deliberately metric-specific. Lower distance/change means greater resistance; a composite score will not be introduced unless supported empirically.

## Aim 2 — Characteristics associated with resistance
Candidate baseline ecological predictors:
- richness
- Shannon diversity
- Simpson diversity
- evenness (to add)
- dominance
- baseline community composition / ordination coordinates
- abundance of dominant/core taxa

Analysis must respect n=5 cows. Prioritize within-cow contrasts, effect sizes, permutation/bootstrap uncertainty, and exploratory inference over high-dimensional predictive ML.

The untreated quarter estimates background ecological change associated with dry-off, calving, and early lactation. For each treated quarter, calculate excess displacement relative to that cow's untreated quarter at the same timepoint.

## Aim 3 — Dimensions of ecological resistance
Assess whether compositional, diversity, retention, and dominance dimensions covary or dissociate. Exploratory correlation/PCA may be used only after inspecting metric distributions and with the limited effective sample size clearly acknowledged.

## Essential robustness analyses
1. Bray-Curtis vs Aitchison.
2. T1→T2 and T1→T3 separately.
3. ASV filtering sensitivity.
4. Potential low-biomass contaminants: prevalence/abundance sensitivity filters, with explicit acknowledgment that sequenced extraction blanks are unavailable.
5. Read-depth sensitivity / rarefaction only where the estimand requires it; do not rarefy by default for compositional distance.
6. Report individual cow/quarter trajectories, not only group centroids.

## Scope boundary
This paper addresses **ecological resistance (magnitude of perturbation)**, not true ecological resilience/recovery kinetics. T2 and T3 can show persistence of displacement, but with only two post-baseline observations the study should not overclaim recovery rate.
