# MMER

**Mammary Metataxonomy Ecological Resistance**

MMER is a reproducible reanalysis framework for quantifying **ecological resistance of the bovine mammary microbiome** to dry-cow interventions using public longitudinal 16S rRNA amplicon data.

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

Cow order repeats within each block: **270 → 355 → 366 → 321 → 365**, four quarters per cow.

The original aliases `270AS_S1`, `270AD_S2`, `270PD_S3`, and `270PS_S4` provide the anatomical key used to reconstruct the repeated four-position quarter order. The inferred mapping is generated and validated by `scripts/00_build_metadata.py` rather than being manually maintained.

## Core scientific questions

1. **Q1 — How far does each mammary community move from its own baseline?**
   Quantify quarter-level ecological displacement from T1 using Bray–Curtis and Aitchison distances, with the untreated quarter as a within-cow temporal reference.
2. **Q2 — What ecological dimensions constitute resistance/perturbation?**
   Decompose the response into community membership retention/turnover, diversity restructuring, richness, evenness, dominance, and related dimensions.
3. **Q3 — Does baseline ecological state predict later perturbability?**
   Test whether baseline richness, Shannon diversity, evenness, dominance, and multivariate community composition predict subsequent displacement.

## Key distinction from the original analysis

The original study primarily asked whether treatment groups differed in microbiome composition/diversity. MMER asks a different question: **how far did each individual quarter move from its own baseline, what ecological dimensions underlie that response, and which baseline ecological features predict greater resistance or perturbability?**

For quarter `q` of cow `i` at post-baseline time `t`:

`D(i,q,t) = d(M(i,q,t), M(i,q,T1))`

Lower displacement indicates greater ecological resistance. MMER retains the raw distance metrics as primary outcomes rather than forcing all distances into a single arbitrary resistance score.

## Production status

The full 60-sample production workflow has now completed end-to-end.

### DADA2

- 8,816,802 input read pairs
- 6,492,116 filtered reads (73.6%)
- 5,633,741 merged reads (86.8% of filtered)
- 5,384,222 non-chimeric reads
- 9,859 non-chimeric ASVs
- median ASV length: 403 bp

### Taxonomy

- 7,736 bacterial ASVs
- 5,074,031 bacterial reads (94.2% of retained reads)
- ~81% of retained reads assigned to genus

### Current ecological findings

- **Q1:** Aitchison displacement shows the strongest treatment/quarter-condition signal. Cloxacillin-assigned quarters showed ~+20 Aitchison units of excess displacement relative to untreated controls; the 20-quarter aggregated contrast remained significant after BH correction (`q = 0.019`). Bray–Curtis displacement showed the same general direction but did not show a significant treatment effect.
- **Q2:** Cloxacillin-associated trajectories show a coherent multidimensional perturbation profile (reduced core retention, increased turnover, larger richness/dominance restructuring), but no individual component survives full-family BH-FDR correction. Q2 is therefore interpreted as ecological decomposition rather than independent treatment-effect evidence for each component.
- **Q3:** Baseline richness, Shannon diversity, evenness, dominance, CLR-PC1, and CLR-PC2 robustly predict later Bray–Curtis displacement. These relationships survive repeated-measures correction, 20-quarter aggregation, T3-only sensitivity analysis, and BH-FDR correction. No single genus explains the effect robustly across all sensitivity analyses, supporting an emergent whole-community interpretation.

See [`docs/production_analysis_status_2026-08-12.md`](docs/production_analysis_status_2026-08-12.md) for the current full results summary and sensitivity-analysis details.

## Repository layout

```text
MMER/
├── config/                 # analysis parameters
├── data/
│   ├── metadata/           # ENA metadata + reconstructed sample map
│   └── raw/                # FASTQ files (gitignored)
├── docs/                   # design, metadata, and analysis notes
├── manuscript/             # manuscript working files
├── results/
│   ├── qc/
│   ├── asv/
│   ├── resistance/
│   ├── figures/
│   └── tables/
├── scripts/                # numbered reproducible analysis steps
├── tests/                  # metadata integrity checks
├── workflow/               # Snakemake entry point
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
- Raw ENA run report is preserved at `data/metadata/PRJEB38332_ena_report.tsv`.

## Reproducibility principles

- Never infer treatment directly from a FASTQ filename without passing the metadata validator.
- Preserve ENA accessions and original aliases.
- Treat **quarter nested within cow** as the biological unit for longitudinal displacement.
- Use the untreated quarter as the within-cow reference for background dry-off/calving-associated ecological change.
- Separate **ecological resistance** from antimicrobial resistance terminology.
- Report low-biomass limitations explicitly; the source study did not include sequenced extraction-negative controls.
- Treat treatment and anatomical quarter as confounded in the source design because each treatment is fixed to the same quarter position across cows.
- Treat the multidimensional resistance score as exploratory and the individual Q2 dimensions as descriptive unless supported by inferential sensitivity analyses.
- Treat baseline-predictor analyses as discovery/hypothesis-generating because the discovery cohort contains only five cows.

## Status

**Production pipeline:** complete for all 60 samples.

**Q1:** inferential hardening complete.

**Q2:** multidimensional decomposition and sensitivity analysis complete.

**Q3:** repeated-measures, 20-quarter, T3-only, depth-diagnostic, PCA-loading, direct-genus, and FDR sensitivity analyses complete.

**Next:** finalize manuscript figures/tables and write the Results/Discussion around ecological resistance as a baseline-dependent community phenotype.
