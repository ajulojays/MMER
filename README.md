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

## Core scientific aims

1. **Quantify quarter-level mammary ecological resistance nested within cow.**
   Measure the magnitude of within-quarter ecological displacement from the pre-intervention baseline.
2. **Identify ecological characteristics associated with resistance.**
   Test whether baseline diversity, evenness, dominance, composition, and individual taxa are associated with subsequent displacement, while using the untreated quarter within each cow as the physiological-transition reference.
3. **Resolve the dimensions of mammary ecological resistance.**
   Determine whether compositional resistance, diversity resistance, taxon retention, dominance stability, and related ecological properties represent a common or multidimensional response.

## Key distinction from the original analysis

The original study primarily asked whether treatment groups differed in microbiome composition/diversity. MMER asks a different question: **how far did each individual quarter move from its own baseline, how heterogeneous was that response, and which baseline ecological features were associated with greater resistance?**

For quarter `q` of cow `i` at post-baseline time `t`:

`D(i,q,t) = d(M(i,q,t), M(i,q,T1))`

Lower displacement indicates greater ecological resistance. MMER retains the raw distance metrics as primary outcomes rather than forcing all distances into a single arbitrary resistance score.

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

Then run the ASV workflow after setting the DADA2 trimming parameters in `config/config.yaml`:

```bash
Rscript scripts/02_dada2_asv.R
Rscript scripts/03_resistance_metrics.R
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
- Treat baseline-predictor analyses as exploratory because the discovery cohort contains only five cows.

## Status

**Phase 0:** metadata reconstruction and repository setup — complete.

**Next:** download reads, run QC/DADA2, reproduce selected original-study summaries, then calculate baseline-relative ecological resistance.
