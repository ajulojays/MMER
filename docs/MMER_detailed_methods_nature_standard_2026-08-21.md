# Methods

## Study design and analytical framework

### Overview

This study was designed as a retrospective, hypothesis-driven reanalysis of publicly available longitudinal bovine mammary microbiome datasets. The central analytical objective was to determine whether the microbial community state present **before** a major longitudinal transition was associated with the degree to which that same mammary ecological unit subsequently retained its baseline composition. We refer to this phenotype as **ecological resistance**. The study was organized deliberately into a discovery/development phase and an external validation/generalization phase. Italy (Biscarini *et al.*) and Manitoba (Derakhshani *et al.*) were used to formulate and estimate the primary association in healthy mammary-quarter trajectories, whereas Porcellato and Van Beeck were used as independent 16S external validation cohorts and Patangia was used as a cross-platform shotgun-metagenomic external generalization cohort [1–5].

A key design principle was that the unit of longitudinal comparison was preserved from the source experiment whenever it could be reconstructed unambiguously. Thus, a mammary quarter was followed against **the same quarter** in the Italy, Manitoba and Porcellato analyses; the selected longitudinal milk trajectory was followed within cow in Van Beeck; and the same front-right quarter/cow trajectory was followed across M0, M2, M4 and M6 in Patangia [1–5]. No cross-sectional sample was allowed to substitute for a missing within-unit follow-up. The data reconstruction therefore preceded phenotype construction and statistical modelling.

The MMER analysis was not intended to reproduce the original publications' inferential questions. Rather, the source studies were first reconstructed so that animal identity, quarter identity, timepoint, treatment/exposure, farm/dairy and sequencing provenance agreed with the source design as closely as the deposited metadata permitted. Only after this reconstruction was validated were the data transformed into a common ecological-resistance framework. This distinction is important: **reconstruction fidelity** refers to recovering the biological sample map and longitudinal design of the source study; **analytical harmonization** refers to the new, common MMER processing and statistical representation used to test a different hypothesis.

### Why the MMER question differs from the source-study questions

The source publications primarily asked whether microbiome diversity or composition differed by treatment, sampling period, farm, mammary niche or lactation stage [1–5]. MMER instead treats heterogeneity among individual longitudinal trajectories as the phenotype. The conceptual contrast is:

**Source-study framing:** exposure/time/farm → average microbiome difference.

**MMER framing:** baseline microbial state → magnitude of subsequent within-unit ecological displacement.

This change in direction of inference is central. A significant time effect establishes that communities change on average. It does not determine whether communities that start from different ecological states are more or less likely to retain their original composition. Likewise, a weak average treatment effect does not imply that all mammary communities respond equally. MMER therefore asks an initial-state-dependence question that is complementary to, rather than a replacement for, the source publications.

## Public data identification, study eligibility and reconstruction strategy

### General eligibility criteria

Candidate studies were screened for the following features: (i) bovine mammary milk or colostrum microbiome data; (ii) repeated sampling of an identifiable biological unit; (iii) recoverable animal identity and, where relevant, quarter identity; (iv) a defensible baseline state that preceded the follow-up state used to define ecological resistance; and (v) sufficient sequence and metadata information to reconstruct a taxonomic community profile. Studies were not required to use identical sequencing platforms, amplicon regions, farms, treatment regimes or sampling intervals, because the purpose of external validation was to test phenomenon-level generalization under realistic heterogeneity rather than to create a synthetic single experiment.

For every study, reconstruction followed four audit stages. First, the source publication and supplementary material were used to define the intended animal population, longitudinal unit, treatment/exposure structure, sampling schedule and sequencing protocol. Second, public-repository records were reconciled with the paper using accession numbers, sample aliases, run identifiers and sample titles. Third, a locked trajectory manifest was constructed, requiring the same biological unit at the necessary timepoints. Fourth, sample counts and grouping structure were checked against the source paper whenever the source paper reported a directly comparable number. Disagreements caused by independent sequence reprocessing or by the additional requirement for complete longitudinal pairs were retained and documented rather than silently altered to force numerical agreement.

### Study-role assignment

The final evidentiary structure comprised five biological cohorts. Italy and Manitoba constituted the discovery/development analysis. Porcellato and Van Beeck constituted independent 16S external validations. Patangia constituted a shotgun-metagenomic external generalization test. This hierarchy was frozen before interpreting cohort-specific secondary analyses. In particular, a null external result was not replaced by a later time window simply because the later window produced a smaller *P* value.

## Italy cohort: source design, public-data reconstruction and selection

### Source experiment

The Italy dataset originated from Biscarini *et al.* [1], BioProject **PRJEB38332**. The source experiment enrolled five healthy second-parity Holstein-Friesian cows from a commercial herd in northern Italy. Cows were selected as clinically healthy and had somatic cell counts (SCC) below 200,000 cells ml−1. Each cow contributed all four mammary quarters. The design intentionally assigned a different dry-cow condition to each anatomical quarter: front-left (FL), untreated control; front-right (FR), internal teat sealant (bismuth subnitrate); rear-right (RR), cephalonium; and rear-left (RL), benzathine cloxacillin [1]. Samples were collected at three longitudinal stages spanning dry-off through the subsequent lactation: dry-off baseline, calving/colostrum and approximately 5 days in milk. The source study amplified the V3–V4 region of the bacterial 16S rRNA gene and sequenced paired-end 2 × 250-bp reads on an Illumina MiSeq [1].

The original publication focused on whether antibiotic and non-antibiotic dry-cow treatments measurably altered the healthy milk microbiome and reported stronger temporal/lactation-stage effects than treatment effects [1]. Because treatment was fixed to anatomical quarter in the deposited design, treatment and anatomical quarter cannot be separated perfectly in a retrospective reanalysis; MMER therefore does not treat the Italy data as a clean randomized estimate of a treatment effect.

### Reconstruction from ENA metadata

The PRJEB38332 public sample map was reconstructed directly from the ENA report and sample aliases. The repository contained 60 sequencing runs, corresponding exactly to 5 cows × 4 quarters × 3 timepoints. Reconstruction was implemented in `scripts/00_build_metadata.py` and was made deliberately explicit rather than inferred from row order without validation.

Deposited aliases encoded a sample-number sequence. Samples **S1–S20** were mapped to the first longitudinal block (T1; dry-off baseline), **S21–S40** to T2 (calving), and **S41–S60** to T3 (5 DIM). Within each 20-sample block, cows occurred in the repeated order 270, 355, 366, 321 and 365. Within each cow, four successive records mapped to FL, FR, RR and RL. The anatomical key was cross-checked against aliases encoding `AS`, `AD`, `PD` and `PS`; for example, the first four aliases beginning with cow 270 agreed with the expected `270AS`, `270AD`, `270PD` and `270PS` pattern. The reconstruction script also checked that the numeric cow prefix embedded in every alias agreed with the inferred cow identifier.

Integrity checks required: exactly 60 records; 60 unique sample numbers covering S1–S60 with no gaps; all five expected cow identifiers; exactly three timepoints for every cow-quarter combination; exactly four quarters for every cow-timepoint combination; and agreement between inferred cow and deposited alias prefix. Paired FASTQ URLs were separated into R1 and R2 fields and written to a frozen sample map. This approach allowed every run accession to be traced to one cow, one anatomical quarter, one treatment/quarter condition and one timepoint.

### Independent sequence processing and taxonomic reconstruction

The full 60-sample Italy sequence set was reprocessed independently rather than reusing the source publication's feature table. The production DADA2 workflow yielded 8,816,802 input read pairs; 6,492,116 pairs passed filtering (73.6%); 5,633,741 reads merged; and 5,384,222 non-chimeric reads remained, corresponding to approximately 61.1% of input. A total of 9,859 non-chimeric ASVs were inferred. Taxonomy was assigned against SILVA 138.2; 7,736 ASVs were classified as bacterial and 5,074,031 reads were retained as bacterial reads. DADA2 was used because it models amplicon-sequencing errors and resolves exact sequence variants rather than clustering all reads into fixed-identity OTUs [9].

For the cross-study discovery analysis, the ecological representation was harmonized at bacterial-family level. Organelle-associated labels (mitochondrial, chloroplast, plastid or other organelle annotations) were removed. The final discovery input was not the complete 60-sample Italy experiment, but the healthy quarter trajectories meeting the pre-specified phenotype-comparability rules described below.

### Italy trajectory selection for the discovery cohort

The discovery question was intentionally restricted to healthy mammary-quarter trajectories so that the initial association was not defined by overt disease. A trajectory required a valid baseline family profile and the harmonized follow-up state used for the pooled Italy–Manitoba phenotype. Twenty Italy quarter trajectories from five cows met the frozen primary criteria. Because all four anatomical quarters from a cow can contribute, these 20 trajectories are not 20 independent animals; cow identity was therefore retained explicitly for cluster-robust inference.

## Manitoba cohort: source design, reconstruction and selection

### Source experiment

The Manitoba dataset originated from Derakhshani *et al.* [2]. Eleven late-pregnant Holstein cows were initially selected from a commercial Manitoba dairy under stringent healthy-quarter criteria: no clinical mastitis or antibiotic treatment during the preceding four months, four functional quarters, no abnormal milk, California Mastitis Test score of 0, no udder edema/redness, and quarter SCC below 200,000 cells ml−1 at dry-off [2]. At dry-off, all quarters received blanket dry-cow therapy containing 200,000 units penicillin G plus 400 mg novobiocin per 10-ml tube together with an internal teat sealant. Samples were collected at dry-off and immediately after calving; the mean dry period was 58.77 ± 8.35 days [2].

The original study sampled two mammary niches: teat-canal swabs and corresponding mammary secretions. The milk/colostrum component consisted of quarter milk before dry-cow therapy and quarter colostrum immediately postpartum. V1–V2 16S rRNA amplicons were generated with modified F27/R357 primers and sequenced using an Illumina MiSeq 600-cycle kit. Public reads were deposited under SRA accessions **SRR7267363–SRR7267478** [2]. The source bioinformatics pipeline used FLASH and UPARSE at 97% OTU identity, and its primary ecological questions concerned compositional shifts and persistence across the dry period after antimicrobial dry-cow therapy [2].

### Reconstruction of mammary-secretions trajectories

For MMER, teat-canal swabs were excluded because they represent a distinct ecological niche. Only mammary-secretions samples were eligible. The baseline state was defined as the pre-DCT quarter milk sample at dry-off and the follow-up state as the corresponding postpartum colostrum from the **same mammary quarter**. This mapping follows the original study's paired sampling design and prevents comparisons between different quarters or between teat-canal and intramammary niches.

The original paper reported that, after its quality filtering, 58 teat-canal swabs and their corresponding mammary-secretions samples from nine cows were retained for downstream analysis; the mammary-secretions component comprised **29 milk samples and 29 colostrum samples** [2]. The MMER reconstruction recovered **29 complete same-quarter milk→colostrum trajectories from nine cows**, matching the source publication's retained mammary-secretion structure. This numerical agreement was used as an important reconstruction check: each trajectory required one pre-DCT milk profile, one postpartum colostrum profile, the same cow identifier and the same quarter identifier.

### Harmonization for the discovery analysis

The source Manitoba paper used a V1–V2 OTU workflow, whereas Italy used V3–V4 data. Direct pooling of sequence-level ASVs across these variable regions would be biologically and technically inappropriate. Therefore, MMER harmonized both discovery datasets at **bacterial-family level**, a deliberately coarser taxonomic representation that is more comparable across amplicon regions than sequence variants. The analysis used the locked family count profiles in the harmonized MMER data layer. Samples were normalized into the common `count_3k` analysis representation before construction of the pooled family matrices. Organelle-associated taxonomic labels were excluded using the same rule applied to Italy.

The exact upstream source-paper OTU pipeline was not copied into MMER, because doing so would preserve incompatible feature definitions between studies. Instead, the source design and biological pairing were reconstructed faithfully, while downstream taxonomic representation was harmonized to answer the new cross-study ecological question. This separation between biological reconstruction and analytical harmonization is intentional and should not be interpreted as an attempt to reproduce the original paper's OTU-level numerical tables.

## Construction of the Italy–Manitoba discovery cohort

### Harmonized timepoint definition

For each eligible trajectory, MMER defined a baseline state (`T1`) and a phenotype-defining follow-up (`T2`) using study-specific sampling occasions that represented the intended dry-off→early-postpartum transition. The trajectory identifier combined study, cow and quarter so that records from different studies or different quarters could never be joined accidentally. The frozen discovery cohort contained **49 complete quarter trajectories from 14 cows**: 20 trajectories from 5 Italy cows and 29 trajectories from 9 Manitoba cows.

### Family-profile construction

For each study, family counts for the same trajectory and harmonized timepoint were aggregated into a wide trajectory × family matrix. Missing families in a sample were represented as zero counts. A common family universe was used so that T1 and T2 matrices had identical columns and order. If replicate records existed at a harmonized timepoint, the locked workflow summarized the normalized counts by their mean before matrix construction. A trajectory entered the primary analysis only if both T1 and T2 existed after harmonization.

## Ecological resistance phenotype

### Bray–Curtis displacement

The primary response variable was based on Bray–Curtis dissimilarity [6]. For two non-negative family-abundance vectors, **a** and **b**, with elements indexed by family *j*, Bray–Curtis dissimilarity was calculated as

$$BC(\mathbf{a},\mathbf{b})=\frac{\sum_{j=1}^{p}|a_j-b_j|}{\sum_{j=1}^{p}(a_j+b_j)}.$$

Here, *p* is the number of bacterial families in the harmonized family universe; *a* is the family profile of a trajectory at baseline; and *b* is the family profile of that same trajectory at follow-up. Bray–Curtis is 0 when the two abundance profiles are identical and approaches 1 as their abundance-weighted composition becomes increasingly dissimilar. It emphasizes changes in relative abundance and does not award similarity merely because the same number of taxa are present.

For trajectory *i*, baseline-to-follow-up displacement was

$$D_i=BC(\mathbf{x}_{i,T1},\mathbf{x}_{i,T2}).$$

The symbols have a simple interpretation. *i* identifies one longitudinal mammary unit (for the discovery analysis, one specific quarter); $\mathbf{x}_{i,T1}$ is that quarter's family-level community at baseline; $\mathbf{x}_{i,T2}$ is the same quarter later; and $D_i$ is how far the community moved compositionally.

### Resistance transformation

Ecological resistance was defined as

$$R_i=1-D_i=1-BC(\mathbf{x}_{i,T1},\mathbf{x}_{i,T2}).$$

Thus, $R_i$ is the outcome used in the primary regression. A value near 1 means the follow-up community retained a composition very similar to its own baseline (high resistance); a value near 0 means large baseline-relative displacement (low resistance, or high susceptibility). This transformation changes only the direction of interpretation; it does not alter the information in the Bray–Curtis distance.

This outcome is fundamentally different from Shannon diversity. Shannon diversity describes richness/evenness **within one sample**, whereas $R_i$ compares two observations from the same biological unit. Two samples may have identical Shannon diversity but very different taxa or abundance assignments; such a change can yield substantial Bray–Curtis displacement even when alpha diversity is unchanged.

## Baseline ecological architecture

### Why baseline features were restricted to T1

All predictor construction used **baseline data only**. No follow-up abundance, follow-up diversity or outcome information was allowed to influence family prevalence filtering, CLR representation or PCA scores. This temporal separation preserves the direction of the hypothesis: starting state → subsequent resistance. It also prevents leakage in which information from the response period contributes to the predictor definition.

### Baseline prevalence filter

For each family *j*, baseline prevalence was

$$Prev_j=\frac{1}{N}\sum_{i=1}^{N}I(x_{ij,T1}>0),$$

where *N* is the number of baseline trajectories and $I(\cdot)$ is an indicator taking value 1 when the family is observed. Families present in at least **10% of baseline trajectories** were retained. This threshold reduces the influence of extremely sparse features whose ordination coordinates would otherwise be driven largely by zero patterns. The discovery analysis retained **161 bacterial families** after this step.

### Centered log-ratio transformation

Microbiome sequencing data are compositional: observed counts are constrained by sequencing depth and represent relative information rather than unconstrained absolute abundances [7,8]. To place the baseline family vectors in an appropriate log-ratio geometry, a centered log-ratio (CLR) transformation was applied. In the discovery analysis a pseudocount of 0.5 was added before taking logarithms. For baseline vector $\mathbf{x}_i=(x_{i1},...,x_{ip})$, the CLR coordinate for family *j* was

$$CLR(x_{ij})=\log(x_{ij}+c)-\frac{1}{p}\sum_{k=1}^{p}\log(x_{ik}+c),$$

where $c=0.5$ in the frozen discovery analysis. Equivalently, each transformed abundance is the logarithm of the family's abundance relative to the geometric mean of all retained families in that sample. Positive CLR values indicate that a family is relatively enriched compared with the sample's geometric mean; negative values indicate relative depletion. CLR transformation does not create absolute-abundance information and does not remove the limitations of low-biomass microbiome data; it simply represents relative compositions in a log-ratio space appropriate for multivariate analysis [7,8].

### Removing study centroids before pooled PCA

Italy and Manitoba differ in geography, source study, 16S variable region, laboratory workflow, exact longitudinal interval and host population. A naïve pooled PCA could therefore allocate its largest axes to between-study technical or biological separation rather than to within-study ecological heterogeneity. To reduce that problem, the mean CLR vector for each study was subtracted from every baseline sample in that study before pooled PCA.

Let $\mathbf{z}_{is}$ denote the CLR vector for trajectory *i* in study *s*, and let

$$\bar{\mathbf{z}}_s=\frac{1}{N_s}\sum_{i\in s}\mathbf{z}_{is}$$

be the study-specific centroid. The centered vector was

$$\mathbf{z}^{*}_{is}=\mathbf{z}_{is}-\bar{\mathbf{z}}_s.$$

In simple terms, this asks: *where does this quarter sit relative to other quarters in the same study?* rather than allowing the first axis merely to encode *Italy versus Manitoba*. Numerical QC confirmed that retained family means were approximately zero within each study after centering. Zero-variance families after centering were removed.

### Principal-component analysis

PCA was performed on the within-study-centered CLR matrix without additional scaling of individual family columns. If $\mathbf{Z}^{*}$ is the centered $N\times p$ matrix, PCA decomposes its covariance structure into orthogonal directions. Each principal component (PC) is a weighted linear combination of the retained family coordinates. The PC score for trajectory *i* on component *k* can be written as

$$PC_{ik}=\sum_{j=1}^{p}w_{jk}z^{*}_{ij},$$

where $w_{jk}$ is the loading (weight) of family *j* on component *k*. The loading describes how strongly and in what direction that family contributes to the axis; it is not, by itself, an estimate of causal or independent taxon effect.

The proportion of baseline variance captured by PC *k* was

$$VE_k=\frac{\lambda_k}{\sum_{m}\lambda_m},$$

where $\lambda_k$ is the eigenvalue associated with PC *k*. In the frozen discovery analysis PC1 explained **10.03%** and PC2 **9.26%** of the within-study-centered baseline variance, for **19.29% cumulative variance**. PC1 and PC2 were pre-specified as a low-dimensional representation of baseline architecture in the frozen primary analysis. Variance explained describes how much baseline compositional heterogeneity an axis captures; it is not equivalent to how strongly that axis predicts resistance.

For regression, PC1 and PC2 scores were standardized across the discovery trajectories:

$$zPC_{ik}=\frac{PC_{ik}-\overline{PC_k}}{SD(PC_k)}.$$

Consequently, a regression coefficient for `z_PC1` or `z_PC2` represents the expected change in ecological resistance associated with a **one-standard-deviation** difference along that baseline ecological axis.

## Discovery hypothesis and primary regression model

### Primary scientific hypothesis

The discovery hypothesis was: **among healthy mammary-quarter trajectories, does multivariate baseline ecological architecture contain information about how strongly the same quarter retains its baseline composition over the dry-off→early-postpartum transition?**

The primary model was

$$R_{iq}=\beta_0+\beta_1zPC1_{iq}+\beta_2zPC2_{iq}+\beta_3Study_{iq}+\varepsilon_{iq}.$$

Each term has a direct interpretation:

- $R_{iq}$: ecological resistance of quarter trajectory *q* from cow *i*; this is the dependent/outcome variable.
- $\beta_0$: intercept, the model's reference expected resistance when continuous predictors are at zero and the categorical study variable is at its reference level.
- $zPC1_{iq}$: standardized baseline position of that quarter along the first ecological axis.
- $\beta_1$: expected difference in resistance per 1-SD increase in PC1, conditional on PC2 and study.
- $zPC2_{iq}$: standardized baseline position along the second ecological axis.
- $\beta_2$: expected difference in resistance per 1-SD increase in PC2, conditional on PC1 and study.
- $Study_{iq}$: indicator for source study (Italy versus Manitoba), included to account for mean differences in resistance not represented by the within-study-centered PCs.
- $\beta_3$: coefficient for the study contrast.
- $\varepsilon_{iq}$: residual variation in resistance not explained by the modeled predictors.

The **primary inferential target was joint**, not a search for whichever single PC had the smallest *P* value:

$$H_0:\beta_1=\beta_2=0.$$

Under this null hypothesis, baseline architecture represented by PC1 and PC2 provides no additional information about ecological resistance after study adjustment. Rejection of the joint null supports the broader architecture–resistance association even if one component is individually weak.

### Dependence among quarters from the same cow

The 49 discovery trajectories arose from 14 cows, and multiple quarters from one cow share host genetics, physiology, management and environmental exposures. Treating all quarters as independent observations would therefore underestimate uncertainty. The linear mean model was fitted by ordinary least squares, but inference used **cow-clustered CR2** variance estimation with small-sample corrections, implemented using `clubSandwich`. CR2 is a bias-reduced cluster-robust covariance estimator designed for settings with a modest number of independent clusters [12]. Individual coefficient tests used Satterthwaite degrees of freedom; the two-parameter joint PC1+PC2 hypothesis used the HTZ small-sample Wald/F approximation.

The clustering assumption is that dependence may be arbitrary among quarters within a cow but that cows are the independent sampling clusters. CR2 does not eliminate bias from a misspecified mean model, unmeasured confounding or incorrect reconstruction of longitudinal units; it addresses covariance estimation under within-cluster dependence.

### Descriptive model fit

Ordinary $R^2$ and adjusted $R^2$ were reported as descriptive measures of variance accounted for by the full linear model. They were not treated as cluster-robust inferential statistics. Because the full discovery model contains PC1, PC2 and study, its $R^2$ should not be described as the variance explained by PC1+PC2 alone unless a separate partial-$R^2$ calculation is performed.

## Discovery reproducibility model across Italy and Manitoba

A second, pre-specified analysis tested whether the architecture–resistance relationship differed detectably between the two discovery cohorts. The model was

$$R_{iq}=\beta_0+\beta_1zPC1_{iq}+\beta_2zPC2_{iq}+\beta_3Study_{iq}+\beta_4(zPC1\times Study)_{iq}+\beta_5(zPC2\times Study)_{iq}+\varepsilon_{iq}.$$

The joint interaction null was

$$H_0:\beta_4=\beta_5=0.$$

This asks whether the slopes relating baseline PC1/PC2 to ecological resistance are statistically distinguishable between Italy and Manitoba. A non-significant interaction does **not** prove identical biological effects; it indicates insufficient evidence, given the available sample size and clustering, to conclude that the multivariate relationship differs.

## External validation philosophy

### Phenomenon-level validation versus exact-axis transport

The primary external-validation goal was **phenomenon-level replication**: in an independently processed cohort, does that cohort's own low-dimensional baseline ecological architecture predict subsequent within-unit resistance? External cohorts therefore underwent their own baseline prevalence filtering, CLR transformation, cohort-structure centering and PCA. This design tests whether the ecological principle generalizes, not whether an identical numerical PC1 or PC2 from the discovery data exists in every cohort.

Because PCA axes are sign-indeterminate and depend on the feature covariance structure, the sign or taxonomic loading of independently fitted external PC1/PC2 cannot be compared literally with the discovery axes. Exact axis transport would require a different analysis: harmonization to an identical feature set followed by projection of frozen Italy–Manitoba loading vectors into external samples. That stronger transportability test was not claimed here.

### External validation decision rule

For each external cohort, the main question remained the joint contribution of the first two baseline ecological axes after adjustment for design variables required by that cohort. The external model and covariance estimator were chosen **before** interpreting the result based on the longitudinal unit and dependence structure. Cohorts with multiple quarter trajectories per cow used cow-clustered inference; cohorts with one trajectory per cow used independent-observation inference with heteroskedasticity-robust covariance where appropriate.

## Porcellato external validation: reconstruction and analysis

### Source study and public-data structure

Porcellato *et al.* [3] sampled 60 Norwegian Red cows from two Norwegian herds, 30 cows at farm A (NMBU Centre for Livestock Production) and 30 at farm K (Kalnes). Cows were sampled at two occasions during the same lactation, with mean DIM approximately 38–44 days during the first period and 213–222 days during the second. Eleven cows left the study before the second sampling. Cows were considered healthy if no clinical mastitis was observed, no mastitis was recorded in the preceding two weeks and no antimicrobial treatment was ongoing [3]. Quarter milk samples were obtained during regular milking. The source study amplified the V3–V4 16S rRNA region with Uni340F/Bac806R and performed short-read amplicon sequencing. The deposited study is **PRJEB35792** [3].

The source publication included more than 400 quarter samples and reported 403 samples after its own quality filtering, 4,832,201 high-quality sequences and 10,010 sequence variants, of which 8,759 were positively assigned at family level [3]. Its scientific focus was the healthy-udder core microbiota, farm and sampling-period differences, and dysbiosis rather than a baseline-resistance phenotype.

### Reconstructing cow × quarter × sampling-period trajectories

Public sample aliases were parsed to recover cow identifier, anatomical quarter, sampling period and farm. Identifiers such as `Cow6_LF_Sampling2_K` contain the information needed to link a late-lactation sample to the corresponding first-sampling LF quarter in the same cow and farm. A trajectory was defined only when the identical cow-quarter combination was present at Sampling1 and Sampling2. This reconstruction yielded **189 complete Sampling1→Sampling2 quarter trajectories** before the MMER sequencing-depth lock.

The source paper's 403-sample post-QC count and the MMER reprocessed sample count are not expected to be identical because MMER independently denoised the deposited sequence data and then imposed a same-quarter longitudinal pairing requirement. The fidelity criterion was therefore biological: recovered farm, cow, quarter and sampling-period identities had to agree with the source design, and the pairing algorithm could not cross cows, quarters or farms.

### Sequence processing and locked depth criterion

Independent DADA2/taxonomy processing yielded 378 samples, 15,567 ASVs and 213 bacterial families overall. Before viewing the external regression result, the primary analysis required at least **1,000 DADA2 non-chimeric reads at both Sampling1 and Sampling2**. This retained **154 complete trajectories**, 93 from farm A and 61 from farm K, corresponding to 308 samples. The paired resistance matrix contained 211 bacterial families; baseline prevalence filtering retained 50 families for PCA.

Stricter post-taxonomy assigned-family-read thresholds were evaluated only as robustness analyses after the primary set was locked: at least 1,000 assigned family reads at both timepoints retained 143 trajectories, and at least 2,000 retained 121. These sensitivity sets were not substituted for the primary 154-trajectory cohort.

### Porcellato ecological architecture and model

Sampling1 was the baseline and Sampling2 the follow-up. For trajectory *iq*, resistance was

$$R_{iq}=1-BC(\mathbf{x}_{iq,S1},\mathbf{x}_{iq,S2}).$$

Baseline family compositions were CLR transformed after zero handling and centered within farm before an independent Porcellato PCA. The first two PC scores were standardized. The primary model was

$$R_{iq}=\alpha_0+\alpha_1zPC1_{iq}+\alpha_2zPC2_{iq}+\alpha_3Farm_{iq}+\epsilon_{iq}.$$

The primary external-validation null was $H_0:\alpha_1=\alpha_2=0$. Farm was included because the source paper showed strong farm-associated compositional differences [3], and the MMER analysis sought architecture information beyond this mean farm structure. Cows could contribute more than one quarter trajectory, so inference used cow-clustered CR2 with Satterthwaite coefficient tests and an HTZ joint PC1+PC2 test. The independent Porcellato PCA explained 19.779% of baseline variance with PC1 and 8.823% with PC2.

## Van Beeck external validation: reconstruction and analysis

### Source study

Van Beeck *et al.* [5] sampled milk and teat skin at three California dairies across three timepoints: baseline at dry-off before treatment, 7 days later, and 55–75 DIM in the next lactation. Cows were categorized using quarter SCC measured at baseline. Low-SCC cows had SCC <100,000 cells ml−1 in all four quarters; cows entered the high-SCC stratum if at least one quarter exceeded 200,000 cells ml−1. Treatment groups comprised low-SCC controls, high-SCC controls, high-SCC cows treated with cephapirin benzathine (CB), and high-SCC cows treated with ceftiofur hydrochloride (CH) [5]. The source study collected 372 milk and teat-skin samples. It generated V4 16S amplicons sequenced on an Ion Torrent S5; demultiplexed single-end reads were processed with QIIME2/DADA2, trimmed and truncated to 290 bp, and classified against SILVA v138. Raw data are deposited under **PRJEB63336** [5].

The original study primarily evaluated the effects of dairy, timepoint, SCC and cephalosporin use on milk and teat-skin diversity. Dairy and timepoint explained larger fractions of beta-diversity than treatment [5]. MMER used the **milk** longitudinal structure only and did not combine milk and teat-skin communities.

### Reconstruction and locked complete trajectories

Public metadata and the previously processed Van Beeck project layer were reconciled into a locked trajectory table. The canonical external-validation metadata contained **60 complete longitudinal milk trajectories and 180 samples**, with exactly one Baseline, one 7 Days and one 55–75 DIM observation per trajectory. Treatment totals in the locked 60-trajectory set were 14 low-SCC untreated controls, 12 high-SCC untreated controls, 18 CB and 16 CH. All trajectories retained dairy identity.

The locked dataset was required to satisfy: 180 metadata rows; 60 unique trajectory identifiers (`core`); exactly three observations per trajectory; one and only one of each required timepoint; valid dairy; and valid treatment. These checks prevented the duplicate-treatment inconsistency present in an earlier unlocked metadata file from entering the canonical external validation.

### Van Beeck family matrix

The external validation used an existing bacterial ASV table rarefied to 3,000 reads per retained sample. Metadata sample names had to match ASV-table row names, and every retained row was required to sum to exactly 3,000 before organelle filtering. ASVs with taxonomy containing mitochondrial, chloroplast, plastid or organelle labels were excluded. ASVs were aggregated to bacterial family; where a Family annotation was absent, a conservative hierarchical label based on Order, then Class, then Phylum, and finally `Unclassified_Bacteria` was used so that unassigned ASVs were not silently discarded merely because family nomenclature was unavailable. Family counts were converted to within-sample relative abundances before Bray–Curtis calculation.

### Primary Van Beeck outcome

The three-timepoint structure permitted two baseline-relative resistance measurements:

$$R_{i,early}=1-BC(\mathbf{x}_{i,Baseline},\mathbf{x}_{i,7d})$$

and

$$R_{i,long}=1-BC(\mathbf{x}_{i,Baseline},\mathbf{x}_{i,55-75DIM}).$$

To avoid selecting a follow-up according to which gave the strongest result, the frozen external-validation phenotype summarized both prespecified post-baseline states:

$$R_{i,overall}=\frac{R_{i,early}+R_{i,long}}{2}.$$

Thus, the Van Beeck outcome measures mean baseline compositional retention across the early and later observations. The temporal component outcomes were not used as additional primary validation hypotheses in the canonical external-validation-only workflow.

### Van Beeck baseline architecture

Only baseline family profiles entered predictor construction. Families present in at least 10% of baseline trajectories were retained. Zeros were replaced using a data-derived pseudocount equal to one-half the smallest positive baseline relative abundance, after which each composition was re-closed to sum to one and CLR transformed. The CLR matrix was centered within dairy before PCA so that the dominant axes would not merely reproduce the three-dairy separation. Zero-variance features were removed. PCA was fitted independently in Van Beeck and PC1/PC2 scores were standardized.

### Van Beeck primary external-validation model

The model was

$$R_{i,overall}=\gamma_0+\gamma_1zPC1_i+\gamma_2zPC2_i+\mathbf{\gamma}_T^T Treatment_i+\mathbf{\gamma}_D^T Dairy_i+\epsilon_i.$$

Here, $R_{i,overall}$ is overall resistance for trajectory *i*; $zPC1_i$ and $zPC2_i$ are standardized baseline ecological coordinates; `Treatment` is a categorical design covariate with four levels; `Dairy` is a categorical design covariate with three levels; $\mathbf{\gamma}_T$ and $\mathbf{\gamma}_D$ denote the sets of dummy-variable coefficients required to represent those categorical factors; and $\epsilon_i$ is unexplained trajectory-level variation. Treatment and dairy were **adjustment variables**, not the primary hypothesis.

The external-validation null was

$$H_0:\gamma_1=\gamma_2=0.$$

Because each outcome row represented one independent trajectory and there was not more than one outcome row per cow/core in the endpoint-specific model, cow-clustered CR2 was not used. The canonical inferential stack comprised a classical nested partial-F test comparing the covariate-only model (`treatment + dairy`) with the full model (`PC1 + PC2 + treatment + dairy`), heteroskedasticity-consistent HC3 covariance for individual coefficients and a joint two-degree-of-freedom Wald test, a design-aware permutation test, and a stratified trajectory bootstrap [14].

### Freedman–Lane permutation in Van Beeck

The permutation analysis followed the Freedman–Lane residual-permutation principle [13]. First, the reduced model containing treatment and dairy but not PC1/PC2 was fitted. Let $\hat{y}^{(0)}_i$ be its fitted value and $e_i^{(0)}$ its residual. Residuals were shuffled **within dairy × treatment strata**, thereby retaining the main design structure while breaking the association between baseline ecological coordinates and resistance under the null. A pseudo-outcome was constructed:

$$y_i^{*(b)}=\hat{y}^{(0)}_i+e_{\pi_b(i)}^{(0)},$$

where $\pi_b$ denotes the within-stratum permutation for iteration *b*. The reduced and full models were refitted to each pseudo-outcome and the partial-F statistic recorded. With 9,999 permutations, the Monte Carlo *P* value was

$$p_{perm}=\frac{1+\sum_{b=1}^{B}I(F_b\geq F_{obs})}{B+1}.$$

The addition of 1 in numerator and denominator prevents a reported *P* value of exactly zero and yields a valid randomization-style estimate under the implemented permutation scheme.

### Stratified bootstrap in Van Beeck

Uncertainty in PC coefficients and partial $R^2$ was assessed using 5,000 complete-trajectory bootstrap resamples. Trajectories were sampled with replacement within dairy × treatment strata, maintaining the observed group structure. The Van Beeck PCA coordinates were held fixed during this bootstrap, so the interval estimates quantify regression uncertainty conditional on the independently learned ecological coordinate system. Percentile 95% intervals were calculated from the 2.5th and 97.5th percentiles of the valid bootstrap distribution. Because nested-model partial $R^2$ is non-negative by construction, its bootstrap interval was interpreted as effect-size uncertainty and not as an additional null-hypothesis test.

## Patangia shotgun-metagenomic external generalization

### Source study

Patangia *et al.* [4] enrolled 24 healthy dairy cows from five farms in County Cork, Ireland and followed them from colostrum to months 2, 4 and 6 of lactation. Cows belonged to three dry-cow-therapy groups: no-antibiotic dry-off plus teat sealant (NOAB; n=9), Cephaguard/cefquinome (CEF; n=8), or Ubro Red (UBRO; n=7). Samples were collected from the **front-right quarter** of each cow. Colostrum was collected within the first hour postpartum (M0), followed by milk at M2, M4 and M6, yielding **96 samples = 24 cows × 4 timepoints** [4]. The source study used whole-genome shotgun sequencing to profile both microbiota and resistome.

The publication reported 494,107,934 raw reads and 108,774,030 microbial reads after trimming and bovine-host removal, with median microbial depth 310,854 reads per sample. Its bioinformatics workflow included read QC/trimming, host subtraction and Kraken2-based taxonomic classification [4,11]. The source scientific question concerned longitudinal effects of antibiotic versus non-antibiotic dry-cow therapy on microbiota diversity and antimicrobial-resistance genes, not ecological-resistance prediction.

### Reconstruction from Kraken2 reports

MMER reconstructed family profiles directly from the Kraken2 report outputs. For every report, the run accession was derived from the filename. Only rows with Kraken rank code `F` were retained. The family name was stripped of indentation and clade counts were parsed as non-negative integer counts. The resulting long table contained **45,868 run × family rows**, **96 unique runs**, **24 cows**, **four timepoints** and **852 unique family labels**. No duplicated run × family combinations were present.

Family rows were joined to a locked analysis manifest by run accession to recover `cow_id`, `timepoint` and sample alias. Each cow was required to contribute exactly one sample at M0, M2, M4 and M6. A complete sample × family matrix was then constructed, filling absent families with zero. Obvious organellar labels containing mitochondrial, chloroplast, plastid or organelle terms were excluded. Counts were converted to relative abundance by sample; all 96 rows were required to have positive depth, finite values and row sums of one.

The reconstructed family-assigned depths were substantial but heterogeneous. At M0 the minimum, median and maximum family depths were approximately 100,293, 213,892 and 1,503,773; corresponding ranges were 69,847–1,519,538 at M2, 184,572–2,194,727 at M4 and 160,959–2,291,040 at M6. These depths refer to the family-level counts represented in the Kraken-derived matrix, not necessarily to all microbial reads reported in the source publication.

### Primary Patangia window and architecture

The external generalization window was frozen as **M0→M2** before later intervals were inspected as secondary/sensitivity analyses. M0 is postpartum colostrum and M2 the first subsequent milk timepoint, providing a clear baseline-to-follow-up trajectory. M0 contained 780 observed families; 680 families passed the ≥10% M0 prevalence filter. After zero replacement and CLR transformation, an independent M0 PCA was fitted. PC1 and PC2 explained 8.958% and 6.580% of baseline variation, respectively.

For cow *i*, the primary resistance was

$$R_i^{M0\rightarrow M2}=1-BC(\mathbf{x}_{i,M0},\mathbf{x}_{i,M2}).$$

The primary model was

$$R_i=\delta_0+\delta_1zPC1_i+\delta_2zPC2_i+\epsilon_i,$$

with joint null $H_0:\delta_1=\delta_2=0$. There was one trajectory-level outcome per cow, so a cow-clustered covariance estimator would create 24 singleton clusters and add no meaningful within-cow correlation adjustment. Accordingly, ordinary model inference and HC3 heteroskedasticity-robust covariance were used. Later M2→M4, M4→M6 and M2→M6 intervals were treated as secondary/exploratory and were not allowed to replace the prespecified M0→M2 primary test.

## Statistical assumptions and interpretation guardrails

### Linearity and additivity

The primary regression models assume that, over the observed range, the conditional mean of ecological resistance can be approximated as an additive linear function of standardized PC1 and PC2 plus required design covariates. This does not assert that ecological systems are mechanistically linear. It defines the primary low-dimensional association model. Quadratic, spline or other nonlinear models can be scientifically interesting but should be treated as secondary unless specified prospectively.

### Outcome scale

Ecological resistance is bounded approximately between 0 and 1 because it is derived as 1 minus Bray–Curtis dissimilarity. Linear regression can in principle predict outside this range. The primary analyses use linear models because coefficients, cluster-robust inference and joint architecture tests are directly interpretable, and observed responses lie within the feasible range. Residual diagnostics and robustness analyses should be examined; alternative bounded-response models may be evaluated as sensitivity analyses but should not be selected post hoc solely to improve significance.

### Independence and clustering

Independence is required at the level of the sampling cluster represented by the covariance estimator. In Italy, Manitoba and Porcellato, multiple quarter trajectories may arise from the same cow, so cow is treated as the independent cluster. In Van Beeck and endpoint-specific Patangia analyses, one trajectory-level response is modeled per retained cow/core, so trajectory/cow observations are treated as independent after design adjustment.

### Heteroskedasticity

CR2 and HC3 reduce reliance on homoskedastic residual variance. They do not protect against incorrect biological pairing, omitted nonlinearity that materially changes the mean structure, or unmeasured confounding.

### PCA and identifiability

PCA is an unsupervised representation of baseline variation. It is fitted without using ecological-resistance outcomes. PC signs are arbitrary: multiplying a loading vector and all corresponding scores by −1 produces the same PCA solution. Therefore, signs of PC coefficients cannot be compared naively across independently fitted cohorts. Likewise, high-loading families define an axis but are not individually tested causal biomarkers.

### Compositionality and zeros

CLR analysis requires strictly positive components; pseudocount replacement is therefore necessary for zeros. The chosen replacement affects samples with many zeros and is a source of analytic sensitivity. The discovery analysis used a fixed 0.5 pseudocount on its normalized count representation; Van Beeck used half the smallest positive baseline relative abundance. External cohorts were processed independently, so pseudocount details were kept within cohort rather than used to claim exact axis identity. These choices should be reproduced exactly from the frozen scripts for final analysis.

### Low-biomass mammary microbiome limitations

Healthy bovine milk is a low-biomass material, making results potentially sensitive to contamination, extraction/library effects, reagent background and sequencing depth. The present study relies primarily on relative taxonomic profiles and cannot distinguish viable from nonviable organisms unless the source experiment supplied complementary culture information. Cross-study differences in 16S region, platform, laboratory, geography, breed, lactation timing and treatment context are therefore treated as sources of heterogeneity rather than ignored.

### Causality

The architecture–resistance association is temporal but observational. Baseline ecological state precedes the defined response, but this does not establish that manipulating PC1/PC2 or any loading family would causally alter resistance. Host physiology, unmeasured environment and technical factors may influence both baseline composition and longitudinal change. We therefore use language such as “associated with,” “predicts” in the statistical sense of contains prospective information, or “phenomenon-level validation,” rather than causal terms such as “confers resistance” unless supported by future experimental intervention.

## Reproducibility and provenance

All MMER cohort assignments, reconstruction rules, analysis manifests, frozen hypotheses and executable scripts are maintained in the public repository `ajulojays/MMER`. Italy reconstruction is implemented in `scripts/00_build_metadata.py`; the frozen Italy–Manitoba discovery analysis in `scripts/22_primary_healthy_twoAims.R`; the Van Beeck canonical external validation in `scripts/26j_vanbeeck_external_validation_only.R`; and publication visualization for Van Beeck in `scripts/26k_vanbeeck_external_validation_plots.R`. The repository also contains frozen evidence documents recording cohort counts, model definitions and inferential guardrails.

For publication, reproducibility is defined at two levels. **Biological provenance reproducibility** requires that accession→cow→quarter→timepoint→treatment/farm mappings can be regenerated and checked against the source publication. **Analytical reproducibility** requires that the frozen community matrices, prevalence rules, CLR/PCA procedures, outcome definitions and statistical tests reproduce the reported estimates. Where exact upstream processing commands for an older external cohort are not currently committed, the manuscript should report only parameters verified from the source publication and frozen MMER outputs rather than infer unrecorded trimming or denoising settings. A final submission audit should confirm software versions and command-line parameters against the archived execution environment before journal submission.

## Study-specific reconstruction audit summary

| Cohort | Source biological design | MMER longitudinal unit | Reconstruction check | Frozen analysis role |
|---|---|---|---|---|
| Italy / PRJEB38332 | 5 cows × 4 quarters × 3 timepoints; fixed quarter-specific dry-cow conditions | same cow-quarter, T1→T2 | 60/60 runs reconstructed; complete S1–S60; cow/quarter/time/treatment alias checks | Discovery |
| Manitoba | initially 11 cows; healthy quarters; pre-DCT milk and postpartum colostrum; source retained 9 cows | same cow-quarter milk→colostrum | 29 paired milk + 29 paired colostrum trajectories from 9 cows, matching source retained mammary-secretions structure | Discovery |
| Porcellato / PRJEB35792 | 60 Norwegian Red cows, 2 farms, 2 sampling periods, four quarters | same cow-quarter Sampling1→Sampling2 | 189 complete reconstructed trajectories; locked ≥1000 DADA2-read criterion retained 154 (93 A, 61 K) | External 16S validation |
| Van Beeck / PRJEB63336 | 3 dairies; baseline, 7 d, 55–75 DIM; four SCC/treatment groups | complete milk trajectory across all 3 states | locked 60 trajectories/180 samples; exact 3-timepoint completeness; treatment totals 14/12/18/16 | External 16S validation |
| Patangia | 24 cows, front-right quarter, M0/M2/M4/M6; 3 DCT groups; 5 farms | same cow longitudinal front-right trajectory | 96 samples, 24 cows × 4 states; 45,868 family rows; 852 family labels; complete finite matrix | Shotgun external generalization |

## References

1. Biscarini, F. *et al.* A randomized controlled trial of teat-sealant and antibiotic dry-cow treatments for mastitis prevention shows similar effect on the healthy milk microbiome. **Front. Vet. Sci. 7**, 581 (2020). https://doi.org/10.3389/fvets.2020.00581.
2. Derakhshani, H., Plaizier, J. C., De Buck, J., Barkema, H. W. & Khafipour, E. Composition of the teat canal and intramammary microbiota of dairy cows subjected to antimicrobial dry cow therapy and internal teat sealant. **J. Dairy Sci. 101**, 10191–10205 (2018). https://doi.org/10.3168/jds.2018-14858.
3. Porcellato, D. *et al.* A core microbiota dominates a rich microbial diversity in the bovine udder and may indicate presence of dysbiosis. **Sci. Rep. 10**, 21608 (2020). https://doi.org/10.1038/s41598-020-77054-6.
4. Patangia, D. V., Grimaud, G., Linehan, K., Ross, R. P. & Stanton, C. Microbiota and resistome analysis of colostrum and milk from dairy cows treated with and without dry cow therapies. **Antibiotics 12**, 1315 (2023). https://doi.org/10.3390/antibiotics12081315.
5. Van Beeck, W. *et al.* Variations in cow milk and teat skin microbiota across the lactation cycle with intramammary cephalosporin use at dry-off. **Appl. Environ. Microbiol. 92**, e02312-25 (2026). https://doi.org/10.1128/aem.02312-25.
6. Bray, J. R. & Curtis, J. T. An ordination of the upland forest communities of southern Wisconsin. **Ecol. Monogr. 27**, 325–349 (1957).
7. Aitchison, J. **The Statistical Analysis of Compositional Data** (Chapman & Hall, 1986).
8. Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V. & Egozcue, J. J. Microbiome datasets are compositional: and this is not optional. **Front. Microbiol. 8**, 2224 (2017). https://doi.org/10.3389/fmicb.2017.02224.
9. Callahan, B. J. *et al.* DADA2: high-resolution sample inference from Illumina amplicon data. **Nat. Methods 13**, 581–583 (2016). https://doi.org/10.1038/nmeth.3869.
10. Quast, C. *et al.* The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. **Nucleic Acids Res. 41**, D590–D596 (2013). https://doi.org/10.1093/nar/gks1219.
11. Wood, D. E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. **Genome Biol. 20**, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0.
12. Pustejovsky, J. E. & Tipton, E. Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. **J. Bus. Econ. Stat. 36**, 672–683 (2018). https://doi.org/10.1080/07350015.2016.1247004.
13. Freedman, D. & Lane, D. A nonstochastic interpretation of reported significance levels. **J. Bus. Econ. Stat. 1**, 292–298 (1983).
14. MacKinnon, J. G. & White, H. Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. **J. Econometrics 29**, 305–325 (1985).
15. Oksanen, J. *et al.* vegan: Community Ecology Package. R package (version used recorded in the frozen analysis environment).
