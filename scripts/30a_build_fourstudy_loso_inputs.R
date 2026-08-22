#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

ROOT <- "/home/samuelajulo/MeteG/MMER_repo"
OUT <- file.path(ROOT, "results/fourstudy_loso/input")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

EXPECTED <- c(
  Italy=20L,
  Manitoba=29L,
  Porcellato=154L,
  VanBeeck=60L
)

clean_family <- function(x) {
  x <- trimws(as.character(x))
  x[x==""] <- NA_character_
  x
}

remove_organelle_labels <- function(d) {
  d %>% filter(
    !grepl("mitochond|chloroplast|plastid|organelle", family, ignore.case=TRUE)
  )
}

normalize_within_trajectory <- function(d) {
  d %>%
    group_by(study, trajectory_id, cow_id) %>%
    mutate(
      total_abundance=sum(abundance, na.rm=TRUE),
      abundance=ifelse(total_abundance>0, abundance/total_abundance, NA_real_)
    ) %>%
    ungroup() %>%
    select(-total_abundance) %>%
    filter(is.finite(abundance), abundance>=0)
}

cat("\n====================================\n")
cat("BUILD FOUR-STUDY LOSO INPUTS\n")
cat("====================================\n")

# =============================================================================
# 1. ITALY + MANITOBA
#    Source: existing harmonized family long table.
#    Baseline is harmonized T1. Only the frozen 49 primary discovery trajectories
#    (20 Italy + 29 Manitoba) are retained.
# =============================================================================

HM_FILE <- file.path(ROOT, "results/harmonized/three_study_family_long.csv")
DISC_OUT <- file.path(ROOT, "results/primary_healthy_twoAims/tables/quarter_ecological_resistance.csv")

stopifnot(file.exists(HM_FILE), file.exists(DISC_OUT))

hm <- read_csv(HM_FILE, show_col_types=FALSE)
disc_out0 <- read_csv(DISC_OUT, show_col_types=FALSE)

stopifnot(
  all(c("study","trajectory_global","cow_global","harmonized_timepoint","family","abundance") %in% names(hm)),
  all(c("study","trajectory_global","cow_global","ecological_resistance") %in% names(disc_out0))
)

disc_out <- disc_out0 %>%
  filter(study %in% c("Italy","Manitoba")) %>%
  transmute(
    study=as.character(study),
    trajectory_id=as.character(trajectory_global),
    cow_id=as.character(cow_global),
    ecological_resistance=as.numeric(ecological_resistance)
  )

stopifnot(
  sum(disc_out$study=="Italy")==EXPECTED[["Italy"]],
  sum(disc_out$study=="Manitoba")==EXPECTED[["Manitoba"]],
  !anyDuplicated(disc_out[,c("study","trajectory_id")])
)

im_family <- hm %>%
  filter(
    study %in% c("Italy","Manitoba"),
    harmonized_timepoint=="T1",
    trajectory_global %in% disc_out$trajectory_id
  ) %>%
  transmute(
    study=as.character(study),
    trajectory_id=as.character(trajectory_global),
    cow_id=as.character(cow_global),
    family=clean_family(family),
    abundance=as.numeric(abundance)
  ) %>%
  filter(!is.na(family), is.finite(abundance), abundance>=0) %>%
  group_by(study,trajectory_id,cow_id,family) %>%
  summarise(abundance=sum(abundance),.groups="drop") %>%
  remove_organelle_labels() %>%
  normalize_within_trajectory()

im_ids <- im_family %>% distinct(study,trajectory_id,cow_id)
stopifnot(
  nrow(im_ids %>% filter(study=="Italy"))==EXPECTED[["Italy"]],
  nrow(im_ids %>% filter(study=="Manitoba"))==EXPECTED[["Manitoba"]],
  nrow(anti_join(disc_out,im_ids,by=c("study","trajectory_id","cow_id")))==0
)

cat("Italy baseline trajectories:",sum(im_ids$study=="Italy"),"\n")
cat("Manitoba baseline trajectories:",sum(im_ids$study=="Manitoba"),"\n")

# =============================================================================
# 2. PORCELLATO
#    Frozen primary external cohort: 154 trajectories.
#    Baseline sample is Sampling1, deterministically reconstructed from
#    farm/cow/quarter: Cow{cow}_{quarter}_Sampling1_{farm}.
# =============================================================================

P_FAM <- file.path(ROOT,"results/study7_porcellato/taxonomy/family_counts_long.csv")
P_OUT <- file.path(ROOT,"results/study7_porcellato/external_resistance/ecological_resistance_primary154.csv")
P_QC  <- file.path(ROOT,"results/study7_porcellato/trajectory_depth_QC_locked.csv")

stopifnot(file.exists(P_FAM),file.exists(P_OUT),file.exists(P_QC))

p_counts <- read_csv(P_FAM,show_col_types=FALSE)
p_out0 <- read_csv(P_OUT,show_col_types=FALSE)
p_qc <- read_csv(P_QC,show_col_types=FALSE)

stopifnot(
  all(c("sample","family_final","count") %in% names(p_counts)),
  all(c("trajectory_id","cow","quarter","farm","ecological_resistance") %in% names(p_out0)),
  all(c("trajectory_id","primary_keep") %in% names(p_qc))
)

p_out <- p_out0 %>%
  inner_join(p_qc %>% select(trajectory_id,primary_keep),by="trajectory_id") %>%
  filter(primary_keep) %>%
  mutate(
    baseline_sample=paste0("Cow",cow,"_",quarter,"_Sampling1_",farm),
    study="Porcellato",
    cow_id=paste0("Porcellato__",farm,"_",cow)
  ) %>%
  transmute(
    study,
    trajectory_id=as.character(trajectory_id),
    cow_id=as.character(cow_id),
    baseline_sample,
    ecological_resistance=as.numeric(ecological_resistance)
  )

stopifnot(
  nrow(p_out)==EXPECTED[["Porcellato"]],
  !anyDuplicated(p_out$trajectory_id)
)

p_family <- p_counts %>%
  transmute(
    baseline_sample=as.character(sample),
    family=clean_family(family_final),
    abundance=as.numeric(count)
  ) %>%
  filter(!is.na(family),is.finite(abundance),abundance>=0) %>%
  inner_join(p_out %>% select(study,trajectory_id,cow_id,baseline_sample),by="baseline_sample") %>%
  select(study,trajectory_id,cow_id,family,abundance) %>%
  group_by(study,trajectory_id,cow_id,family) %>%
  summarise(abundance=sum(abundance),.groups="drop") %>%
  remove_organelle_labels() %>%
  normalize_within_trajectory()

p_ids <- p_family %>% distinct(study,trajectory_id,cow_id)
stopifnot(
  nrow(p_ids)==EXPECTED[["Porcellato"]],
  nrow(anti_join(p_out,p_ids,by=c("study","trajectory_id","cow_id")))==0
)

cat("Porcellato baseline trajectories:",nrow(p_ids),"\n")

# =============================================================================
# 3. VAN BEECK
#    Frozen all-60 cohort. Baseline sample identifiers and resistance outcome are
#    taken from family_ecological_resistance_60.csv. Baseline family composition
#    is reconstructed from the locked rarefied 3k bacterial ASV table and the
#    bacterial taxonomy object, using the same family fallback hierarchy used in
#    the Van Beeck family-level analyses.
#
#    IMPORTANT: the available locked files do not expose a separate animal ID.
#    Therefore core/trajectory is used as the clustering ID. This is explicit in
#    the QC output and should be reported as a limitation/sensitivity issue.
# =============================================================================

VB_OUT <- file.path(ROOT,"results/vanbeeck/all60_ecological_resistance/family_ecological_resistance_60.csv")
VB_SEQ <- file.path(ROOT,"results/vanbeeck/analysis_3k/seqtab_bacterial_rarefied_3000.rds")
VB_TAX <- file.path(ROOT,"results/vanbeeck/taxonomy/taxonomy_bacterial.rds")

stopifnot(file.exists(VB_OUT),file.exists(VB_SEQ),file.exists(VB_TAX))

vb0 <- read_csv(VB_OUT,show_col_types=FALSE)
seqtab <- readRDS(VB_SEQ)
tax <- readRDS(VB_TAX)

stopifnot(
  all(c("core","baseline_sample","overall_resistance") %in% names(vb0)),
  nrow(vb0)==EXPECTED[["VanBeeck"]]
)

vb_out <- vb0 %>%
  transmute(
    study="VanBeeck",
    trajectory_id=as.character(core),
    cow_id=as.character(core),
    baseline_sample=as.character(baseline_sample),
    ecological_resistance=as.numeric(overall_resistance)
  )

# Ensure orientation is sample x ASV.
if (!all(vb_out$baseline_sample %in% rownames(seqtab))) {
  if (all(vb_out$baseline_sample %in% colnames(seqtab))) {
    seqtab <- t(seqtab)
  } else {
    missing <- setdiff(vb_out$baseline_sample,rownames(seqtab))
    stop("Van Beeck baseline samples missing from seqtab: ",paste(head(missing,10),collapse=", "))
  }
}

# Align taxonomy to ASV columns.
if (is.null(rownames(tax))) stop("Van Beeck taxonomy object has no rownames")
common_asv <- intersect(colnames(seqtab),rownames(tax))
if (length(common_asv)==0) stop("No shared ASV IDs between Van Beeck seqtab and taxonomy")
seqtab <- seqtab[,common_asv,drop=FALSE]
tax <- tax[common_asv,,drop=FALSE]

# Remove organellar labels defensively.
tax_text <- apply(tax,1,function(z) paste(z,collapse=" "))
org <- grepl("mitochond|chloroplast|plastid|organelle",tax_text,ignore.case=TRUE)
seqtab <- seqtab[,!org,drop=FALSE]
tax <- tax[!org,,drop=FALSE]

get_tax_col <- function(mat,nm) {
  if (nm %in% colnames(mat)) as.character(mat[,nm]) else rep(NA_character_,nrow(mat))
}

fam_lab <- clean_family(get_tax_col(tax,"Family"))
ord <- clean_family(get_tax_col(tax,"Order"))
cls <- clean_family(get_tax_col(tax,"Class"))
phy <- clean_family(get_tax_col(tax,"Phylum"))

for (i in seq_along(fam_lab)) {
  if (is.na(fam_lab[i])) {
    if (!is.na(ord[i])) fam_lab[i] <- paste0("Unclassified_",ord[i])
    else if (!is.na(cls[i])) fam_lab[i] <- paste0("Unclassified_",cls[i])
    else if (!is.na(phy[i])) fam_lab[i] <- paste0("Unclassified_",phy[i])
    else fam_lab[i] <- "Unclassified_Bacteria"
  }
}

# Aggregate sample x ASV -> sample x family.
vb_mat <- t(rowsum(t(as.matrix(seqtab)),group=fam_lab,reorder=FALSE))
vb_mat <- vb_mat[vb_out$baseline_sample,,drop=FALSE]

vb_long <- as.data.frame(vb_mat) %>%
  rownames_to_column("baseline_sample") %>%
  pivot_longer(-baseline_sample,names_to="family",values_to="abundance") %>%
  filter(abundance>0) %>%
  inner_join(vb_out %>% select(study,trajectory_id,cow_id,baseline_sample),by="baseline_sample") %>%
  select(study,trajectory_id,cow_id,family,abundance) %>%
  group_by(study,trajectory_id,cow_id,family) %>%
  summarise(abundance=sum(abundance),.groups="drop") %>%
  remove_organelle_labels() %>%
  normalize_within_trajectory()

vb_ids <- vb_long %>% distinct(study,trajectory_id,cow_id)
stopifnot(
  nrow(vb_ids)==EXPECTED[["VanBeeck"]],
  nrow(anti_join(vb_out,vb_ids,by=c("study","trajectory_id","cow_id")))==0
)

cat("Van Beeck baseline trajectories:",nrow(vb_ids),"\n")

# =============================================================================
# 4. COMBINE + HARD QC
# =============================================================================

baseline_family_long <- bind_rows(im_family,p_family,vb_long) %>%
  mutate(
    study=factor(study,levels=names(EXPECTED)),
    family=clean_family(family)
  ) %>%
  arrange(study,trajectory_id,family) %>%
  mutate(study=as.character(study))

outcomes <- bind_rows(
  disc_out,
  p_out %>% select(study,trajectory_id,cow_id,ecological_resistance),
  vb_out %>% select(study,trajectory_id,cow_id,ecological_resistance)
) %>%
  mutate(study=factor(study,levels=names(EXPECTED))) %>%
  arrange(study,trajectory_id) %>%
  mutate(study=as.character(study))

# Expected study sizes.
counts <- outcomes %>% count(study,name="trajectories")
for (s in names(EXPECTED)) {
  got <- counts$trajectories[counts$study==s]
  if (length(got)!=1 || got!=EXPECTED[[s]]) {
    stop("Unexpected trajectory count for ",s,": got ",paste(got,collapse=",")," expected ",EXPECTED[[s]])
  }
}
stopifnot(nrow(outcomes)==sum(EXPECTED))
stopifnot(!anyDuplicated(outcomes[,c("study","trajectory_id")]))
stopifnot(!anyNA(outcomes$ecological_resistance),all(is.finite(outcomes$ecological_resistance)))

# Exact trajectory agreement between family and outcome tables.
f_ids <- baseline_family_long %>% distinct(study,trajectory_id,cow_id)
o_ids <- outcomes %>% distinct(study,trajectory_id,cow_id)
miss_family <- anti_join(o_ids,f_ids,by=c("study","trajectory_id","cow_id"))
extra_family <- anti_join(f_ids,o_ids,by=c("study","trajectory_id","cow_id"))
if (nrow(miss_family)>0 || nrow(extra_family)>0) {
  print(miss_family)
  print(extra_family)
  stop("Family/outcome trajectory mismatch")
}

# Row-sum QC should be exactly ~1 after within-trajectory normalization.
rowsums <- baseline_family_long %>%
  group_by(study,trajectory_id,cow_id) %>%
  summarise(sum_abundance=sum(abundance),families_positive=sum(abundance>0),.groups="drop")

if (any(abs(rowsums$sum_abundance-1)>1e-8)) {
  stop("At least one baseline trajectory does not sum to 1 after normalization")
}

input_qc <- outcomes %>%
  group_by(study) %>%
  summarise(
    trajectories=n(),
    cows_or_clusters=n_distinct(cow_id),
    mean_resistance=mean(ecological_resistance),
    sd_resistance=sd(ecological_resistance),
    min_resistance=min(ecological_resistance),
    max_resistance=max(ecological_resistance),
    .groups="drop"
  ) %>%
  left_join(
    baseline_family_long %>%
      group_by(study) %>%
      summarise(
        families_observed=n_distinct(family),
        family_rows=n(),
        .groups="drop"
      ),by="study"
  ) %>%
  left_join(
    rowsums %>%
      group_by(study) %>%
      summarise(
        min_positive_families=min(families_positive),
        median_positive_families=median(families_positive),
        max_positive_families=max(families_positive),
        min_row_sum=min(sum_abundance),
        max_row_sum=max(sum_abundance),
        .groups="drop"
      ),by="study"
  ) %>%
  mutate(
    cluster_note=case_when(
      study=="VanBeeck" ~ "No separate animal ID in locked source table; trajectory core used as cluster ID",
      TRUE ~ "Animal/cow identifier available"
    )
  )

family_harmonization_qc <- baseline_family_long %>%
  distinct(study,family) %>%
  mutate(present=1L) %>%
  pivot_wider(names_from=study,values_from=present,values_fill=0L) %>%
  mutate(
    n_studies_present=Italy+Manitoba+Porcellato+VanBeeck
  ) %>%
  arrange(desc(n_studies_present),family)

pairwise_overlap <- expand_grid(
  study_A=names(EXPECTED),
  study_B=names(EXPECTED)
) %>%
  rowwise() %>%
  mutate(
    families_A=list(unique(baseline_family_long$family[baseline_family_long$study==study_A])),
    families_B=list(unique(baseline_family_long$family[baseline_family_long$study==study_B])),
    n_A=length(families_A[[1]]),
    n_B=length(families_B[[1]]),
    intersection=length(intersect(families_A[[1]],families_B[[1]])),
    union=length(union(families_A[[1]],families_B[[1]])),
    jaccard=intersection/union
  ) %>%
  ungroup() %>%
  select(study_A,study_B,n_A,n_B,intersection,union,jaccard)

# =============================================================================
# 5. WRITE FROZEN STANDARDIZED INPUTS + QC
# =============================================================================

write_csv(baseline_family_long,file.path(OUT,"baseline_family_long.csv"))
write_csv(outcomes,file.path(OUT,"outcomes.csv"))
write_csv(input_qc,file.path(OUT,"input_QC_by_study.csv"))
write_csv(family_harmonization_qc,file.path(OUT,"family_harmonization_QC.csv"))
write_csv(pairwise_overlap,file.path(OUT,"family_pairwise_overlap_QC.csv"))
write_csv(rowsums,file.path(OUT,"baseline_row_sum_QC.csv"))

capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo_builder.txt"))

cat("\n====================================\n")
cat("FOUR-STUDY STANDARDIZED INPUT QC\n")
cat("====================================\n")
print(input_qc,width=Inf)

cat("\nFamilies observed in all 4 studies:",sum(family_harmonization_qc$n_studies_present==4),"\n")
cat("Families observed in >=3 studies:",sum(family_harmonization_qc$n_studies_present>=3),"\n")
cat("Families observed in >=2 studies:",sum(family_harmonization_qc$n_studies_present>=2),"\n")
cat("Total standardized trajectories:",nrow(outcomes),"\n")

cat("\n====================================\n")
cat("[PASS] FOUR-STUDY LOSO INPUT BUILD COMPLETE\n")
cat("====================================\n")
cat("baseline_family_long.csv:",file.path(OUT,"baseline_family_long.csv"),"\n")
cat("outcomes.csv:",file.path(OUT,"outcomes.csv"),"\n")
cat("Expected total trajectories:",sum(EXPECTED),"\n")
cat("Observed total trajectories:",nrow(outcomes),"\n")
