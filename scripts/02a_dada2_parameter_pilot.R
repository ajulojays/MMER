#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ============================================================
# MMER — DADA2 parameter pilot
#
# Goal:
# Compare alternative filtering/truncation strategies using
# real retention through:
#
# input -> filtered -> denoised -> merged -> non-chimeric
#
# Production DADA2 parameters should NOT be chosen until this
# pilot has been evaluated.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

ROOT <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

META <- file.path(
  ROOT,
  "data",
  "metadata",
  "PRJEB38332_sample_map.csv"
)

RAW_DIR <- file.path(ROOT, "data", "primer_trimmed")

OUT_ROOT <- file.path(
  ROOT,
  "results",
  "dada2_parameter_pilot_primer_trimmed"
)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("============================================================\n")
cat(" MMER DADA2 PARAMETER PILOT\n")
cat("============================================================\n")
cat("Project root :", ROOT, "\n")
cat("Metadata     :", META, "\n")
cat("Raw reads    :", RAW_DIR, "\n")
cat("Output       :", OUT_ROOT, "\n\n")

# ------------------------------------------------------------
# 1. Load and validate metadata
# ------------------------------------------------------------

if (!file.exists(META)) {
  stop("Metadata file not found: ", META)
}

meta <- read_csv(
  META,
  show_col_types = FALSE,
  col_types = cols(
    cow_id = col_character()
  )
)

required_cols <- c(
  "run_accession",
  "cow_id",
  "quarter",
  "treatment",
  "timepoint"
)

missing_cols <- setdiff(required_cols, names(meta))

if (length(missing_cols) > 0) {
  stop(
    "Metadata missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

meta <- meta %>%
  mutate(
    R1 = file.path(
      RAW_DIR,
      paste0(run_accession, "_R1.fastq.gz")
    ),
    R2 = file.path(
      RAW_DIR,
      paste0(run_accession, "_R2.fastq.gz")
    )
  )

missing_r1 <- meta$R1[!file.exists(meta$R1)]
missing_r2 <- meta$R2[!file.exists(meta$R2)]

if (length(missing_r1) > 0 || length(missing_r2) > 0) {
  stop(
    "FASTQ files are missing. R1 missing=",
    length(missing_r1),
    "; R2 missing=",
    length(missing_r2)
  )
}

cat("[OK] All 60 paired FASTQ samples located.\n")

# ------------------------------------------------------------
# 2. Select representative pilot samples
#
# First take one representative sample from every:
# treatment x timepoint combination = 12 samples
#
# Then add additional samples deterministically until N=20.
#
# This keeps the pilot reasonably fast while spanning the
# experimental design.
# ------------------------------------------------------------

PILOT_N <- 20

core <- meta %>%
  arrange(
    treatment,
    timepoint,
    cow_id,
    quarter
  ) %>%
  group_by(
    treatment,
    timepoint
  ) %>%
  slice(1) %>%
  ungroup()

remaining <- meta %>%
  filter(
    !run_accession %in% core$run_accession
  ) %>%
  arrange(
    cow_id,
    treatment,
    timepoint
  )

pilot <- bind_rows(
  core,
  head(
    remaining,
    max(0, PILOT_N - nrow(core))
  )
) %>%
  distinct(run_accession, .keep_all = TRUE) %>%
  head(PILOT_N)

write_csv(
  pilot %>%
    select(
      run_accession,
      cow_id,
      quarter,
      treatment,
      timepoint
    ),
  file.path(
    OUT_ROOT,
    "pilot_samples.csv"
  )
)

cat(
  "[OK] Pilot samples selected:",
  nrow(pilot),
  "\n"
)

print(
  pilot %>%
    count(
      treatment,
      timepoint
    )
)

# ------------------------------------------------------------
# 3. Parameter grid
#
# R2 quality deteriorates strongly after ~175-200 bp, BUT this
# is V3-V4 and aggressive truncation can destroy paired-end
# overlap. Therefore we explicitly test both fixed-length and
# no-fixed-truncation strategies.
# ------------------------------------------------------------

PARAMETERS <- tribble(
  ~strategy, ~truncF, ~truncR, ~maxEEF, ~maxEER,
  "A_225_205_EE2_3", 225, 205, 2, 3,
  "B_225_200_EE2_3", 225, 200, 2, 3,
  "C_230_200_EE2_3", 230, 200, 2, 3,
  "D_noTrunc_EE2_3",   0,   0, 2, 3,
  "E_225_200_EE2_5", 225, 200, 2, 5
)

write_csv(
  PARAMETERS,
  file.path(
    OUT_ROOT,
    "parameter_grid.csv"
  )
)

print(PARAMETERS)

# ------------------------------------------------------------
# Utility
# ------------------------------------------------------------

getN <- function(x) {
  sum(getUniques(x))
}

safe_pct <- function(num, den) {
  ifelse(
    den > 0,
    100 * num / den,
    NA_real_
  )
}

all_strategy_summaries <- list()

# ------------------------------------------------------------
# 4. Run strategies
# ------------------------------------------------------------

for (i in seq_len(nrow(PARAMETERS))) {

  p <- PARAMETERS[i, ]

  strategy <- p$strategy

  cat("\n")
  cat("============================================================\n")
  cat(" STRATEGY:", strategy, "\n")
  cat("============================================================\n")

  cat(
    "truncLen = c(",
    p$truncF,
    ", ",
    p$truncR,
    ")\n",
    sep = ""
  )

  cat(
    "maxEE    = c(",
    p$maxEEF,
    ", ",
    p$maxEER,
    ")\n",
    sep = ""
  )

  strategy_dir <- file.path(
    OUT_ROOT,
    strategy
  )

  filtered_dir <- file.path(
    strategy_dir,
    "filtered"
  )

  dir.create(
    filtered_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  sample_ids <- pilot$run_accession

  fnFs <- pilot$R1
  fnRs <- pilot$R2

  filtFs <- file.path(
    filtered_dir,
    paste0(
      sample_ids,
      "_F_filt.fastq.gz"
    )
  )

  filtRs <- file.path(
    filtered_dir,
    paste0(
      sample_ids,
      "_R_filt.fastq.gz"
    )
  )

  names(fnFs) <- sample_ids
  names(fnRs) <- sample_ids
  names(filtFs) <- sample_ids
  names(filtRs) <- sample_ids

  # ----------------------------------------------------------
  # Filtering
  # ----------------------------------------------------------

  cat("[STEP] filterAndTrim\n")

  filt_out <- filterAndTrim(
    fnFs,
    filtFs,
    fnRs,
    filtRs,

    truncLen = c(
      p$truncF,
      p$truncR
    ),

    maxN = 0,

    maxEE = c(
      p$maxEEF,
      p$maxEER
    ),

    truncQ = 2,

    rm.phix = TRUE,

    compress = TRUE,

    multithread = TRUE,

    verbose = TRUE
  )

  # ----------------------------------------------------------
  # Learn error models
  # ----------------------------------------------------------

  cat("[STEP] learnErrors forward\n")

  errF <- learnErrors(
    filtFs,
    multithread = TRUE,
    randomize = TRUE,
    verbose = TRUE
  )

  cat("[STEP] learnErrors reverse\n")

  errR <- learnErrors(
    filtRs,
    multithread = TRUE,
    randomize = TRUE,
    verbose = TRUE
  )

  saveRDS(
    errF,
    file.path(
      strategy_dir,
      "error_model_forward.rds"
    )
  )

  saveRDS(
    errR,
    file.path(
      strategy_dir,
      "error_model_reverse.rds"
    )
  )

  # ----------------------------------------------------------
  # Dereplicate
  # ----------------------------------------------------------

  cat("[STEP] dereplication\n")

  derepFs <- derepFastq(
    filtFs,
    verbose = FALSE
  )

  derepRs <- derepFastq(
    filtRs,
    verbose = FALSE
  )

  names(derepFs) <- sample_ids
  names(derepRs) <- sample_ids

  # ----------------------------------------------------------
  # Denoising
  # ----------------------------------------------------------

  cat("[STEP] denoise forward\n")

  dadaFs <- dada(
    derepFs,
    err = errF,
    multithread = TRUE,
    pool = FALSE
  )

  cat("[STEP] denoise reverse\n")

  dadaRs <- dada(
    derepRs,
    err = errR,
    multithread = TRUE,
    pool = FALSE
  )

  # ----------------------------------------------------------
  # Merge paired reads
  # ----------------------------------------------------------

  cat("[STEP] mergePairs\n")

  mergers <- mergePairs(
    dadaFs,
    derepFs,
    dadaRs,
    derepRs,
    minOverlap = 12,
    maxMismatch = 0,
    verbose = TRUE
  )

  # ----------------------------------------------------------
  # Sequence table
  # ----------------------------------------------------------

  cat("[STEP] makeSequenceTable\n")

  seqtab <- makeSequenceTable(
    mergers
  )

  # ----------------------------------------------------------
  # Chimera removal
  # ----------------------------------------------------------

  if (
    nrow(seqtab) > 0 &&
    ncol(seqtab) > 0
  ) {

    cat("[STEP] removeBimeraDenovo\n")

    seqtab_nochim <- removeBimeraDenovo(
      seqtab,
      method = "consensus",
      multithread = TRUE,
      verbose = TRUE
    )

  } else {

    warning(
      strategy,
      ": no merged sequences available."
    )

    seqtab_nochim <- seqtab
  }

  saveRDS(
    seqtab,
    file.path(
      strategy_dir,
      "seqtab_merged.rds"
    )
  )

  saveRDS(
    seqtab_nochim,
    file.path(
      strategy_dir,
      "seqtab_nochim.rds"
    )
  )

  # ----------------------------------------------------------
  # Tracking
  # ----------------------------------------------------------

  input_reads <- filt_out[, 1]
  filtered_reads <- filt_out[, 2]

  denoised_F <- sapply(
    dadaFs,
    getN
  )

  denoised_R <- sapply(
    dadaRs,
    getN
  )

  merged_reads <- sapply(
    mergers,
    getN
  )

  if (
    nrow(seqtab_nochim) > 0 &&
    ncol(seqtab_nochim) > 0
  ) {

    nonchim_reads <- rowSums(
      seqtab_nochim
    )

    nonchim_reads <- nonchim_reads[
      sample_ids
    ]

  } else {

    nonchim_reads <- rep(
      0,
      length(sample_ids)
    )

    names(nonchim_reads) <- sample_ids
  }

  track <- tibble(
    run_accession = sample_ids,

    input = as.numeric(
      input_reads
    ),

    filtered = as.numeric(
      filtered_reads
    ),

    denoised_F = as.numeric(
      denoised_F[sample_ids]
    ),

    denoised_R = as.numeric(
      denoised_R[sample_ids]
    ),

    merged = as.numeric(
      merged_reads[sample_ids]
    ),

    nonchim = as.numeric(
      nonchim_reads[sample_ids]
    )
  ) %>%

    left_join(
      pilot %>%
        select(
          run_accession,
          cow_id,
          quarter,
          treatment,
          timepoint
        ),
      by = "run_accession"
    ) %>%

    mutate(
      filtered_pct = safe_pct(
        filtered,
        input
      ),

      merged_pct_input = safe_pct(
        merged,
        input
      ),

      merged_pct_filtered = safe_pct(
        merged,
        filtered
      ),

      nonchim_pct_input = safe_pct(
        nonchim,
        input
      ),

      nonchim_pct_merged = safe_pct(
        nonchim,
        merged
      )
    )

  write_csv(
    track,
    file.path(
      strategy_dir,
      "read_tracking.csv"
    )
  )

  # ----------------------------------------------------------
  # ASV length diagnostics
  # ----------------------------------------------------------

  if (ncol(seqtab_nochim) > 0) {

    asv_lengths <- nchar(
      colnames(
        seqtab_nochim
      )
    )

  } else {

    asv_lengths <- numeric(0)
  }

  if (length(asv_lengths) > 0) {

    length_summary <- tibble(
      min_asv_length = min(
        asv_lengths
      ),

      median_asv_length = median(
        asv_lengths
      ),

      mean_asv_length = mean(
        asv_lengths
      ),

      max_asv_length = max(
        asv_lengths
      )
    )

  } else {

    length_summary <- tibble(
      min_asv_length = NA_real_,
      median_asv_length = NA_real_,
      mean_asv_length = NA_real_,
      max_asv_length = NA_real_
    )
  }

  write_csv(
    tibble(
      ASV_length = asv_lengths
    ),
    file.path(
      strategy_dir,
      "asv_lengths.csv"
    )
  )

  # ----------------------------------------------------------
  # Strategy-level summary
  # ----------------------------------------------------------

  summary_row <- tibble(

    strategy = strategy,

    truncF = p$truncF,
    truncR = p$truncR,

    maxEEF = p$maxEEF,
    maxEER = p$maxEER,

    n_samples = nrow(track),

    total_input = sum(
      track$input
    ),

    total_filtered = sum(
      track$filtered
    ),

    total_merged = sum(
      track$merged
    ),

    total_nonchim = sum(
      track$nonchim
    ),

    filtered_pct = safe_pct(
      sum(track$filtered),
      sum(track$input)
    ),

    merged_pct_input = safe_pct(
      sum(track$merged),
      sum(track$input)
    ),

    merged_pct_filtered = safe_pct(
      sum(track$merged),
      sum(track$filtered)
    ),

    nonchim_pct_input = safe_pct(
      sum(track$nonchim),
      sum(track$input)
    ),

    nonchim_pct_merged = safe_pct(
      sum(track$nonchim),
      sum(track$merged)
    ),

    n_ASVs_merged = ncol(
      seqtab
    ),

    n_ASVs_nonchim = ncol(
      seqtab_nochim
    )
  ) %>%
    bind_cols(
      length_summary
    )

  write_csv(
    summary_row,
    file.path(
      strategy_dir,
      "strategy_summary.csv"
    )
  )

  print(summary_row)

  all_strategy_summaries[[
    strategy
  ]] <- summary_row
}

# ------------------------------------------------------------
# 5. Combined comparison
# ------------------------------------------------------------

comparison <- bind_rows(
  all_strategy_summaries
) %>%
  arrange(
    desc(
      nonchim_pct_input
    )
  )

write_csv(
  comparison,
  file.path(
    OUT_ROOT,
    "parameter_comparison.csv"
  )
)

cat("\n")
cat("============================================================\n")
cat(" PILOT COMPLETE\n")
cat("============================================================\n\n")

print(
  comparison %>%
    select(
      strategy,
      filtered_pct,
      merged_pct_filtered,
      nonchim_pct_input,
      n_ASVs_nonchim,
      median_asv_length
    )
)

cat("\nPrimary comparison written to:\n")
cat(
  file.path(
    OUT_ROOT,
    "parameter_comparison.csv"
  ),
  "\n"
)

cat("\nPer-strategy tracking available under:\n")
cat(
  OUT_ROOT,
  "\n"
)

cat("\n")
cat("IMPORTANT:\n")
cat(
  "Do not choose the winning strategy solely from filtering retention.\n"
)
cat(
  "Prioritize successful merging, non-chimeric retention, and biologically\n"
)
cat(
  "plausible V3-V4 ASV-length distributions.\n"
)

cat("\n[DONE] MMER DADA2 parameter pilot complete.\n")
