#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dada2)
  library(readr)
  library(dplyr)
})

root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))] |> sub("--file=", "", x = _)), ".."), mustWork = TRUE)
meta_path <- file.path(root, "data/metadata/PRJEB38332_sample_map.csv")
raw_dir <- file.path(root, "data/raw")
out_dir <- file.path(root, "results/asv")
qc_dir <- file.path(root, "results/qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(meta_path, show_col_types = FALSE)
fnFs <- file.path(raw_dir, paste0(meta$run_accession, "_R1.fastq.gz"))
fnRs <- file.path(raw_dir, paste0(meta$run_accession, "_R2.fastq.gz"))
names(fnFs) <- meta$run_accession
names(fnRs) <- meta$run_accession
stopifnot(all(file.exists(fnFs)), all(file.exists(fnRs)))

# IMPORTANT: these truncation values are conservative starting points only.
# Inspect quality profiles before final production analysis.
filt_dir <- file.path(out_dir, "filtered")
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
filtFs <- file.path(filt_dir, paste0(names(fnFs), "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(names(fnRs), "_R_filt.fastq.gz"))

filt <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      truncLen = c(240, 200), maxN = 0,
                      maxEE = c(2, 2), truncQ = 2,
                      rm.phix = TRUE, compress = TRUE,
                      multithread = TRUE)
write.csv(filt, file.path(qc_dir, "dada2_filtering.csv"))

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
derepFs <- derepFastq(filtFs); names(derepFs) <- names(fnFs)
derepRs <- derepFastq(filtRs); names(derepRs) <- names(fnRs)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 20, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

saveRDS(seqtab.nochim, file.path(out_dir, "asv_table.rds"))
write.csv(seqtab.nochim, file.path(out_dir, "asv_table.csv"))

track <- cbind(filt,
  denoisedF = sapply(dadaFs, getUniques) |> sapply(sum),
  merged = sapply(mergers, getUniques) |> sapply(sum),
  nonchim = rowSums(seqtab.nochim)
)
write.csv(track, file.path(qc_dir, "dada2_read_tracking.csv"))
cat("[OK] DADA2 ASV table:", nrow(seqtab.nochim), "samples x", ncol(seqtab.nochim), "ASVs\n")
