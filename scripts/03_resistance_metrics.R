#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(zCompositions)
  library(compositions)
})

args <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args[grep("--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
meta <- read_csv(file.path(root, "data/metadata/PRJEB38332_sample_map.csv"), show_col_types = FALSE)
seqtab <- readRDS(file.path(root, "results/asv/asv_table.rds"))
out_dir <- file.path(root, "results/resistance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- as.matrix(seqtab)
stopifnot(all(meta$run_accession %in% rownames(counts)))
counts <- counts[meta$run_accession, , drop = FALSE]
rel <- counts / rowSums(counts)

# Bray-Curtis is bounded [0,1]. Primary interpretation: lower displacement = higher resistance.
bray <- as.matrix(vegdist(rel, method = "bray"))

# Aitchison sensitivity analysis using multiplicative zero replacement.
# cmultRepl expects samples x components.
clr_mat <- clr(cmultRepl(counts, label = 0, method = "CZM", output = "prop"))
ait <- as.matrix(dist(clr_mat))

alpha <- tibble(
  run_accession = rownames(counts),
  richness = specnumber(counts),
  shannon = diversity(counts, index = "shannon"),
  simpson = diversity(counts, index = "simpson"),
  dominance = apply(rel, 1, max)
)

baseline <- meta |> filter(timepoint == "T1") |>
  select(cow_id, quarter, baseline_run = run_accession)
post <- meta |> filter(timepoint %in% c("T2", "T3")) |>
  select(cow_id, quarter, treatment, timepoint, phase, post_run = run_accession)
pairs <- left_join(post, baseline, by = c("cow_id", "quarter"))

metric_for_pair <- function(base, post) {
  bvec <- rel[base, ]; pvec <- rel[post, ]
  present_b <- bvec > 0; present_p <- pvec > 0
  retention <- if (sum(present_b) == 0) NA_real_ else sum(present_b & present_p) / sum(present_b)
  tibble(
    bray_curtis = bray[base, post],
    bray_resistance = 1 - bray[base, post],
    aitchison = ait[base, post],
    taxon_retention = retention,
    baseline_dominance = max(bvec),
    post_dominance = max(pvec),
    dominance_change_abs = abs(max(pvec) - max(bvec))
  )
}

metrics <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  bind_cols(pairs[i, ], metric_for_pair(pairs$baseline_run[i], pairs$post_run[i]))
})) |>
  left_join(alpha |> rename(baseline_run = run_accession, baseline_richness = richness,
                             baseline_shannon = shannon, baseline_simpson = simpson,
                             baseline_alpha_dominance = dominance), by = "baseline_run") |>
  left_join(alpha |> rename(post_run = run_accession, post_richness = richness,
                             post_shannon = shannon, post_simpson = simpson,
                             post_alpha_dominance = dominance), by = "post_run") |>
  mutate(
    richness_change_abs = abs(post_richness - baseline_richness),
    shannon_change_abs = abs(post_shannon - baseline_shannon),
    simpson_change_abs = abs(post_simpson - baseline_simpson)
  )

write_csv(metrics, file.path(out_dir, "quarter_level_resistance_metrics.csv"))
write_csv(alpha, file.path(out_dir, "alpha_diversity_by_sample.csv"))
cat("[OK] wrote", nrow(metrics), "baseline-to-post quarter trajectories\n")
