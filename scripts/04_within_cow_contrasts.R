#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})
args <- commandArgs(trailingOnly = FALSE)
script <- sub("--file=", "", args[grep("--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
x <- read_csv(file.path(root, "results/resistance/quarter_level_resistance_metrics.csv"), show_col_types = FALSE)

# Control-adjusted excess displacement within the same cow and timepoint.
control <- x |> filter(treatment == "untreated_control") |>
  select(cow_id, timepoint,
         control_bray = bray_curtis,
         control_aitchison = aitchison,
         control_retention = taxon_retention)

out <- x |> left_join(control, by = c("cow_id", "timepoint")) |>
  mutate(
    excess_bray_vs_control = bray_curtis - control_bray,
    excess_aitchison_vs_control = aitchison - control_aitchison,
    retention_difference_vs_control = taxon_retention - control_retention
  )
write_csv(out, file.path(root, "results/resistance/control_adjusted_resistance.csv"))
cat("[OK] within-cow control-adjusted contrasts written\n")
