#!/usr/bin/env bash
set -euo pipefail

# Canonical Van Beeck external-validation runner.
# Executes the frozen analysis and then generates publication-ready figures.

ROOT="${MMER_ROOT:-$HOME/MeteG/MMER_repo}"
cd "$ROOT"

mkdir -p logs/external_generalization

Rscript scripts/26j_vanbeeck_external_validation_only.R \
  2>&1 | tee logs/external_generalization/vanbeeck_external_validation_only.log

Rscript scripts/26k_vanbeeck_external_validation_plots.R \
  2>&1 | tee logs/external_generalization/vanbeeck_external_validation_plots.log

echo "[PASS] Van Beeck external validation + figures complete"
echo "Results: $ROOT/results/vanbeeck/external_validation_only"
echo "Figures: $ROOT/figures/vanbeeck_external_validation"
