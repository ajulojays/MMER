#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
META="$ROOT/data/metadata/PRJEB38332_sample_map.csv"
OUT="$ROOT/data/raw"
mkdir -p "$OUT"

if [[ ! -f "$META" ]]; then
  python "$ROOT/scripts/00_build_metadata.py"
fi

python - "$META" "$OUT" <<'PY'
import csv, pathlib, subprocess, sys
meta, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
with meta.open() as fh:
    for row in csv.DictReader(fh):
        run = row["run_accession"]
        for read in ("r1", "r2"):
            url = row[f"fastq_{read}"]
            dest = out / f"{run}_{read.upper()}.fastq.gz"
            if dest.exists() and dest.stat().st_size > 0:
                print(f"[SKIP] {dest.name}")
                continue
            print(f"[GET] {run} {read.upper()}")
            subprocess.run(["wget", "-c", "-O", str(dest), url], check=True)
PY
