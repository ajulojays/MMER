#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENA_REPORT="${ROOT_DIR}/data/metadata/PRJEB38332_ena_report.tsv"
OUTDIR="${ROOT_DIR}/data/raw"
N_PARALLEL=16

mkdir -p "${OUTDIR}"

if [[ ! -f "${ENA_REPORT}" ]]; then
    echo "[ERROR] ENA report not found:"
    echo "        ${ENA_REPORT}"
    exit 1
fi

echo "============================================================"
echo " MMER parallel FASTQ downloader"
echo "============================================================"
echo "[INFO] ENA report : ${ENA_REPORT}"
echo "[INFO] Output dir : ${OUTDIR}"
echo "[INFO] Workers    : ${N_PARALLEL}"
echo

DOWNLOAD_LIST="$(mktemp)"
trap 'rm -f "${DOWNLOAD_LIST}"' EXIT

python - "${ENA_REPORT}" "${OUTDIR}" <<'PY' > "${DOWNLOAD_LIST}"
import sys
from pathlib import Path
import pandas as pd

ena_report = Path(sys.argv[1])
outdir = Path(sys.argv[2])

df = pd.read_csv(ena_report, sep="\t", dtype=str)

required = {"run_accession", "fastq_ftp"}

missing = required - set(df.columns)
if missing:
    raise SystemExit(
        f"[ERROR] ENA report missing required columns: {sorted(missing)}"
    )

records = []

for _, row in df.iterrows():
    run = str(row["run_accession"]).strip()
    ftp = str(row["fastq_ftp"]).strip()

    if not run or run == "nan":
        continue
    if not ftp or ftp == "nan":
        continue

    urls = [u.strip() for u in ftp.split(";") if u.strip()]

    if len(urls) != 2:
        raise SystemExit(
            f"[ERROR] Expected paired FASTQs for {run}, found {len(urls)}: {urls}"
        )

    for read_number, url in enumerate(urls, start=1):

        # HTTPS is generally more robust than anonymous FTP.
        if url.startswith("ftp://"):
            url = "https://" + url[len("ftp://"):]
        elif not url.startswith(("http://", "https://")):
            url = "https://" + url

        outfile = outdir / f"{run}_R{read_number}.fastq.gz"

        records.append((run, read_number, url, str(outfile)))

# Ensure no accidental duplicates.
seen = set()
unique = []

for rec in records:
    key = (rec[0], rec[1])
    if key not in seen:
        seen.add(key)
        unique.append(rec)

if len(unique) != 120:
    raise SystemExit(
        f"[ERROR] Expected 120 FASTQ entries (60 samples × 2), found {len(unique)}"
    )

for run, read_number, url, outfile in unique:
    print(f"{run}\t{read_number}\t{url}\t{outfile}")
PY

N_EXPECTED="$(wc -l < "${DOWNLOAD_LIST}")"

echo "[INFO] FASTQ entries identified: ${N_EXPECTED}"

if [[ "${N_EXPECTED}" -ne 120 ]]; then
    echo "[ERROR] Expected 120 FASTQ entries."
    exit 1
fi

download_one() {
    line="$1"

    IFS=$'\t' read -r run read_number url outfile <<< "${line}"

    echo "[GET] ${run} R${read_number}"

    # wget --continue resumes interrupted downloads.
    wget \
        --continue \
        --tries=8 \
        --timeout=60 \
        --read-timeout=60 \
        --retry-connrefused \
        --waitretry=5 \
        --no-verbose \
        --output-document="${outfile}" \
        "${url}"

    if [[ ! -s "${outfile}" ]]; then
        echo "[ERROR] Empty output: ${outfile}" >&2
        return 1
    fi
}

export -f download_one

echo
echo "[INFO] Starting ${N_PARALLEL} parallel downloads..."
echo

xargs \
    -P "${N_PARALLEL}" \
    -d '\n' \
    -I {} \
    bash -c 'download_one "$1"' _ "{}" \
    < "${DOWNLOAD_LIST}"

echo
echo "============================================================"
echo " Download complete — validating dataset"
echo "============================================================"

N_DOWNLOADED="$(
    find "${OUTDIR}" \
        -maxdepth 1 \
        -type f \
        -name 'ERR*_R[12].fastq.gz' \
        | wc -l
)"

echo "[INFO] FASTQ files present: ${N_DOWNLOADED}/120"

if [[ "${N_DOWNLOADED}" -ne 120 ]]; then
    echo "[ERROR] Expected 120 FASTQ files, found ${N_DOWNLOADED}"
    exit 1
fi

echo "[INFO] Checking for zero-byte files..."

EMPTY="$(
    find "${OUTDIR}" \
        -maxdepth 1 \
        -type f \
        -name 'ERR*_R[12].fastq.gz' \
        -size 0 \
        | wc -l
)"

if [[ "${EMPTY}" -ne 0 ]]; then
    echo "[ERROR] Found ${EMPTY} zero-byte FASTQ files."
    exit 1
fi

echo "[INFO] Running gzip integrity checks with ${N_PARALLEL} workers..."

find "${OUTDIR}" \
    -maxdepth 1 \
    -type f \
    -name 'ERR*_R[12].fastq.gz' \
    -print0 \
    | xargs -0 -P "${N_PARALLEL}" -n 1 gzip -t

echo
echo "[PASS] 120/120 FASTQ files present."
echo "[PASS] All gzip integrity checks passed."
echo
du -sh "${OUTDIR}"

echo
echo "============================================================"
echo " MMER raw sequencing dataset READY"
echo "============================================================"
