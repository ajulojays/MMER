#!/usr/bin/env python3
"""Reconstruct and validate PRJEB38332 cow × quarter × timepoint metadata.

The inference is intentionally explicit and auditable:
- S1-20 = T1; S21-40 = T2; S41-60 = T3.
- Within each 20-sample block, cows occur 270,355,366,321,365.
- Four positions per cow repeat as FL,FR,RR,RL.
- The anatomical key is supported by deposited aliases 270AS/AD/PD/PS.
"""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/metadata/PRJEB38332_ena_report.tsv"
OUT = ROOT / "data/metadata/PRJEB38332_sample_map.csv"

COW_ORDER = ["270", "355", "366", "321", "365"]
QUARTERS = [
    ("FL", "AS", "untreated_control"),
    ("FR", "AD", "teat_sealant"),
    ("RR", "PD", "cephalonium"),
    ("RL", "PS", "cloxacillin"),
]
TIMEPOINTS = {
    1: ("T1", "dry_off_baseline"),
    2: ("T2", "calving"),
    3: ("T3", "5_DIM"),
}

def extract_s(alias: str) -> int:
    m = re.search(r"_S(\d+)$", str(alias))
    if not m:
        raise ValueError(f"Cannot extract S-number from alias: {alias}")
    return int(m.group(1))

def infer(s: int):
    if not 1 <= s <= 60:
        raise ValueError(f"Unexpected S-number: {s}")
    block = (s - 1) // 20 + 1
    within = (s - 1) % 20
    cow_idx = within // 4
    q_idx = within % 4
    timepoint, phase = TIMEPOINTS[block]
    cow = COW_ORDER[cow_idx]
    quarter, anatomical_code, treatment = QUARTERS[q_idx]
    return cow, quarter, anatomical_code, treatment, timepoint, phase

def main():
    df = pd.read_csv(SRC, sep="\t", dtype=str)
    required = {"run_accession", "sample_accession", "fastq_ftp", "sample_alias", "sample_title"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing ENA columns: {sorted(missing)}")

    df["sample_number"] = df["sample_alias"].map(extract_s)
    inferred = df["sample_number"].map(infer)
    df[["cow_id", "quarter", "anatomical_code", "treatment", "timepoint", "phase"]] = pd.DataFrame(
        inferred.tolist(), index=df.index
    )

    # Strong integrity checks.
    assert len(df) == 60, f"Expected 60 runs, found {len(df)}"
    assert df["sample_number"].nunique() == 60
    assert set(df["sample_number"]) == set(range(1, 61))
    assert set(df["cow_id"]) == set(COW_ORDER)
    assert (df.groupby(["cow_id", "quarter"])["timepoint"].nunique() == 3).all()
    assert (df.groupby(["cow_id", "timepoint"])["quarter"].nunique() == 4).all()

    # Alias prefix should agree with inferred cow for all deposited samples.
    alias_cow = df["sample_alias"].str.extract(r"^(\d+)", expand=False)
    if not (alias_cow == df["cow_id"]).all():
        bad = df.loc[alias_cow != df["cow_id"], ["sample_alias", "sample_number", "cow_id"]]
        raise AssertionError(f"Alias/cow inference mismatch:\n{bad}")

    # Direct anatomical-key sanity check where encoded in the first four aliases.
    expected = {1: "270AS", 2: "270AD", 3: "270PD", 4: "270PS"}
    for s, prefix in expected.items():
        alias = df.loc[df.sample_number.eq(s), "sample_alias"].iloc[0]
        assert alias.startswith(prefix), (s, alias, prefix)

    # Split paired FASTQ URLs for downstream reproducibility.
    fastqs = df["fastq_ftp"].str.split(";", expand=True)
    df["fastq_r1"] = "ftp://" + fastqs[0].str.replace(r"^ftp://", "", regex=True)
    df["fastq_r2"] = "ftp://" + fastqs[1].str.replace(r"^ftp://", "", regex=True)

    keep = [
        "sample_number", "run_accession", "sample_accession", "secondary_sample_accession",
        "experiment_accession", "sample_alias", "sample_title", "cow_id", "quarter",
        "anatomical_code", "treatment", "timepoint", "phase", "fastq_r1", "fastq_r2",
    ]
    df = df.sort_values("sample_number")[keep]
    df.to_csv(OUT, index=False)
    print(f"[OK] wrote {len(df)} validated samples -> {OUT}")
    print(df.groupby(["timepoint", "treatment"]).size().unstack(fill_value=0))

if __name__ == "__main__":
    main()
