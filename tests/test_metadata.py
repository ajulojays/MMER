#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
p = ROOT / "data/metadata/PRJEB38332_sample_map.csv"
assert p.exists(), "Run scripts/00_build_metadata.py first"
df = pd.read_csv(p, dtype={"cow_id": str})
assert len(df) == 60
assert df.run_accession.nunique() == 60
assert df.sample_number.tolist() == list(range(1,61))
assert sorted(df.cow_id.unique()) == sorted(["270","355","366","321","365"])
assert set(df.quarter) == {"FL","FR","RR","RL"}
assert set(df.treatment) == {"untreated_control","teat_sealant","cephalonium","cloxacillin"}
assert set(df.timepoint) == {"T1","T2","T3"}
assert (df.groupby(["cow_id","quarter"]).size() == 3).all()
assert (df.groupby(["cow_id","timepoint"]).size() == 4).all()
print("[PASS] metadata integrity checks")
