from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

META = (
    ROOT
    / "data/pilot_cow270/metadata/cow270_sample_map.csv"
)

RAW = ROOT / "data/pilot_cow270/raw"


def load():
    assert META.exists()
    return pd.read_csv(META, dtype={"cow_id": str})


def test_twelve_samples():
    assert len(load()) == 12


def test_single_cow():
    assert set(load()["cow_id"]) == {"270"}


def test_four_quarters():
    assert set(load()["quarter"]) == {"FL", "FR", "RR", "RL"}


def test_three_timepoints():
    assert set(load()["timepoint"]) == {"T1", "T2", "T3"}


def test_four_treatments():
    assert set(load()["treatment"]) == {
        "untreated_control",
        "teat_sealant",
        "cephalonium",
        "cloxacillin",
    }


def test_each_quarter_is_longitudinal():
    df = load()
    counts = df.groupby("quarter").size()
    assert (counts == 3).all()


def test_each_timepoint_has_four_quarters():
    df = load()
    counts = df.groupby("timepoint").size()
    assert (counts == 4).all()


def test_all_fastqs_exist():
    df = load()

    for run in df["run_accession"]:
        assert (RAW / f"{run}_R1.fastq.gz").exists()
        assert (RAW / f"{run}_R2.fastq.gz").exists()


def test_expected_fastq_count():
    assert len(list(RAW.glob("*.fastq.gz"))) == 24
