#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data/metadata/PRJEB38332_sample_map.csv"


def load_metadata():
    assert METADATA.exists(), "Run scripts/00_build_metadata.py first"
    return pd.read_csv(METADATA, dtype={"cow_id": str})


def test_sample_count_is_60():
    df = load_metadata()
    assert len(df) == 60


def test_run_accessions_are_unique():
    df = load_metadata()
    assert df["run_accession"].nunique() == 60


def test_sample_numbers_are_complete():
    df = load_metadata()
    assert sorted(df["sample_number"].tolist()) == list(range(1, 61))


def test_five_expected_cows_present():
    df = load_metadata()
    expected = {"270", "355", "366", "321", "365"}
    assert set(df["cow_id"].unique()) == expected


def test_four_expected_quarters_present():
    df = load_metadata()
    assert set(df["quarter"].unique()) == {"FL", "FR", "RR", "RL"}


def test_expected_treatments_present():
    df = load_metadata()
    assert set(df["treatment"].unique()) == {
        "untreated_control",
        "teat_sealant",
        "cephalonium",
        "cloxacillin",
    }


def test_three_timepoints_present():
    df = load_metadata()
    assert set(df["timepoint"].unique()) == {"T1", "T2", "T3"}


def test_each_cow_quarter_has_three_timepoints():
    df = load_metadata()
    counts = df.groupby(["cow_id", "quarter"]).size()
    assert len(counts) == 20
    assert (counts == 3).all()


def test_each_cow_timepoint_has_four_quarters():
    df = load_metadata()
    counts = df.groupby(["cow_id", "timepoint"]).size()
    assert len(counts) == 15
    assert (counts == 4).all()


def test_each_treatment_timepoint_has_five_samples():
    df = load_metadata()
    counts = df.groupby(["treatment", "timepoint"]).size()
    assert len(counts) == 12
    assert (counts == 5).all()


def test_each_treatment_has_fifteen_samples():
    df = load_metadata()
    counts = df.groupby("treatment").size()
    assert (counts == 15).all()


def test_quarter_treatment_assignment_is_fixed():
    df = load_metadata()

    expected = {
        "FL": "untreated_control",
        "FR": "teat_sealant",
        "RR": "cephalonium",
        "RL": "cloxacillin",
    }

    for quarter, treatment in expected.items():
        observed = set(df.loc[df["quarter"] == quarter, "treatment"])
        assert observed == {treatment}


def test_longitudinal_trajectories_are_complete():
    df = load_metadata()

    expected_timepoints = {"T1", "T2", "T3"}

    for (_, _), group in df.groupby(["cow_id", "quarter"]):
        assert set(group["timepoint"]) == expected_timepoints
