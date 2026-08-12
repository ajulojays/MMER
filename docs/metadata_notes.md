# PRJEB38332 metadata reconstruction notes

## Deposited sample pattern
The ENA report contains 60 runs. Sample aliases contain a numeric cow prefix and Illumina sample number (`_S1` ... `_S60`).

The first cow is unusually informative:

- `270AS_S1`
- `270AD_S2`
- `270PD_S3`
- `270PS_S4`

These encode the four mammary-quarter positions and establish the within-cow quarter order. Subsequent aliases preserve the cow prefix and sequential S number even when the anatomical suffix is absent.

## Inferred repeated design
Each 20-sample block contains the same five cows in the same four-sample order:

`270, 355, 366, 321, 365`

Blocks:
- S1–S20 → T1
- S21–S40 → T2
- S41–S60 → T3

Within each cow:
1. FL/AS → untreated control
2. FR/AD → teat sealant
3. RR/PD → cephalonium
4. RL/PS → cloxacillin

## Audit rule
Do not silently change this mapping. Any revision must:
1. cite an external primary-source key or author-provided metadata;
2. update `scripts/00_build_metadata.py`;
3. pass `tests/test_metadata.py`;
4. document the change here.
