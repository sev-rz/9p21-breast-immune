# Immune-cold breast cancer and loss of 9p21.3

Course project (CBIO 7102, Fall 2026) by Sevinch Rakhmon-Zoda and Jenna Khanfar.

**Question.** Does losing the chromosome region 9p21.3 (CDKN2A/B, MTAP and the
type-I interferon genes) make breast tumours "immune cold", or is that just
tumour purity and general DNA damage?

**Short answer.** In 915 TCGA breast tumours, purity explains most of the variation
in immune share (0.55 of R² 0.67). After adjusting for subtype, DNA damage and
purity, 9p21.3 loss has no detectable effect (0.06 logit units, 95% range
-0.02 to 0.13; underpowered, not proven zero). No rival explanation improves the
model, 9p21.3 ranks 104th of 756 regions, and the same pipeline finds drops in
head and neck cancer and melanoma. Details: `docs/progress_report.pdf`.

## Folder layout

```
README.md              this file
LICENSE, CITATION.cff  MIT license (code), how to cite
metadata/              where the data come from (sources.tsv) and checksums
data/raw/              downloaded files, never edited (not in git)
data/processed/        files the scripts build from raw (not in git)
scripts/               analysis code, run in the order below
results/               tables written by the scripts (not in git)
docs/                  pre-registration, progress report, final figures
```

Rule of thumb: commit what a person wrote; do not commit what a program can
download or rebuild. That is why `data/` and `results/` are empty in git.

## How to reproduce

Python 3 with pandas, numpy, scipy, statsmodels, scikit-learn and matplotlib.
Run everything from the project folder.

1. **Get the raw data.** Follow `metadata/README.md`: download the files listed in
   `metadata/sources.tsv` into `data/raw/`, then check them with
   `shasum -a 256 -c metadata/checksums.sha256`.
2. **Pack the big copy-number file** (once): `python scripts/pack_gistic.py`
3. **Build the analysis table:** `python scripts/x00_assemble.py`
4. **Run the experiments** (each writes CSV files to `results/`):

| Script | What it does |
|---|---|
| `x01_brca_core.py` | breast: variance split, adjustment ladder, rivals, dose and gain, doubly robust estimate |
| `x02_pancancer.py`, `x02b_classes_pos.py`, `x02c_arm_pos.py` | other cancers, positive controls, deletion classes, gene vs arm loss |
| `x03_scan.py`, `x03b_symmetry.py` | genome scan with 1,000 shuffles, loss vs gain |
| `x04_external.py` | tumour RNA, METABRIC and cell lines (needs internet, cBioPortal API) |
| `x05_robust.py` | the 15 variations of the main estimate |
| `make_figs.py`, `make_tables.py`, `make_numbers.py` | figures and tables for the report |

`audit.py`, `pilot.py`, `cbp.py`, `cbp2.py` are the early exploration and
download helpers, kept for the record. Random seeds are fixed in the scripts.

## Method in four lines

- Outcome: logit of leukocyte fraction (LF, from DNA methylation). Exposure: CDKN2A
  copy-number call of -1 or lower (320 of 915 tumours).
- Model: `logit(LF) = b0 + b1*loss + subtype + burden + purity`, least squares with robust (HC3) errors.
- Rules fixed before the final fits (`docs/PREREGISTRATION.md`): smallest effect of
  interest 0.10 logit units; four positive-control cancers (LUAD, LUSC, HNSC, SKCM).
- Four rival explanations: R1 interferon genes, R2 CDKN2A/MTAP, R3 whole 9p arm, R4 burden and purity only.

## Data and credits

All data are public (TCGA PanCanAtlas via the GDC, cBioPortal, METABRIC, CCLE);
sources and access dates are in `metadata/`. Thorsson et al. 2018 (Immunity)
provides the leukocyte fraction. An AI assistant helped with data processing,
analysis and drafting; the authors are responsible for the interpretation.
