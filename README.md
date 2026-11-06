# Immune-cold breast cancer and loss of 9p21.3

Course project (CBIO 7102, Fall 2026), Sevinch Rakhmon-Zoda and Jenna Khanfar.

**Question.** Does losing 9p21.3 make breast tumours immune-cold, or is it just tumour purity and DNA damage?

**Answer.** In 915 TCGA breast tumours, purity explains most of the variation in immune share. After adjusting for subtype, DNA damage and purity, 9p21.3 loss has no detectable effect: 0.06 logit units, 95% range -0.02 to 0.13. Full write-up: `docs/progress_report.pdf`.

## Folders

```
data/raw/        downloaded files, never edited (not in git)
data/processed/  built by the scripts (not in git)
metadata/        sources.tsv (where data come from), checksums.sha256
scripts/         analysis code
results/         tables written by the scripts (not in git)
docs/            progress report, pre-registration
```

Commit what a person wrote, not what a program can download or rebuild.

## Run it

Python 3 with pandas, numpy, scipy, statsmodels, scikit-learn.

1. Download the files in `metadata/sources.tsv` into `data/raw/` (GDC link pattern: `https://api.gdc.cancer.gov/data/<file id>`; subtype with `scripts/cbp.py` and `cbp2.py`), then check them: `shasum -a 256 -c metadata/checksums.sha256`
2. `python scripts/pack_gistic.py` (once; packs the big copy-number file)
3. `python scripts/x00_assemble.py` (builds the analysis table: 915 breast tumours, one row each)
4. Run the experiments; each writes CSV files to `results/`:

| Script | What it does |
|---|---|
| `x01_brca_core.py` | breast: adjustment ladder, rival models, dose and gain, doubly robust estimate |
| `x02*.py` | other cancers, positive controls, deletion classes, gene vs arm loss |
| `x03*.py` | genome scan with 1,000 shuffles |
| `x04_external.py` | tumour RNA, METABRIC, cell lines (needs internet) |
| `x05_robust.py` | the 15 variations of the main estimate |

## Method

- Outcome: logit of leukocyte fraction (from DNA methylation). Exposure: CDKN2A copy-number call of -1 or lower.
- Model: `logit(LF) = b0 + b1*loss + subtype + burden + purity`, least squares with robust (HC3) errors.
- Rules fixed in advance (`docs/PREREGISTRATION.md`): smallest effect of interest 0.10; positive-control cancers LUAD, LUSC, HNSC, SKCM.
- Four rivals: R1 interferon genes, R2 CDKN2A/MTAP, R3 whole 9p arm, R4 burden and purity only.

Data: TCGA PanCanAtlas (GDC), cBioPortal, METABRIC, CCLE; all public, accessed 16 October 2026. Code: MIT license.
