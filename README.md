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
docs/                  pre-registration, progress report, final figures, decisions
workflow/, config/     Snakemake pipeline and its settings
containers/            Dockerfile and Apptainer recipe
pixi.toml, pixi.lock   software environment (exact versions)
```

Rule of thumb: commit what a person wrote; do not commit what a program can
download or rebuild. That is why `data/` and `results/` are empty in git.

## How to reproduce

You need [Pixi](https://pixi.prefix.dev/latest/) and an internet connection.

```bash
pixi install          # builds the exact environment recorded in pixi.lock
pixi run reproduce    # downloads the raw data and rebuilds every result (Snakemake)
```

A GitHub Action (`.github/workflows/reproduce.yml`) checks this on a fresh clone.
Or run it in Docker: see `containers/README.md`. `pixi run dag` draws the workflow.

The Snakemake rules in `workflow/Snakefile` run these scripts in order
(`metadata/README.md` explains the raw files and checksums):

| Script | What it does |
|---|---|
| `download_data.py`, `cbp.py`, `cbp2.py` | download raw data (GDC, cBioPortal) and verify checksums |
| `pack_gistic.py` | pack the big copy-number file (once) |
| `x00_assemble.py` | build the analysis table |
| `x01_brca_core.py` | breast: variance split, adjustment ladder, rivals, dose and gain, doubly robust estimate |
| `x02_pancancer.py`, `x02b_classes_pos.py`, `x02c_arm_pos.py` | other cancers, positive controls, deletion classes, gene vs arm loss |
| `x03_scan.py`, `x03b_symmetry.py` | genome scan with 1,000 shuffles, loss vs gain |
| `x04_external.py` | tumour RNA, METABRIC and cell lines (cBioPortal API) |
| `x05_robust.py` | the 15 variations of the main estimate |
| `make_figs.py`, `make_tables.py`, `make_numbers.py` | figures, tables and numbers for the report (`report/`) |

`audit.py`, `pilot.py` are early exploration, kept for the record. Random seeds are fixed in the scripts.

### Hardware and run time

Any laptop or Linux machine with about 8 GB of memory and 3 GB of free disk (the raw
GISTIC file is 589 MB and is packed to a smaller array). Needs an internet connection for
the downloads. The genome scan (1,000 label shuffles) is the slow step; to test quickly,
lower `n_permutations` in `config/config.yaml` and restore 1000 for the real results.

### Expected outputs

After `pixi run reproduce`:

- `results/*.csv`: 26 result tables (`x1a_shapley.csv`, `x2_ladder.csv`, `x5_scan.csv`, `x14_robust.csv`, ...)
- `report/figures/`, `report/tables/`, `report/numbers.tex`: figures, tables and numbers for the LaTeX report
- `reports/9p21-report.html`: the Quarto report, computed from the result tables
- headline numbers to compare: R² 0.67 of which purity 0.55; adjusted 9p21.3 effect 0.06 logit
  units (95% interval -0.02 to 0.13); 9p21.3 ranks 104th of 756 regions

## Method in four lines

- Outcome: logit of leukocyte fraction (LF, from DNA methylation). Exposure: CDKN2A
  copy-number call of -1 or lower (320 of 915 tumours).
- Model: `logit(LF) = b0 + b1*loss + subtype + burden + purity`, least squares with robust (HC3) errors.
- Rules fixed before the final fits (`docs/PREREGISTRATION.md`): smallest effect of
  interest 0.10 logit units; four positive-control cancers (LUAD, LUSC, HNSC, SKCM).
- Four rival explanations: R1 interferon genes, R2 CDKN2A/MTAP, R3 whole 9p arm, R4 burden and purity only.

## Data availability

All data are public and none are stored in this repository: TCGA PanCanAtlas files from the GDC
(downloaded by `scripts/download_data.py`), PAM50 subtype from cBioPortal
(`scripts/download_cbioportal.py`), METABRIC and CCLE through the cBioPortal API. Sources, access
dates and checksums are in `metadata/`; where the files sit on the laptop and on the cluster is in
`metadata/data_locations.md`.

## Data and credits

All data are public (TCGA PanCanAtlas via the GDC, cBioPortal, METABRIC, CCLE);
sources and access dates are in `metadata/`. Thorsson et al. 2018 (Immunity)
provides the leukocyte fraction. An AI assistant helped with data processing,
analysis and drafting; the authors are responsible for the interpretation.
