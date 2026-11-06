# metadata

Files that describe the raw data. They are small, written by hand and committed.

## sources.tsv

Where each raw file comes from. Columns: `name`; `source` (provider); `location`
(download link or study page); `local_file` (where it goes in the project);
`md5` (checksum published by the provider, `NA` when the provider gives none);
`note` (original file name and what we use it for).

Data were accessed on 16 October 2026. The GDC files are open access and are
downloaded from `https://api.gdc.cancer.gov/data/<file id>`. The cBioPortal
studies are read through its public API (`scripts/cbp.py`, `scripts/cbp2.py`,
`scripts/x04_external.py`).

## checksums.sha256

SHA-256 checksum of each raw file we keep. Check with:

```
shasum -a 256 -c metadata/checksums.sha256
```

(`sha256sum -c` on Linux.) The big GISTIC file is not listed here: its md5 is in
`sources.tsv`.

## Raw files

| File | What it holds |
|---|---|
| leuk.tsv | cancer type, TCGA barcode, leukocyte fraction (no header row) |
| mastercalls.txt | ABSOLUTE purity, ploidy and genome doublings per sample |
| absolute_scores.tsv | aneuploidy score and LOH fraction per sample |
| seg_scores.tsv | number of copy-number segments and fraction genome altered |
| cbp_patient.csv | PAM50 subtype per patient (cBioPortal) |
| cbp_clinical.csv | extra per-sample clinical values (cBioPortal) |

Tables are matched on the first 15 characters of the TCGA barcode (the tumour
sample); subtype is matched on the first 12 (the patient).
