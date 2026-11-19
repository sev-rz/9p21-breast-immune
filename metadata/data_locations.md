# Where the data live

Raw data are public. They are never committed to Git. They exist in two places, and the
laptop copy is rebuilt by `pixi run download`.

| Place | Path | Used by |
|---|---|---|
| Laptop (this repo) | `data/raw/` (downloads), `data/processed/` (built by scripts) | `scripts/`, Snakemake |
| MBZUAI cluster (sc-hub project `brca-loci-vs-cna`) | `/l/users/sevinch.zoda/schub/projects/brca-loci-vs-cna/work/data/` | the cluster analysis, which this repo's pipeline mirrors |

Cluster paths follow sc-hub's `projects/<project>/work/` layout (checked through the sc-hub
file listing, 2026-11-19); the data stay on the cluster and are read where they are.

## Raw files shared by both (same source, different file names)

| Repo file (`data/raw/`) | Cluster file (`work/data/tcga/`) | Source |
|---|---|---|
| `leuk.tsv` | `leuk_estimate.tsv` | GDC PanCanAtlas, Thorsson et al. 2018 |
| `mastercalls.txt` | `absolute_calls.txt` | GDC PanCanAtlas, ABSOLUTE |
| `absolute_scores.tsv` | `absolute_scores.tsv` | GDC PanCanAtlas |
| `seg_scores.tsv` | `seg_based_scores.tsv` | GDC PanCanAtlas |
| `gistic_thr.tsv` (589 MB) | `gistic_thresholded.tsv` | GDC PanCanAtlas, GISTIC2 thresholded calls |
| `cbp_patient.csv`, `cbp_clinical.csv` | not copied | cBioPortal API (`brca_tcga_pan_can_atlas_2018`) |

Download links and checksums: `sources.tsv`, `checksums.sha256`.

## Cluster-only raw files (not used by this repo's scripts)

Under `work/data/tcga/`: `cibersort.tsv`, `clinical_pancan.tsv`, `mitcr_stats.tsv`,
`neoantigens.tsv`, `scores_160.tsv.gz`, `tcga_subtype.tsv`, `til_*` (tumour-infiltrating
lymphocyte maps and Table S1). Under `work/data/raw/`: `data_cna.txt` (METABRIC copy number),
`thorsson_supp2-4.xlsx` (small files, likely not the real supplements: check before use).
Also: `data/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz` (Wu et al. 2021, GEO GSE176078),
`work/ref/` (UCSC cytoBand and refGene, GENCODE v44 annotation, `gene_coords_hg38.csv`), and CPTAC and
METABRIC extracts in `work/` (`cptac_*`, `mb_*`).
Exact download links for these are in the sc-hub journal cells, not recorded here yet.

## Processed on the cluster (`work/`)

`gistic_brca.pkl`, `gistic_5cancers.pkl` (packed GISTIC calls), `e0_brca_model_table_n915.pkl`
(analysis table, n = 915), `e*_*.csv` (results). On the laptop the equivalents are
`data/processed/thr.npy`, `panc.pkl` and `results/`.

## Not found

The C1-C6 immune-subtype table and per-tumour HLA LOH are not on the cluster or here
(see the sc-hub handover of 2026-10-18).
