# Decision log

- 2026-11-19: project reorganised like the LTEE project: Pixi environment (`pixi.toml`,
  `pixi.lock`), Snakemake pipeline (`workflow/Snakefile`), Docker/Apptainer recipes
  (`containers/`), GitHub Action that runs `pixi run reproduce` from a fresh clone.
- 2026-11-19: `scripts/x03_scan.py` used the name `P` for both the processed-data path
  (from `common.py`) and the number of shuffles, which fails on a clean run. The shuffle
  count is now `NPERM` (default 1000, unchanged); results are not affected.
