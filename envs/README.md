# envs

Conda environment files for people who do not use Pixi, exported from `pixi.toml` and `pixi.lock`
(these remain the source of truth).

- `environment.yml`: version ranges. `conda env create -p ./conda-env -f envs/environment.yml`
  (regenerate: `pixi workspace export conda-environment envs/environment.yml`).
- `explicit/`: exact package lists per platform, from `pixi.lock`
  (`pixi workspace export conda-explicit-spec envs/explicit`). Build with
  `conda create -p ./env --file envs/explicit/default_<platform>_conda_spec.txt`, your own platform only.
