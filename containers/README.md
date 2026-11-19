# containers

Recipes for container images of this project's environment:

- `Dockerfile`: builds an image whose software comes from `pixi.lock`, so the
  container and your laptop run exactly the same packages.
- `apptainer.def`: the same for HPC clusters that use Apptainer.

```bash
docker build --platform linux/amd64 -f containers/Dockerfile -t 9p21-breast-immune:0.1.0 .
docker run --rm -v "$PWD":/work -w /work 9p21-breast-immune:0.1.0 snakemake --cores 4
```

Data (`data/`) and results are never baked into the image; the project folder
is mounted at `/work`. Saved images (`*.tar.gz`, `*.sif`) are not committed.
