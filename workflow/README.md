# workflow

`Snakefile` runs the analysis steps in the right order, one rule per script
(download → pack GISTIC → assemble → x01…x05 → figures, tables, numbers).
Settings are in `config/config.yaml`.

```bash
pixi run run-snakemake        # everything (downloads ~600 MB first time)
pixi run dag                  # draw the graph to results/dag.pdf
pixi run snakemake -n         # dry run: what would be rebuilt
```

The genome scan (`x03_scan`) is the slow step. For a quick test lower
`n_permutations` in `config/config.yaml`; restore 1000 for the real results.
