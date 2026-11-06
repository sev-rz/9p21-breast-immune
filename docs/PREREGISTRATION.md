# Analysis pre-registration (written before the final model fits)

Timestamp: 2026-10-16 18:43 UTC

Disclosure: an exploratory pilot (scripts/pilot.py: BRCA, shallow loss, 13 loci, no multiplicity control) was run before this note. Its results informed the choice of tests, not the thresholds below.

- Primary outcome: logit(LF), PanImmune methylation leukocyte fraction (TCGA primary tumours, sample type 01, one per sample).
- Primary exposure: 9p21.3 loss = GISTIC thresholded call <= -1 at CDKN2A (proxy for the CDKN2A/B + MTAP + IFN-cluster unit; deep-loss and per-gene versions are secondary).
- Primary adjustment set: PAM50 subtype (BRCA only; cancer-type fixed effects pan-cancer), ABSOLUTE purity, leave-arm-out burden (fraction of genes with |call| >= 1 outside the 9p arm).
- Smallest effect size of interest (SESOI): 0.10 logit units. Equivalence declared if the 90% CI lies inside +/-0.10.
- Outcome classes: Locus-specific (adjusted 95% CI excludes 0 and estimate keeps >=50% of the subtype-only estimate); Burden/purity-driven (subtype-only effect CI excludes 0, adjusted CI includes 0, ratio <50%); Equivalent to null (90% CI inside SESOI); Underpowered (CI includes 0 and extends beyond SESOI); Suppressed (sign flips on adjustment).
- Multiplicity: Benjamini-Hochberg within each experiment family; everything outside the primary contrast is exploratory.
- A-priori positive-control cancers (published 9p21/9p loss immune-cold effects): LUAD, LUSC, HNSC, SKCM. All other cancers are reported as a family without selection.
- Feasibility gate for breakpoint classes: >= 30 tumours in an informative class; otherwise reported as non-identifiable.
- Burden uses gene-level GISTIC calls (gene-count fraction), not base-pair-weighted FGA; global FGA/AS/ploidy from ABSOLUTE/seg files are covariates in sensitivity models.
- SKCM: only -06 metastatic samples allowed (deviation: prereg said sample type 01).
