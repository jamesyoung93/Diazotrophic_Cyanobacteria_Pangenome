# Manuscript data and analysis

This directory contains the data, code, figures and methods notes for *Beyond nif*. The primary panel contains 426 assembly versions (112 nifHDK-positive and 314 comparison genomes) and 2,286 retained protein families. The inventories contain 476 residual-association and 981 high-purity families. The experimental shortlist contains 22 Tier A, 140 Tier B and 5 high-purity families.

## Contents

- [Annotation results](annotations/README.md): current scores, candidate workbook, reconstructed study annotations, mapping sensitivities and Figures 3/5.
- [Prediction analysis](models/README.md): fixed-matrix prediction, genome-length baselines, saved predictions and model refits.
- [Figures](figures/): the five manuscript figures and captions, including the Figure 4 comparison source tables.
- [Supporting tables](tables/README.md): the related-atlas comparison, nif-keyword audit and proteomics source index used in the manuscript.
- [Cohort and taxonomy sensitivity](docs/COHORT_AND_TAXONOMY_SENSITIVITY.md): dated assembly-status, retained-profile, genus and prevalence analyses.
- [Functional annotation](docs/FUNCTIONAL_ANNOTATION.md): category definitions and limits of product-name matching and co-occurrence tests.

## Reproduce

Use Python 3.10.16 and the pinned dependencies. From the repository root:

```sh
python -m pip install -r requirements.txt
python manuscript_release_2026_10/scripts/run_release_checks.py --out release_verification
```

Choose a new output directory outside this directory. The runner verifies hashes before and after execution, reconstructs annotations against frozen validation targets, and checks saved model metrics and predictions in a copy of the reference outputs. Inputs include dated NCBI responses; replay needs no live input downloads.

This check does not fit models or rebuild sequence clusters from raw genomes. Separate model-refit commands cover 460 fits plus a 20-fit training-prevalence analysis. A clean-environment primary refit preserved masks, class decisions, reported AUC/AP and importance ranks but failed the unchanged strict probability tolerance for 140 logistic-regression predictions (maximum difference 4.67e-10). The [model validation record](models/validation/VALIDATION.md) documents that limit.

## Interpretation and provenance

Clustering and the full-panel vocabulary remain fixed in the prediction sensitivities. Current NCBI status is not backdated to collection, and a verified original retrieval timestamp is unavailable. Three ubiquitous residual-association families remain in the inventory but are excluded from prioritization. Tier labels describe heuristic score bands.

`RELEASE_SHA256.json` covers this directory. Input and reference manifests preserve the exact data used by the analyses. Baseline inputs and comparison outputs remain where a current manuscript result or verification depends on them. Acquisition code revisions and data hashes are identified in the cohort note and model input manifest.

[Terminology and data fields](docs/TERMINOLOGY.md) define the two inventories, model rankings and prioritization tiers. [Related-atlas comparison](docs/RELATED_ATLAS.md) provides the versioned external source, reported performance and limits of the comparison. [Analysis provenance](docs/ANALYSIS_PROVENANCE.md) records filtering and scoring revisions.
