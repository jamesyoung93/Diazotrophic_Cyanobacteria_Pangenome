# October 2026 manuscript release

This release contains the frozen data, analysis code, figures and detailed notes for *Beyond nif*. Protein-family content predicts the nifHDK label across held-out genus groups. The two historical candidate inventories serve different purposes. External annotation layers help choose experiments; they do not independently validate genetic requirement or convergent biological support.

## Start here

- [Cohort and taxonomy sensitivity](docs/COHORT_AND_TAXONOMY_SENSITIVITY.md): the 426-genome panel, dated 401-current analysis, retained-profile identity concern, taxonomy and prevalence checks.
- [Functional annotation](docs/FUNCTIONAL_ANNOTATION.md): the meaning and limitations of the legacy categories, product-name patterns and category co-occurrence test.
- [Annotation replay](annotations/README.md): reconstruct source annotations, reproduce the corrected 167-family shortlist and regenerate Figures 3 and 5.
- [Model audit](models/README.md): fixed-matrix prediction, genome-length baselines, saved predictions and optional model refits.
- [Figures](figures/): manuscript figure assets. Figure 4 retains its historical analysis scope; see its source files.
- [Release notes](RELEASE_NOTES.md): changes from the June release and limits of verification.

The primary panel contains 426 assembly versions (112 nifHDK-positive and 314 comparison genomes) and 2,286 retained protein families. The frozen inventories contain 476 historically Model-Supported and 981 Highly Pure families. The current shortlist contains 22 Tier A, 140 Tier B and 5 Highly Pure families. Three ubiquitous Model-Supported families remain in the inventory but are excluded from prioritization. Tier names are heuristic score bands, not levels of validation.

## Reproduce

Use Python 3.10.16 and the pinned dependencies. From the repository root:

```sh
python -m pip install -r manuscript_release_2026_10/requirements.txt
python manuscript_release_2026_10/scripts/run_release_checks.py --out release_verification
```

Choose a new output directory outside the release. The runner verifies release hashes before and after execution, reconstructs annotations against frozen validation targets, and independently checks saved model metrics and predictions in a temporary copy of the reference outputs. It does not fit models. See the model README for the separate 460-fit replay and 20-fit training-prevalence check. Inputs include the dated NCBI responses, so replay does not require live data downloads.

GitHub Actions runs release-integrity checks, the annotation replay and saved-model verification on Windows. These checks do not certify a full end-to-end rebuild from downloaded genomes. In a separate clean environment, primary model fitting preserved masks, class decisions, reported AUC/AP and importance ranks but failed the unchanged strict probability tolerance for 140 logistic-regression predictions (maximum difference 4.67e-10). That failure remains documented in the model bundle; tolerances have not been relaxed.

## Scope and provenance

This directory supersedes `manuscript_release_2026_06` for manuscript interpretation. The June files remain historical records. The original clustering and full-panel vocabulary remain fixed in all model sensitivities. Current NCBI status is not backdated to collection, and a verified original retrieval timestamp is unavailable.

`RELEASE_SHA256.json` covers the released files. The annotation and model bundles also retain their input and reference manifests. Machine-readable legacy fields and filenames remain where needed for reproducibility, including older tier descriptions and proposed-revision status strings. Their historical wording does not override the interpretation above.
