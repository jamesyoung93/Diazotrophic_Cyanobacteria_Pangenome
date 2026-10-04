# Beyond nif

Data and code for **Predicting cyanobacterial diazotrophy and prioritizing candidate protein families beyond Nif**.

The analysis uses 426 cyanobacterial assemblies and 2,286 protein families. It evaluates prediction across held-out genus groups and supplies two candidate inventories: 476 Model-Supported and 981 Highly Pure families. The current experimental shortlist contains **167 families: 22 Tier A, 140 Tier B and 5 Highly Pure**. External annotations guide experimental selection; they do not independently validate biological function.

## Read the results

| Resource | Contents |
|---|---|
| [Candidate workbook](manuscript_release_2026_10/annotations/results/current167/candidate_families.xlsx) | Formatted candidate tables |
| [167-family shortlist](manuscript_release_2026_10/annotations/results/current167/corrected_TopFamilies.csv) | Current scores and tiers |
| [Complete inventory](manuscript_release_2026_10/annotations/results/current167/corrected_AllComposite.csv) | All 1,457 families and annotations |
| [Figures and captions](manuscript_release_2026_10/figures/) | Figures 1–5 and source material |
| [Prediction results](manuscript_release_2026_10/models/reference/) | Fold metrics, predictions, feature masks and sensitivities |
| [Cohort sensitivity](manuscript_release_2026_10/docs/COHORT_AND_TAXONOMY_SENSITIVITY.md) | Assembly status, sample profiles, taxonomy and prevalence |
| [Annotation methods](manuscript_release_2026_10/docs/FUNCTIONAL_ANNOTATION.md) | Product-name patterns, keyword bins and interpretation limits |

## Reproduce the release

Use Python 3.10.16. From the repository root:

```sh
python -m pip install -r requirements.txt
python manuscript_release_2026_10/scripts/run_release_checks.py --out release_verification
```

Choose a new output directory. This command checks file integrity, reconstructs annotations, reproduces current scores and figure tables, and independently checks saved model metrics. It performs no model refits. GitHub Actions runs these checks on Windows.

For model refits, use the [model instructions](manuscript_release_2026_10/models/README.md). The [release guide](manuscript_release_2026_10/README.md) describes the data and reproducibility limits. The separate [sequence-processing workflow](unified_pipeline_clean/README.md) supplies acquisition, HMM scanning and clustering code.

The working tree contains the current release and its required source material. Prior versions are available through Git history.
