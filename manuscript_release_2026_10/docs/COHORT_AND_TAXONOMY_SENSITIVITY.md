# Cohort and taxonomy sensitivity analysis

The primary analysis uses a fixed historical panel of 426 assembly versions. Dated checks of assembly status, sample identity, taxonomy and training-set prevalence qualify its performance without replacing its candidate inventory. Family-content prediction persists across these checks, but the size and direction of changes depend on the classifier and the metric.

This note accompanies the complete [portable model audit](../models/README.md). The `models` directory contains the original inputs, dated public snapshots, source, saved reference predictions and verification reports.

## Original cohort and archive provenance

The original collection screened 489 cyanobacterial records. The archived RefSeq-accession rule retained GCF records and excluded 63; the completeness criterion excluded none further. The resulting 426 versions contain 112 nifHDK-positive and 314 comparison genomes across 75 stored genus-name tokens. Those tokens define the primary split groups; they are not guaranteed to be current genus-rank assignments.

The archived protein-family matrix has 2,286 binary columns, following global sequence clustering and a requirement for at least 40 carrier genomes. The original 981-family Highly Pure and 476-family historically Model-Supported inventories remain fixed. No sensitivity here rewrites their membership, rankings or the downstream candidate score.

The recorded saved-output revision is `ac2e43a6ff605196027d062063cc8903b5d5c964`, whose author and commit timestamps are 12 January 2026 at 22:48:37 −0600, equivalent to 13 January 2026 at 04:48:37 UTC. Its tree includes `unified_pipeline_clean/unified_pipeline_run_public/gene_family_matrix.csv` and `complete_genomes_with_proteins.csv`. The recorded original-code revision is `73bf6c9c8e061603e872a703a26d940cab910fb4`. These identifiers and file hashes are recorded in [input_manifest.json](../models/input_manifest.json).

The saved-output timestamp is archive provenance, not a verified download date or proof of the public status of the assemblies at collection. All 426 archived `annotation_date` cells are blank, and an original retrieval timestamp was not captured. The archived extractor also did not record assembly status and looked for notes at the record top level, whereas the dated public schema places relevant notes in nested records. Blank archived notes therefore cannot establish that warnings were absent when the data were collected.

The 416-assembly overlap with a related 449-genome atlas belongs to a separate descriptive cross-atlas comparison. It is not another filtering stage in the cohort sequence below.

## Shared modeling and metric definitions

The four classifiers retain the original parameters, with no tuning. Exact carrier-purity eligibility (10% < positive carrier fraction < 90%) and breadth in at least three groups are fitted using training rows only. Logistic-regression scaling is fitted on training data only. Each fit uses one CPU/BLAS thread. The original 2,286-family vocabulary and sequence clusters remain fixed on the 426-panel in every arm.

The principal metric is the unweighted mean of five held-out-fold ROC AUCs. The dispersion shown below is the sample standard deviation of those five values, not a confidence interval. Saved summaries also contain fold and pooled average precision, pooled genome-weighted AUC, and genus-weighted metrics. These answer different questions, especially with unequal fold sizes and different probability scales across fitted models; they are not interchangeable estimates. [Metric definitions](../models/metric_interpretation.txt) specify the weighting and population conventions.

All analyses are descriptive sensitivities of the fixed matrix. They are not independently collected external validations, whole-genome deduplication, or an end-to-end rebuild from raw proteins. Changing a cohort and repartitioning it can alter several mechanisms jointly. No inferential p value is assigned to these comparisons.

## Assembly status on 30 September 2026

NCBI Datasets reports were retrieved for each exact accession version on 30 September 2026 UTC; historical versions were requested where needed, without silently substituting newer assemblies. Nine batch snapshots supplied 422 reports and a historical-version request resolved the remaining four. The reconciliation has no unresolved record.

| Dated status | Assemblies | nifHDK positive | Comparison |
|---|---:|---:|---:|
| Current | 401 | 106 | 295 |
| Suppressed | 22 | 6 | 16 |
| Previous version | 3 | 0 | 3 |
| Total | 426 | 112 | 314 |

Suppression reasons are heterogeneous; they must not all be labeled contamination. Previous versions are not demonstrated quality failures. The [exact-version status and reason table](../models/inputs/public_status/cohort_status_reconciled.csv), [labeled reconciliation](../models/reference/current_status/frozen_current_status_with_labels.csv), raw JSON responses and [retrieval manifest](../models/inputs/public_status/retrieval_manifest.json) retain the individual records and provenance. Current status is not backdated to original collection.

Two separate exclusions preserve surviving genomes’ original seed-42 folds and stored labels. Removing the 22 suppressed versions leaves 404 assemblies (106 positive, 298 comparison); restricting to current versions leaves 401 (106 positive, 295 comparison). Both retain 70 stored genera and refit training-only feature filters and models. Their full [design and interpretation](../models/analysis/model_current_status_sensitivity_audit.txt) and [summary](../models/reference/current_status/summary.csv) are retained.

## Retained-profile identity and grouping

Two negative-labeled records, Synechococcus sp. MIT S9506 (`GCF_047302775.1`) and Prochlorococcus marinus str. SS51 (`GCF_047302855.1`), occupy original folds 4 and 2. They have identical 178-family retained binary profiles, 179 retained protein-accession/family rows, nif-hit summaries, and saved length/GC statistics (1,627,281 bp and 31% GC). This raises a train/test independence concern but does not demonstrate identity of their complete genomic sequences. Their BioSamples differ. In the dated snapshot MIT S9506 is suppressed with a contamination reason, while SS51 is current.

Two single-row influence checks remove one record at a time while retaining all surviving fold assignments. Removing SS51 is an influence check, not a conclusion that its assembly is invalid. The [pair evidence and analysis](../models/analysis/model_sample_identity_sensitivity_audit.txt) and [model summaries](../models/reference/sample_identity/summary.csv) preserve both tests.

A separate conservative grouping analysis connects all stored Synechococcus, Prochlorococcus and [Synechococcus] labels; Leptolyngbya and [Leptolyngbya]; and Phormidium and [Phormidium]. Its 71 connected groups are repartitioned using the original deterministic shuffled-round-robin procedure with seed 42 and used for the training breadth criterion. Merging entire genera is deliberately broader than resolving the flagged pair. It changes many fold assignments, including one test fold containing 239 of 426 assemblies. The mean AUCs of 0.820–0.852 cannot be attributed solely to the pair, bracketed labels or any one taxonomic error. Bracketed labels are not asserted to be verified taxonomic synonyms.

## Official genus-rank assignments

Dated taxonomy records for 394 unique taxon IDs resolve 422 assemblies into 73 official genus-rank IDs. Four positive-labeled records lack a genus rank: `GCF_029919255.1`, `GCF_003574135.1`, `GCF_002163975.1` and `GCF_000829235.1`. No genus is fabricated from a name token. The all-resolved arm has 108 positive and 314 comparison assemblies. Intersecting current status with genus resolution leaves 398 assemblies (103 positive, 295 comparison) in 70 official genera. One of the four unresolved records was already suppressed, so only three further records are removed from the 401-current cohort.

For each cohort, integer genus IDs are sorted numerically, shuffled with `np.random.default_rng(42)` and allocated round-robin to five folds. This explicitly differs from the original first-occurrence ordering of stored labels. The training breadth filter counts official genus IDs; train and test genus IDs are disjoint. Three official genera span multiple folds under the original stored-label partition. Official taxonomy still assigns the flagged MIT S9506/SS51 pair to different genera (1129 and 1218); they remain separated in the all-resolved arm. The current-and-resolved arm excludes MIT S9506 because of its status.

These arms change cohort, grouping and folds together, so they do not isolate a taxonomy-only effect. The [taxonomy reconciliation](../models/inputs/taxonomy_genus_audit.csv), [raw retrieval manifest](../models/inputs/taxonomy_public/retrieval_manifest.json), [full analysis](../models/analysis/model_official_taxonomy_sensitivity_audit.txt) and [all metric views](../models/reference/official_taxonomy/summary.csv) retain the details. In the current-and-resolved arm, AUC means weighted equally by official genus within each test fold range from 0.758 to 0.796 for the four classifiers, lower than the corresponding unweighted fold means.

## Training-set carrier prevalence

An additional analysis retains all 426 assemblies and original five folds but requires at least 40 carriers within each training fold, alongside the same purity and breadth filters. It retains 449–555 families per fold versus 568–653 in the primary procedure. Here 449 is a family count and has no connection to the related atlas’s genome count.

This nested check is identifiable from the saved matrix because a family with fewer than 40 carriers globally cannot have 40 in a training subset. It still assumes the original global sequence clusters. An absolute threshold of 40 is also stricter in relative terms in the smaller training cohorts (about 10.1–14.2%, versus 9.39% of the full panel). It changes eligibility and does not isolate a causal effect of unsupervised global selection. No inventory reranking is performed. [Methods and masks](../models/analysis/model_training_prevalence_sensitivity_audit.txt), [summary](../models/reference/training_prevalence/summary.csv), and [fold eligibility counts](../models/reference/training_prevalence/prevalence_filter_counts.csv) are retained.

## Model-specific results

Values are mean fold ROC AUC ± sample SD; calculations below are read directly from the saved reference CSVs.

| Analysis | Assemblies | Random forest | Gradient boosting | Logistic regression | XGBoost |
|---|---:|---:|---:|---:|---:|
| Primary panel | 426 | 0.853 ± 0.181 | 0.868 ± 0.124 | 0.878 ± 0.095 | 0.893 ± 0.078 |
| Exclude suppressed | 404 | 0.808 ± 0.271 | 0.877 ± 0.074 | 0.859 ± 0.095 | 0.877 ± 0.099 |
| Current versions | 401 | 0.806 ± 0.260 | 0.874 ± 0.104 | 0.864 ± 0.098 | 0.881 ± 0.101 |
| Remove MIT S9506 | 425 | 0.862 ± 0.173 | 0.866 ± 0.111 | 0.878 ± 0.094 | 0.865 ± 0.115 |
| Remove SS51 | 425 | 0.861 ± 0.173 | 0.876 ± 0.118 | 0.879 ± 0.094 | 0.872 ± 0.117 |
| Conservative connected groups | 426 | 0.820 ± 0.121 | 0.831 ± 0.090 | 0.852 ± 0.096 | 0.844 ± 0.106 |
| Official genus resolved | 422 | 0.890 ± 0.124 | 0.821 ± 0.120 | 0.873 ± 0.100 | 0.871 ± 0.118 |
| Current and genus resolved | 398 | 0.905 ± 0.029 | 0.804 ± 0.151 | 0.844 ± 0.109 | 0.869 ± 0.082 |
| At least 40 training carriers | 426 | 0.808 ± 0.268 | 0.883 ± 0.093 | 0.865 ± 0.127 | 0.891 ± 0.083 |

The two one-predictor logistic-regression baselines use identical folds and unchanged settings within each arm. Their mean fold AUCs are:

| Analysis | Assembly length only LR | Retained family count only LR |
|---|---:|---:|
| Primary panel | 0.816 | 0.749 |
| Exclude suppressed | 0.824 | 0.739 |
| Current versions | 0.822 | 0.740 |
| Remove MIT S9506 | 0.816 | 0.748 |
| Remove SS51 | 0.816 | 0.748 |
| Conservative connected groups | 0.736 | 0.749 |
| Official genus resolved | 0.753 | 0.813 |
| Current and genus resolved | 0.820 | 0.808 |

The 401-current RF mean falls by 0.046904 (0.852807 to 0.805903). Gradient boosting rises by 0.006302, logistic regression falls by 0.014469 and XGBoost falls by 0.012251. Thus a blanket description such as “the result moved by 0.01” is inaccurate. Family-content LR remains above the genome-length-only baseline in this arm (0.863973 versus 0.821722), while RF is below that baseline. Complete fold predictions and the [same-surviving-row original-prediction reference](../models/reference/current_status/original_OOF_same_rows_reference.csv) distinguish evaluation-population change from refitting. The [paired differences](../models/reference/current_status/sensitivity_differences.csv) retain all models and metric views.

The training-prevalence change is similarly model-dependent: RF decreases by about 0.045, GB increases by about 0.016, LR decreases by about 0.014 and XGBoost decreases by about 0.002. These sensitivities qualify the primary estimate without selecting a preferred cohort or partition after inspecting performance.

## Complete outputs and replay

The portable bundle’s `reference/seed42`, `reference/current_status`, `reference/sample_identity`, `reference/official_taxonomy` and `reference/training_prevalence` directories retain summaries, fold metrics, out-of-fold predictions, assignments, selected-family masks, scaler/parameter exports where applicable, and independent verification reports. Raw status/taxonomy inputs and hashes are under `models/inputs`. The [bundle README](../models/README.md), [source/input manifest](../models/input_manifest.json) and [validation record](../models/validation/VALIDATION.md) describe the software versions and checks.

From `models`, reference results can be verified without model fitting:

```sh
python run_audit.py --mode verify --out reference
```

After installing the pinned dependencies as documented in that README, fresh runs can be written to a separate directory:

```sh
python run_audit.py --mode primary --out rerun_outputs
python run_audit.py --mode sensitivities --out rerun_outputs
python run_training_prevalence.py --out extra_prevalence_outputs
```

The primary command performs 20 fits, the three sensitivity modes perform 240 fits, and the separate training-prevalence command adds 20. No live NCBI API is needed to replay the dated inputs. Original-runtime fresh-copy replay and independent metric/mask/prediction checks are recorded in the bundle. An isolated backend reproduced masks, decisions, reported AUC/AP and importance ranks but failed the strict probability tolerance for 140 LR predictions (maximum difference 4.67 × 10⁻¹⁰). That failure is retained rather than silently relaxing tolerance; full cross-backend replay of every mode is not claimed. The supplied reference files preserve the audited results; this documentation revision does not constitute a newly fitted analysis.
