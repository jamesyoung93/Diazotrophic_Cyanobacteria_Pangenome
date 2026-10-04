# Portable model audit v2

For the October repository release, the [release check runner](../scripts/run_release_checks.py)
verifies a copy of the saved model outputs, preserving the frozen references and
their hashes. The direct `verify` command below writes verification reports under
`--out`; use a separate copy when preserving the distributed files is required.

This folder is self-contained. Copy or unzip it anywhere, change into this folder,
and run the commands below. No parent workspace, user-specific path, network API,
large protein database, spreadsheet, or manuscript is needed after installing the
listed Python dependencies. Inputs and source are checked against bundled SHA256s.

## Run

Use Python 3.10.16 with the exact versions in requirements.txt. In a clean environment:

```sh
python -m pip install -r requirements.txt
python run_audit.py --mode primary --out rerun_outputs
python run_audit.py --mode repeated --out rerun_outputs
python run_audit.py --mode sample_identity --out rerun_outputs
python run_audit.py --mode current_status --out rerun_outputs
python run_audit.py --mode official_taxonomy --out rerun_outputs
python run_audit.py --mode verify --out rerun_outputs
```

Alternatively, `python run_audit.py --mode all --out rerun_outputs` runs all five
analyses and verification (460 bounded fits). `primary` is an alias for `seed42`:
four classifiers over five genus folds (20 fits), writing to `rerun_outputs/seed42`.
The repeated run fits four fixed LR specifications over seeds 42-51 (200 fits).
Sample identity uses four scenarios and six models (120 fits); current status
and official taxonomy each use two scenarios and six models (60 fits each).
`--mode sensitivities` runs only those last three analyses (240 fits).
Every fit uses one CPU/BLAS thread. `--out` can be any separate output directory.
Runtime depends on hardware; no training, download, or installation is needed merely
to inspect reference results. To verify saved reference results without fitting:

```sh
python run_audit.py --mode verify --out reference
```

The last command requires the normalized reference columns supplied in this bundle.
The included fresh-copy validation results record the actual reproducibility check.

The existing reference runtime reports SciPy module version 1.13.1 but installed
distribution metadata 1.13.0. This disagreement is disclosed rather than treated
as a clean installation. The reference runtime used threadpoolctl 3.5.0. The v2
installation specification uses threadpoolctl 3.6.0 because its documented Windows
long-library-path fix addresses the suspected cause of a crash observed in a deeply nested isolated virtual
environment; joblib 1.4.2 is explicitly pinned. No scientific model parameters or
comparison tolerances change. See validation evidence for the scope and outcome
of clean-environment testing, separate from same-existing-runtime replay.
Upstream fix: https://github.com/joblib/threadpoolctl/pull/189

The separate isolated primary run completed all 20 fits, but **strict probability
replay failed** for 140 LR predictions (maximum difference 4.67e-10) at the unchanged
rtol1e-10/atol1e-12. Exact masks, class decisions, reported AUC/AP and importance
ranks agree; no tolerance was relaxed. The default CLI preserves that assertion
and may therefore exit nonzero after saving computed results in this backend.
See [the clean-environment comparison](validation/clean_environment/REPORT.txt).
This is distinct from exact fresh-copy replay in the original runtime. Only the
primary run was fitted in the isolated environment; full cross-backend replay of
all modes is not claimed. The saved comparison can be regenerated without fitting:

```sh
python compare_clean_backend.py --reference reference/seed42 --clean validation/clean_environment/seed42 --out cross_backend_rerun
```

Its JSON explicitly retains `strict_probability_comparison: FAIL`; add
`--assert-strict` to reproduce the strict failure as a nonzero exit.

## Inputs and source

`inputs/gene_family_matrix.csv` is the unchanged 426-by-2286 binary family matrix.
`inputs/complete_genomes_with_proteins.csv` is the unchanged released metadata and
contains genus, diazotroph label and genuine total_ungapped_length. The genome-length
column agrees with `inputs/assembly_quality.tsv` for all 426 accessions. No genuine
total-proteome protein count was available. Retained-family protein-map rows were
not treated as total proteins.

`source/04_classify_original.py` and `source/04_classify_proposed.py` are the original
and review-only proposed source; a unified diff is also provided. The CLI extracts
the required pure functions from the proposed source AST, avoiding its legacy main
inventory writer and optional plotting dependencies. The seed42 run invokes its
actual run_genus_cv function. The repeated runner uses an independently checked,
vectorized equivalent exact-count mask. Both use the unchanged original custom
splitter: shuffle the 75 unique stored genus tokens in their first-occurrence
order and allocate round-robin to five folds. This is not sklearn GroupKFold
and is not stratified by class. These are literal stored labels, not a guarantee
of disjoint current official genus-rank taxonomy.

Purity and minimum-three-genus eligibility use training rows only. The exact purity
cutoffs retain 10% < positive carrier fraction < 90%. The lower inclusive exclusion
uses integer counts, fixing the original 1-0.9 floating-boundary issue. StandardScaler
is fit on training rows only. Family clustering and the global >=40-genome vocabulary
remain fixed on the full panel: this is a fixed-matrix audit, not full end-to-end
external validation. No hyperparameter retuning or seed selection was done.

## Results and interpretation

Reference seed42 mean fold ROC-AUC (fold SD): RF 0.8528(0.1810), GB 0.8675(0.1242),
LR 0.8784(0.0949), XGB 0.8931(0.0783). Equal-genus-weighted pooled AUC is 0.7823-0.8154,
and genome-weighted pooled AUC is 0.8564-0.9211. These are different estimands; see
metric_interpretation.txt for exact definitions, genus/class coverage and AP caveats.
Fold SD is not a confidence interval.

The repeated analyses compare exact-boundary family LR, assembly-length-only LR,
retained-family-count-only LR, and a separately labeled joint family-plus-length LR.
All use identical partitions and LR defaults. Across seeds 42-51, family LR exceeds
length-only LR in mean fold AUC by 0.0185-0.1469 (mean 0.0783). Adding length to family
LR changes mean fold AUC by -0.0015 to +0.0007. These dependent partition ranges are
descriptive, not independent experiments or inferential significance tests.

`reference/seed42` and `reference/repeated` hold fold metrics, OOF predictions,
selected family lists, summaries and paired differences. The CLI independently
recomputes AUC by positive-negative pair comparison and AP by threshold increments,
checks class labels and exact split membership, reaggregates fold summaries, and
compares every OOF probability with its reference. Contract checks cover unchanged
splits, exact 10%/90% boundaries and held-out-label invariance of train feature masks.
Reference comparison tolerance is rtol 1e-10/atol 1e-12; the fresh-copy test uses the
recorded Windows environment. Version pinning does not guarantee bitwise results
on every platform/compiler; any mismatch is reported rather than silently accepted.

## Dated identity, assembly-status and official-genus sensitivities

`sample_identity` replays the control, excludes either of two flagged negative
records while retaining original folds, and separately connects explicitly named
genus labels before repartitioning. The pair has identical retained family vectors
and retained protein/family pairs plus matching assembly summary metadata. Complete
genomic-sequence identity was not established. The reduced protein-map input holds
only the 358 relevant rows; its full-source and subset hashes are documented.
One pair member, MIT S9506, is currently suppressed for contamination. Its exclusion
is a suppressed-record sensitivity; excluding SS51 is a paired influence check.

`current_status` removes 22 currently suppressed exact assembly versions (404
survive), then separately retains only 401 current versions (also excludes three
previous versions). The public snapshot has 401 current, 22 suppressed and three
previous records, with no unresolved accessions. Suppressed records have varied
reasons; previous versions are not designated quality failures. Original folds
remain fixed. Current status is not backdated to original collection.

`official_taxonomy` excludes four unresolved genus-rank records (all positive),
then groups 422 records by 73 official genus tax IDs. A second arm intersects
current status and resolved genus (398 records, 70 official genera). Within each
arm, integer tax IDs are sorted numerically before the seed-42 shuffle and
round-robin split. That declared order differs from original first-occurrence
stored-label order. The training breadth filter counts official genus IDs. These
arms change population, grouping and partition together and cannot isolate a
taxonomy-only or leakage-only effect.

Each sensitivity uses all four classifiers and length-only/count-only LR baselines,
exact training-only filters and unchanged parameters. Full fold metrics, predictions,
feature masks, class/group counts, coefficients/scalers and independent checks are
in corresponding `reference` subfolders. Both stored-genus and audit-group weightings
are retained; for official taxonomy, `audit_group_balanced` means current genus IDs.
The same-surviving-row original-OOF comparison retains original partitions and
cannot isolate a repartitioning mechanism. No inferential p values are assigned.

Raw public status and taxonomy responses, retrieval manifests, reconciled CSVs and
SHA256 hashes are included under `inputs`. No live API is needed to replay these
dated analyses. Interpretation and limitations are in the three `analysis/*audit.txt`
reports. A future taxonomy/status snapshot would be a new analysis, not a replacement
of these frozen inputs.

## Additional training-prevalence sensitivity (separate 20 fits)

```sh
python run_training_prevalence.py --out extra_prevalence_outputs
```

This optional entry point adds an absolute requirement of at least 40 carriers in
each training fold to the same purity/breadth filters, on the original 426 rows and
seed-42 folds. Its masks contain 449–555 families; primary masks contain 568–653.
Mean AUCs are RF0.8078, GB0.8834, LR0.8646 and XGB0.8908. This is nested and
identifiable from the fixed matrix because a globally<40 family cannot reach40
training carriers, but the relative-prevalence severity increases in smaller
training cohorts and global clustering remains fixed. It is a separate sensitivity,
not a replacement primary estimate or isolated test of selection bias.

Its source, reference outputs and interpretation are under `analysis` and
`reference/training_prevalence`. The script checks exact masks/assignments,
independently reconstructs metrics and LR predictions, and enforces the unchanged
reference-probability tolerance. It is **additional** to the main `all` command's
460 fits; that command and its original five analyses remain unchanged.

## Frozen inventory

The existing 476-family Model-Supported inventory and 981-family Highly Pure positive
inventory, including released rankings, were frozen. No command in this bundle
invokes their writers or substitutes new importances. Any corrected CV importances
are named AUDIT_ONLY. frozen_inventory_manifest.json records original release hashes.
The three global exact 10% boundary families have Negative direction and belong to
neither positive inventory. Full-data non-pure eligibility changes 627->624; corrected
prediction uses fold-specific masks, not a fixed 624-column training matrix.

The frozen 476 set contains 473 estimable positive full-panel carrier/noncarrier
label-frequency contrasts plus three ubiquitous families whose absent-group
contrast was undefined but treated as positive by the source default. Membership
in importance exports admits zero values and does not further filter that set;
importance percentiles determine ordering. Those three ubiquitous edge cases are
different from the three negative exact-10% cutoff families above. They are retained
for resource traceability and are not newly validated candidates.

Do not run the review-only proposed source main routine over a released resource:
its future model importances would constitute a separately versioned ranking update.
The supported entry point is run_audit.py, which writes only audit artifacts.

Fresh-copy execution evidence: [validation/VALIDATION.md](validation/VALIDATION.md).
