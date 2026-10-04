# Portable annotation and prioritization audit

This reproducible analysis package preserves the frozen inputs and regenerates the corrected scoring rule, annotation diagnostics, mapping sensitivities and Figures 3/5. Analysis requires no input downloads or prediction-model refits. Historical source text is retained for provenance.

The three result scopes are deliberately separate:

| Scope | Tier A | Tier B | Highly Pure shortlist | Total | At least one active-up annotation |
|---|---:|---:|---:|---:|---:|
| Original released shortlist |31|138|5|174|71|
| Current correction, primary |22|140|5|167|65|
| Direct-only condensate sensitivity |22|139|5|166|64|

All use the frozen inventory of 1,457 families:476 Model-Supported and 981 Highly Pure. The original 174 includes169 Model-Supported families and 5 Highly Pure families; it is not174 Model-Supported families. The intermediate 169-family correction is retained in replay outputs for traceability, but it is not the primary corrected shortlist.

**Run from an extracted copy.** Use Python 3.10.16 with the imported scientific-library versions pinned in `requirements-scientific.txt`:

```text
python -m pip install -r requirements-scientific.txt
python run_all.py --out replay
```

The first command installs dependencies only if needed; the audit itself requires no network. `run_all.py` finds inputs relative to its own location, checks packaged hashes, runs each analysis, compares outputs with pre-bundle reference results, and writes `replay/verification.json`, per-step logs and `replay/OUTPUT_HASHES.json`. A successful run ends with `PASS`. Supplying a new `--out` folder keeps prior results intact. All scripts also expose explicit file arguments through `--help`.

The supplied `results/` folder is one completed replay. The separate extraction test is documented in `notes/fresh_extraction_verification.json`. A second replay in a newly created scientific environment is documented in `notes/clean_environment_verification.json`. All 40 generated CSVs and all seven regression summaries are numerically identical across the two environments. Both figures have identical decoded RGBA pixels; their encoded PNG file hashes differ. These tests ran on Windows.

The initial runtime reports SciPy module 1.13.1 versus distribution metadata 1.13.0, and Pillow module 11.1.0 versus metadata 10.4.0. Both values are preserved in `runtime_versions.json`; pins follow the imported module versions. The clean-environment proof records matching module and distribution versions for every pinned scientific dependency. PDF creation metadata and image encoding can affect file hashes without changing scientific data or rendered pixels.

**What runs, and what each input means.**

1. `scripts/reconstruct_annotations.py` rebuilds study-level active-up classifications from normalized source measurements, collapses unique studies/family, recounts genus-proxy breadth from the genome matrix and metadata, reconstructs condensate mapping/ranks from protein rows, and reproduces the released identifier-based atlas scoring bug. It does not import the production scorer or use workbook outcome columns to calculate these components. It opens the original workbook only after reconstruction for comparison. For the context-count comparison alone, a released blank is compared with a verified empty context set of size 0; scores and other missingness semantics are unchanged.
2. `scripts/score_current.py` removes the unverified cross-atlas component and imported HGT penalty, retains the Highly Pure+1 heuristic offset, and excludes three matrix-verified ubiquitous families before assigning tiers. Other weights/thresholds remain fixed. Numeric scores are retained even for excluded controls. It checks all 1,457 scores against the independent reconstruction and writes released, intermediate 169 and current 167 snapshots plus migrations.
3. `scripts/historical_diagnostics.py` reproduces the original D1–D8 calculations using the original workbook and all 1,457 frozen inventory rows. Historical Tier A refers to 31 families; historical shortlist refers to 174. Its editorial headings were qualified without changing the calculations. The carbon product-name pattern returns 29 families, not the 31 stated in the earlier prose. The untouched original script is in `provenance/`.
4. `scripts/independent_regressions.py` repeats the original regressions with odds-ratio intervals, categorical detection-opportunity sensitivity, common-prevalence-support analysis and exploratory Tier A composition models. These regression outputs use the original 31-family Tier A classification. Current 22-family composition is reported separately in the corrected score summary and Figure 5. A nonsignificant coefficient is not evidence of equivalence; the 22.3% figure is a mean of family-level active-up/mapped-study ratios, not a constant independent chance in every study.
5. `scripts/mapping_sensitivity.py` separately excludes the five selected ambiguous-accession mappings and the 804 cluster-only fallback mappings. It preserves the original 1,311-family percentile reference. A conservative analysis makes unsupported released best annotations missing; a second places the best retained score against the fixed original reference. Neither silently reranks unaffected families. Both direct-only analyses remove only GF_02157 from Tier B. All 22 Tier A best condensate annotations and memberships remain unchanged.
6. `rank_audit/run_rank_audit.py` regenerates original-versus-corrected model rank comparisons and 50-fit coefficient stability from supplied saved importance/coefficients. It does not refit models. `scripts/current_shortlist_rank.py` then selects the exact current 22 A and 140 B families. The rank-input manifest documents source/subset hashes. Repeated partitions are dependent; coefficient signs are conditional on correlated predictors and are not biological effect directions.
7. `scripts/plot_figures.py` regenerates compact Figure 3 with all 22 current Tier A rows, and Figure 5 with confidence intervals and denominators. Figure 5 diagnostic strata retain the frozen 476 Model-Supported inventory, including three ubiquitous controls excluded from prioritization. Its Tier A reference is 16/22 product-pattern matches (72.7%). Figure 3 displays released historical importance; the score itself uses a flat 2.0 Model-Supported/1.0 Highly Pure inventory weight and never uses an importance rank.
8. `scripts/verify_reproduction.py` compares current scores/tiers/migrations, independent reconstruction, historical diagnostic results, regression output, mapping sensitivities, rank tables and figure source tables with reference outputs captured before bundle construction.

**Read these outputs first.**

- `results/current167/corrected_TopFamilies.csv`: primary corrected 167-family shortlist.
- `results/current167/score_revision_summary.json`: exact counts, exclusions and score decomposition, with historical scopes labeled.
- `results/current167/tier_migrations.csv`: changes from original release; `control_migrations_vs169.csv` isolates the ubiquitous-control step.
- `results/figures/`: Figure 3/5 PNG/PDF files, source CSVs and captions.
- `results/mapping_sensitivity/summary.json`: separately labeled 167-versus 166 sensitivity and family-level migration tables.
- `results/rank_sensitivity/`: current 22/140 rank/sign tables and interpretation.
- `results/diagnostics/independent_regressions.json` and `results/logs/03_historical_diagnostics.txt`: original frozen-inventory diagnostics.
- `snapshots/current167_author_review.xlsx`: previously verified, static review workbook. The portable replay regenerates CSVs, not this formatted XLSX; its original saved-file QA is in `notes/`.

`inputs/released/proteomics_composite_family_evidence.xlsx` is a byte-identical original source workbook. `expected/` contains separately saved validation targets; replay never updates these. `PACKAGE_MANIFEST.json` hashes immutable package files, while each replay hashes its outputs. `provenance/input_copy_manifest.json` records source hashes and the genome-metadata reduction to three required columns and 426 unique assemblies. The separate rank-input manifest documents subsets of saved model outputs. Original source code in `provenance/` is explanatory and is not run by `run_all.py`.

The named functional patterns are exploratory product-name annotations, not curated pathway membership. Tier labels preserved in CSVs are historical machine-readable strings; they do not establish convergent biological evidence. Read `REPRODUCTION_BOUNDARIES.md` before interpreting any score, association or apparent cross-atlas match.

Figure terminology: Figure3 calls its name-based exclusion "nif-keyword", replacing "core" to avoid implying every excluded name is a canonical structural nif component. This change follows the recorded clean replay and changes no algorithm, data or image. `notes/caption_terminology_revision.json` records the exact byte-level scope and repeated saved-output verification.

The October presentation revision labels Figure 5C as product-name matches across study coverage. Its regex, scientific data and intervals are unchanged. The [functional annotation note](../docs/FUNCTIONAL_ANNOTATION.md) separates legacy bins, context labels and product-name patterns.
