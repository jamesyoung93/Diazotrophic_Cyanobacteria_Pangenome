# Candidate annotation and prioritization

The current shortlist contains **167 families: 22 Tier A, 140 Tier B and 5 high-purity**; 65 have at least one active-phase-up annotation. The complete inventory contains 476 residual-association and 981 high-purity families. Three ubiquitous residual-association controls are excluded from prioritization.

## Results

- [Candidate workbook](results/current167/candidate_families.xlsx): formatted tables for review and experimental selection.
- [Current shortlist](results/current167/corrected_TopFamilies.csv) and [complete inventory](results/current167/corrected_AllComposite.csv): primary machine-readable results.
- [Score definition and counts](results/current167/score_revision_summary.json): components, exclusions and comparison statistics.
- [Mapping sensitivity](results/mapping_sensitivity/summary.json): direct-only condensate analysis gives 166 shortlisted families; all 22 Tier A memberships remain unchanged.
- [Rank sensitivity](results/rank_sensitivity/): model-rank and coefficient-sign summaries for the current 22/140 residual-association shortlist.
- [Figure source tables](results/figures/): values, denominators and intervals for Figures 3 and 5.
- [Regression diagnostics](results/diagnostics/independent_regressions.json): detection-opportunity and prevalence adjustments.

## Replay

Use Python 3.10.16 with `requirements-scientific.txt`. The repository-level [check runner](../scripts/run_release_checks.py) writes a separate output directory and preserves distributed references. To run only this analysis, from the repository root:

```sh
python manuscript_release_2026_10/annotations/run_all.py --out /absolute/path/to/new_annotation_replay
```

The replay checks input hashes, reconstructs annotations from normalized measurements and membership tables, calculates scores and tiers, repeats confound and mapping analyses, summarizes saved model ranks, and regenerates figure sources. It compares results with frozen `expected/` files and reports `PASS` only if checks succeed. It does not refit prediction models or regenerate the formatted workbook.

## Methods and scope

The score omits the incompatible cross-atlas identifier component and its imported HGT penalty, retains the high-purity inventory offset, and excludes the three matrix-verified ubiquitous controls before tier assignment. Its model term is a flat inventory-membership weight, not model importance. Other weights and thresholds are fixed.

Detection/prevalence regressions use the stated baseline inventory and, where explicitly labeled, the baseline 31-family Tier A group. Figure 3 displays all 22 current Tier A families. Figure 5 retains all 476 residual-association families in diagnostic strata and uses a separate 16/22 product-name-pattern reference. These denominators are not interchangeable. The 22.3% proteomics value is a mean of family-level active-up/mapped-study ratios, not an independent per-study probability.

The source workbook, normalized measurements, original scoring code used by reconstruction, and fixed expected results are analysis inputs. They remain because current validation depends on them. Replay can regenerate intermediate comparison tables; the distributed result folder presents current candidate tables and relevant diagnostics.

[Reproduction boundaries](REPRODUCTION_BOUNDARIES.md) describe missing sequence evidence, mapping uncertainty, fixed clustering and environment limits. [Functional annotation](../docs/FUNCTIONAL_ANNOTATION.md) distinguishes product names, keyword bins, context labels and pattern matches. Machine-readable tier strings retain their source wording and do not establish convergent biological evidence.

The [terminology guide](../docs/TERMINOLOGY.md) maps the manuscript labels to the unchanged source identifiers in the CSV files and candidate workbook. Display names do not alter inventory membership, scores, ranks or tier assignments.
