# Filtering and prioritization provenance

The main manuscript reports the current prediction procedure and prioritization score. This note records implementation differences that explain the retained source files and comparison columns.

## Feature filtering and model ranking

The original discovery calculation applied a label-dependent purity filter before cross-validation. It excluded 1,659 of 2,286 families, leaving 627 before subsequent filtering; 981 excluded families were diazotroph-pure. Exact implementation of the intended ≤0.10 and ≥0.90 boundaries excludes three additional comparison-associated families, leaving 624 globally eligible families. None of those three belongs to either positive candidate inventory.

Prediction evaluation now fits purity and genus-breadth filters only to each training fold, retaining 568–653 families. Test labels and test-genus occurrences do not enter those filters. The 981-family high-purity inventory and 476-entry residual-association inventory remain full-panel discovery outputs. Their original importance ranking is retained separately from the training-fold ranking. The latter has consensus Spearman correlation 0.831 across the 476 entries, with 33 of the original top 50 remaining in the top 50.

## Prioritization score

The earlier score used an invalid raw-identifier join across independently clustered atlases. Current scoring removes the cross-atlas bonus and its imported HGT-passenger flag, retains the original inventory, proteomics, condensate and morphotype weights, and excludes three matrix-verified ubiquitous controls before tier assignment. The high-purity offset remains a baseline weight rather than external evidence.

The current shortlist contains 22 Tier A, 140 Tier B and five high-purity-tier families (167 total), with 65 active-phase-up annotations. Relative to the earlier 174-family shortlist, nine Tier A families move to Tier B, five Tier B families fall below threshold and two additional Tier B families are excluded as ubiquitous controls. No family enters the shortlist. The third ubiquitous control was not in the earlier shortlist. These comparisons are recorded in [score migrations](../annotations/results/current167/score_migrations.csv), [tier migrations](../annotations/results/current167/tier_migrations.csv) and [ubiquitous controls](../annotations/results/current167/ubiquitous_controls.csv).

The [annotation replay](../annotations/README.md) reconstructs the scores from archived measurements and mappings and compares them with fixed reference outputs. The [terminology guide](TERMINOLOGY.md) maps manuscript terms to retained data identifiers. Terminology changes do not change scientific values, classifier fits, candidate membership or score rules.
