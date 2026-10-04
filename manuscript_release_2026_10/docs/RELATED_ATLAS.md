# Related protein-family atlas comparison

The related analysis is supplied by the ERISE-BNERC repository [Predictive Comparative Genomics for Diazotrophy](https://github.com/erise-bnerc/predictive-comparative-genomics-for-diazotropy/tree/85277c79a0639129fada046ff2638660603dee89), revision `85277c79a0639129fada046ff2638660603dee89` (9 June 2026). The [retrieval manifest](../tables/related_atlas/public_source_manifest.json) records the exact files, retrieval dates and hashes used for the source checks.

## Cohorts and prediction summaries

| Quantity | Primary analysis | Related analysis |
|---|---|---|
| Assemblies | 426 | 449 |
| Retained protein families | 2,286 | 2,551 |
| Candidate inventories | 476 residual-association entries, including three controls, and 981 high-purity families | 503 source-labeled Tier 1 candidate families |
| Evaluation | Stored genus labels blocked; feature filters fitted within training folds | Stratified genus-blocked evaluation |
| XGBoost fold-mean ROC AUC ± sample SD | 0.893 ± 0.078 | 0.968 ± 0.028 |
| Related-analysis mean average precision | — | 0.972 |

The related analysis shares 416 assemblies with the primary panel. These results are not a controlled performance comparison: inclusion rules, feature spaces, preprocessing and split design differ. Its reported AUC is therefore retained here as source context rather than placed beside the primary AUC in the main manuscript table.

The source [stratified-CV script](https://github.com/erise-bnerc/predictive-comparative-genomics-for-diazotropy/blob/85277c79a0639129fada046ff2638660603dee89/new_analytical_scripts/task5_stratified_cv.py), [fold metrics](https://github.com/erise-bnerc/predictive-comparative-genomics-for-diazotropy/blob/85277c79a0639129fada046ff2638660603dee89/new_analytical_task_results/task5_stratified_cv_fold_metrics.csv) and [summary](https://github.com/erise-bnerc/predictive-comparative-genomics-for-diazotropy/blob/85277c79a0639129fada046ff2638660603dee89/new_analytical_task_results/task5_stratified_cv_summary.csv) were inspected, and summary arithmetic was checked against the five folds. The underlying matrix and prediction-level outputs were unavailable for replay; the source-reported performance was not independently reproduced from those inputs.

## Product-label comparison

Normalized product labels connect 442 of the 503 related-atlas candidates to at least one family in the primary 1,457-family inventory. Of these, 363 have informative labels after excluding generic product names. Among the 442 matched related-atlas families, 402 match at least one residual-association entry and 153 match at least one high-purity family; 113 match both inventories. These counts use related-atlas families as their unit, and product-label matching does not establish one-to-one orthology.

The [product-overlap summary](../tables/related_atlas/aryal_current_product_overlap_summary.csv) supplies the counts. The archived comparison tables retain their source names and values. Statements in the historical [comparison note](../tables/related_atlas/aryal_cv_split_overlap_audit.md) about multi-layer evidence are superseded by the interpretation here: shared genomes, incompatible family identifiers and non-independent annotations preclude independent validation. A product-label match contributes no score.

## Identifier and score checks

All 309 pairs in the earlier raw-GF-identifier join had different normalized products. For example, primary DnaK GF_00614 was paired with a DAHP synthase in the related atlas. Current prioritization omits both the resulting atlas bonus and an imported HGT-passenger penalty. A protein-membership or sequence-based crosswalk is required before transferring family-level evidence. The [analysis provenance note](ANALYSIS_PROVENANCE.md) records the score changes.
