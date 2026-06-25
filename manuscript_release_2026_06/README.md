# June 2026 manuscript release

This folder contains the manuscript-facing release layer for the current "Beyond nif" protein-family analysis of cyanobacterial diazotrophy.

The core pangenome/modeling workflow in `../unified_pipeline_clean/` generated the 426-complete-genome protein-family universe and the two non-overlapping candidate inventories:

- `476` Model-Supported families: classifier-supported accessory candidates.
- `981` Highly Pure families: near-diagnostic `nif`/diazotrophy-context families held out from the accessory-biology narrative.

The files in this release layer add the biological interpretation and prioritization overlays used by the current manuscript.

## Important scoring decisions

- External proteomics evidence is scored only when the protein is up in a nitrogen-fixation-relevant active phase or cell state.
- Local Cyanothece proteomics is not used in the current scoring.
- FOX ensemble probability is not used in the current scoring.
- Curated FOX-context labels may appear as annotation context, not as a primary score driver.
- Related-atlas concordance supports interpretation but is capped so it does not dominate the score.
- Condensate-driver rank is retained as a minor, indirect prioritization signal; it is not treated as proof that a family forms or enters a condensate.
- Highly Pure families are not penalized for missing related-atlas support when the related atlas was scoped around Model-Supported families.

## Tables

Core current tables:

- `tables/composite_proteomics_priority/summary_metrics.csv`
- `tables/composite_proteomics_priority/component_definitions.csv`
- `tables/composite_proteomics_priority/module_story_summary.csv`
- `tables/composite_proteomics_priority/top_story_families.csv`
- `tables/composite_proteomics_priority/high_purity_marker_families.csv`
- `tables/composite_proteomics_priority/family_composite_all.csv`

External proteomics mapping summaries:

- `tables/literature_proteomics_family_mapping/mapping_summary.csv`
- `tables/literature_proteomics_family_mapping/source_map.csv`
- `tables/literature_proteomics_family_mapping/family_summary.csv`
- `tables/literature_proteomics_family_mapping/family_evidence_by_study.csv`

Related-atlas and robustness overlays:

- `tables/related_atlas/aryal_current_shared_family_detail.csv`
- `tables/related_atlas/aryal_current_product_overlap_summary.csv`
- `tables/related_atlas/aryal_cv_split_overlap_summary.tsv`
- `tables/related_atlas/aryal_cv_split_overlap_audit.md`
- `tables/related_atlas/aryal_task7_hgt_proximity_from_github.csv`

Condensate-driver overlays:

- `tables/condensate/condensate_ranking_enrichment_summary.tsv`
- `tables/condensate/condensate_ranking_family_overlay.tsv`
- `tables/condensate/condensate_ranking_mapping_summary.tsv`
- `tables/condensate/condensate_ranking_top_model_supported.tsv`

Phase-alignment and narrative audits:

- `tables/phase_alignment/active_phase_up_top_vs_nontop_summary.csv`
- `tables/phase_alignment/narrative_active_phase_up_impact_summary.csv`
- `tables/phase_alignment/phase_alignment_study_summary.csv`
- `tables/phase_alignment/phase_alignment_narrative_summary.csv`

## Scripts

The `scripts/` folder contains the current overlay builders and audits used to generate the release-layer tables from local/source inputs:

- `build_literature_proteomics_family_mapping.py`
- `build_proteomics_composite_evidence.py`
- `analyze_condensate_ranking_overlay.py`
- `audit_aryal_cv_split_overlap.py`
- `audit_nfix_phase_alignment.py`
- `compare_composite_iteration_to_baseline.py`

These scripts are intended as manuscript-release provenance for the overlay analyses. The core pangenome construction and modeling pipeline remains in `../unified_pipeline_clean/`.

## Figure assets

`figures/` contains the editable current Figure 1 architecture slide and a rendered PNG preview.

## Notes for reviewers

The release layer is intentionally additive. It does not replace the original pangenome workflow; it documents the additional biological prioritization analyses now used in the manuscript.