# Aryal CV-Split Protein-Level Overlap Audit

## Bottom line

Aryal's claimed starting point is closely related to the current protein-family work, but it is not the identical manuscript universe. The CV split package contains 449 unique genomes, including 33 GCA assemblies, while the current manuscript matrix has 426 GCF RefSeq genomes. The overlap is 416 genomes.

The family-level comparison is more cautious: Aryal's Task 7 output gives a 503-family Tier 1 proxy, while the current manuscript inventory has 476 Model-Supported and 981 Highly Pure families. Direct GF-ID overlap exists, but because the Aryal run uses a 2,551-family namespace and the current matrix uses 2,286 families, non-overlap should be interpreted as version/namespace drift unless the exact 449 x 2,551 matrix is supplied.

## Key counts

- Aryal CV label raw rows: 771 (complete_genomes_labeled-james.csv rows)
- Aryal CV unique labels: 449 (deduplicated assembly_accession)
- Aryal CV duplicated accessions: 237 (duplicates had no diazotroph-label conflicts)
- Aryal CV split unique genomes: 449 (union of train/test split files)
- Current manuscript matrix genomes: 426 (Diazotrophic_Cyanobacteria_Pangenome-main/unified_pipeline_clean/unified_pipeline_run_public/gene_family_matrix.csv)
- Current manuscript matrix families: 2286 (Diazotrophic_Cyanobacteria_Pangenome-main/unified_pipeline_clean/unified_pipeline_run_public/gene_family_matrix.csv)
- CV split genomes also in current matrix: 416 (97.65% of current matrix)
- CV split genomes absent from current matrix: 33 (all are GCA accessions)
- Current matrix genomes absent from CV split: 10 (all are GCF accessions)
- Current Model-Supported families: 476 (Diazotrophic_Cyanobacteria_Pangenome-main/unified_pipeline_clean/unified_pipeline_run_public/fox_report/tier1_positive_model_selected.tsv)
- Current Highly Pure families: 981 (Diazotrophic_Cyanobacteria_Pangenome-main/unified_pipeline_clean/unified_pipeline_run_public/fox_report/tier2_pure_positive_heldout.tsv)
- Aryal Task7/Tier1-proxy families: 503 (https://raw.githubusercontent.com/erise-bnerc/predictive-comparative-genomics-for-diazotropy/main/new_analytical_task_results/task7_hgt_proximity.csv)
- Aryal families found in current matrix columns: 464 ()
- Aryal families found in current inventory: 309 ()
- Aryal families found in current Model-Supported: 110 ()
- Aryal families found in current Highly Pure: 199 ()
- Current Model-Supported families recovered by Aryal Tier1 IDs: 110 (23.11% of current Model-Supported)
- Current inventory families recovered by Aryal Tier1 IDs: 309 (21.21% of current 1457-family inventory)
- Aryal-only Tier1 IDs vs current inventory: 194 (GF IDs absent from current Model-Supported/Highly Pure inventory)
- Aryal-only Tier1 IDs vs current matrix: 39 (GF IDs absent from current 2,286-family matrix)
- Shared Aryal/current-inventory exact product matches: 0 (out of 309 shared GF IDs)
- Shared Aryal/current-inventory FOX-layer families: 31 ()
- Shared Aryal/current-inventory proteomics directional families: 33 ()
- Shared Aryal/current-inventory three-layer families: 3 ()
- Aryal families with exact product match to current inventory: 442 (GF-ID independent product-level crosswalk)
- Aryal informative products with exact current match: 363 (excludes hypothetical/DUF/generic family proteins)
- Aryal product matches touching current Model-Supported: 402 ()
- Aryal product matches touching current Highly Pure: 153 ()
- Aryal product matches touching current FOX layer: 47 ()
- Aryal product matches touching current proteomics directional layer: 100 ()
- Aryal product matches touching current three-layer evidence: 9 ()

## Interpretation

The new light helps: this is not an unrelated old gene-level analysis. It is a protein-level analysis built around a genome set that substantially overlaps the manuscript set. However, it still should not be merged as if it were the final 426-genome manuscript analysis, because the genome set, family count, Tier 1 count, and nif-anchor IDs differ.

The most useful direct comparison is therefore not raw GF-ID overlap alone. Raw GF IDs overlap numerically, but the product labels do not match for shared GF IDs, indicating namespace drift between clustering versions. The stronger biological comparison is product/function-level overlap.

At the product/function level, Aryal's 503-family Tier 1 proxy strongly recovers the same biological territory as the current inventory: 442 Aryal families have an exact product-label match in the current 1,457-family inventory, 402 touch current Model-Supported families, 47 touch the FOX layer, 100 touch directional proteomics, and 9 touch the current three-layer model+FOX+proteomics evidence set.

## Strongest product-level overlaps

- cytochrome b6 (GF_01462 -> GF_01297): Respiration & bioenergetics
- iron-sulfur cluster assembly accessory protein (GF_00584 -> GF_02143): Metallocluster & Fe-S biogenesis
- histidinol dehydrogenase (GF_02394 -> GF_00355): Other / unassigned
- NADP-dependent isocitrate dehydrogenase (GF_00332 -> GF_00624): Other / unassigned
- glycolate oxidase subunit GlcD (GF_01767 -> GF_00314): Respiration & bioenergetics
- glucose-6-phosphate dehydrogenase (GF_00859 -> GF_00486;GF_01818): Other / unassigned
- alpha-ketoacid dehydrogenase subunit beta (GF_01787 -> GF_00043): Other / unassigned
- thioredoxin (GF_00635 -> GF_00098;GF_00352;GF_00965): O2 protection & redox homeostasis
- Fe-S cluster assembly protein SufB (GF_01827 -> GF_00968): Metallocluster & Fe-S biogenesis

## Output files

- evidence_tables/aryal_cv_split_overlap_summary.tsv
- evidence_tables/aryal_cv_split_genome_overlap.csv
- evidence_tables/aryal_task7_tier1_vs_current_inventory_detail.csv
- evidence_tables/aryal_current_shared_family_detail.csv
- evidence_tables/aryal_current_product_overlap_long.csv
- evidence_tables/aryal_current_product_overlap_summary.csv
- evidence_tables/aryal_tier1_not_in_current_inventory.csv
- evidence_tables/current_inventory_not_in_aryal_tier1.csv
- evidence_tables/aryal_overlap_by_module.tsv
- evidence_tables/aryal_vs_current_nif_anchor_check.csv
