Top-100 comparison package (Protein-driven vs Gene-driven)
Generated: 2026-01-13

Inputs:
- Protein-driven: Recent_Pangenome_results.zip -> Recent_Pangenome_results/results_tables/Table_Tier1_Top100.tsv
  Enriched with: Recent_Pangenome_results/fox_report/master_family_table.tsv
  Representative protein_id per family from: protein_family_map.tsv
- Gene-driven: random_forest_top_100_features_with_rep.csv
  Annotation source: expanded_categories_annotated.json (modal protein_name across occurrences)

Files:
1) top100_protein_vs_gene_item_level.csv
   - One row per hit (100 Protein-driven + 100 Gene-driven).
   - Columns:
     set_type: "Protein-driven" or "Gene-driven"
     set_rank: 1..100 within each set
     feature_kind: gene_family (protein-driven) vs gene_feature (gene-driven)
     feature_id: GF_XXXXX for gene_family; genome_*_CDS_* for gene_feature
     annotation_primary: product/protein_name used for functional binning
     functional_bin: rule-based high-level bin (see below)
     is_unknown: True if annotation suggests hypothetical/DUF/uncharacterized
     representative_gene: (gene-driven only) locus tag provided as representative
     representative_locus_tag: best-available locus tag
     representative_protein_id: best-available protein accession
     effect_size_mean / consensus_rank_pct_mean / module_bin / fox_status / diazotroph_pct_* / n_member_*: (protein-driven only)
     importance / Rank / n_genome_occurrences_in_json: (gene-driven only)

2) top100_protein_vs_gene_functional_composition.csv
   - Functional-bin counts and percent composition within each top-100 list.

3) top100_protein_vs_gene_diazotrophy_share.csv
   - Share of each list in a predefined "diazotrophy-physiology" bin set.

Functional binning (keyword rules; conservative, single-bin assignment):
- Nitrogenase structural
- Nitrogenase cofactor & maturation
- Heterocyst differentiation & envelope
- Metal/cofactor transport
- Redox & electron transport / respiration
- Photosynthesis & carbon fixation
- Stress response / ROS / proteostasis
- Regulation & signaling
- Translation & ribosome
- DNA replication/repair
- Transporters (non-metal specific)
- Central metabolism (C/N/energy)
- Unknown / uncharacterized
- Other (misc.)

Diazotrophy-physiology bins used in the summary share:
Nitrogenase structural; Nitrogenase cofactor & maturation; Heterocyst differentiation & envelope;
Metal/cofactor transport; Redox & electron transport / respiration; Photosynthesis & carbon fixation;
Stress response / ROS / proteostasis; Regulation & signaling.
