Literature proteomics to protein-family mapping

Mapping is accession/locus based only: old locus tag, RefSeq locus tag, old accession, or WP accession -> NCBI feature/crosswalk bridge -> current genome_protein_family_map.tsv.
Product-name/fuzzy matching is intentionally not used in the mapped family evidence sheets.
Rows are deduplicated to one source protein record before family summarization, and family_evidence_by_study.csv is one row per study + gene_family.
Proteins with a WP bridge but no gene_family are marked mapped_to_wp_not_in_filtered_family_atlas; these are real proteins but absent from the filtered 2,286-family atlas.
Aryal 2011 binary .xls files are recorded in source_map.csv but not parsed here because this runtime lacks a reliable .xls reader.