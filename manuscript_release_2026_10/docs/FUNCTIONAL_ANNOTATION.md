# Functional annotation: definitions and limits

The main manuscript describes named products and uses external annotations to prioritize experiments. It does not treat a keyword category as a curated pathway, a co-occurring module or proof of function. The legacy category census is archived here for traceability.

## Annotation fields and display rules

| Field or scheme | Meaning | Appropriate use |
|---|---|---|
| `product` | Released representative product name | Describe a candidate while retaining uncertainty in its annotation. |
| `module_bin` | Legacy keyword bin assigned from product text | Reproduce the original category census; not a curated functional inventory. |
| `story_role` and FOX-context labels | Broad interpretation/context labels from prior annotation steps | Trace the original grouping and exploratory analyses; not biochemical pathway membership or FOX model probabilities. |
| Figure 5C product-name pattern | A fixed case-insensitive regular expression applied to `product` | Describe pattern-match frequency across proteomics study-coverage strata. |
| Figure 3 display group | A fixed assignment of the 22 Tier A products to five reading groups | Organize the figure without changing score-bearing bins, penalties, selection or tiers. |

The `housekeeping_lineage_flag` uses product text and the legacy bin. It contributes to score penalties and prioritization exclusions. These fields and scoring rules remain frozen. A future systematic reannotation must separately evaluate any changes to scores and tier membership; changing the display alone does not reclassify families.

## Figure 3: candidate display groups

The figure retains all 22 Tier A families and their prioritization scores. Its fixed display order is amino-acid biosynthesis; protein quality control; translation/RNA-associated products; tetrapyrrole-associated products; and other products. These groups are manual readings of product labels, not curated pathway assignments or evidence of co-inherited modules. Within each group, rows use decreasing carrier purity, with family ID as a tie-breaker. The groups do not modify `module_bin`, housekeeping flags or tier selection.

Carriers and nifHDK-positive carriers are recounted from the frozen 426-genome matrix and labels. Carrier purity is the fraction of a family's carriers labeled nifHDK-positive. It describes this panel; it does not measure coverage of diazotroph genomes, adjust for genome size or independently validate a candidate. Only that column is shaded. Study counts and contributing organisms retain source-specific response annotations; they are not detection-matched rates. The condensate percentile and prioritization score are shown as separate annotations.

An asterisk identifies GF_01387 in the supplied gene-level top-100 mapping. This descriptive overlap is not independent confirmation: the original gene-feature construction was not replayed. Full product names, display groups, carrier counts, annotations and marker provenance are in the [Figure 3 source table](../annotations/results/figures/Figure3_source.csv). [Release verification](../scripts/verify_release.py) independently checks the carrier counts against the model inputs, membership and scores against the frozen candidate table, the display order, mapping marker and matching figure copies.

## Archived Table 2

[The legacy census](legacy_table2_keyword_counts.csv) contains 476 residual-association and 981 high-purity families. It is reproduced directly from `module_bin` in the [complete current table](../annotations/results/current167/corrected_AllComposite.csv). The table is unsuitable as a functional composition result:

- Of 287 residual-association families labeled Other / unassigned, 247 (86%) lack the generic-name terms used by the diagnostic. Many have named products. Being outside a keyword bin does not mean the protein is uncharacterized.
- The separate carbon-metabolism product-name pattern matches 29 residual-association families, whereas the legacy carbon bin contains four. The earlier prose count of 31 does not reproduce and is not used.
- The substring `fur` also matches inside `sulfur`, placing six sulfur-related products in Regulation, including cysteine desulfurases. The audit exposes this limitation without silently changing a score-dependent field.
- A separate generic/hypothetical product-name pattern flags 74 of 476 residual-association families. This is a name-based diagnostic, not a measurement of biological novelty.

Exact patterns and per-family outputs are retained in [historical_diagnostics.py](../annotations/scripts/historical_diagnostics.py) and its [saved output](../annotations/results/logs/03_historical_diagnostics.txt). These patterns are operational definitions; they do not establish curated functions.

## Figure 5C: product-name matches across study coverage

The pattern is:

```text
clp|chaperon|dnak|groe|grol|dnaj|grpe|release factor|translation|elongation factor|typa|trna|nyn|ribosom|isopropylmalate|threonine synthase|ketol-acid|dihydroxy-acid|aminomutase|amino acid
```

Among the frozen 476 residual-association families, the match counts are 17/203, 21/115, 19/75, 15/49 and 12/34 for zero through four mapped proteomics studies. The current Tier A reference is 16/22 (72.7%); the historical 31-family Tier A had 18/31 matches. The reference and coverage-stratum counts use different denominators.

This gradient describes how product-name matches vary with study coverage. It does not by itself establish detection bias as the cause, curated functional composition, or an explanation of the Tier A share. Tier assignment also uses active-phase-up annotations. The complete source rows and plotted intervals are in [Figure 5 sources](../annotations/results/figures/). The figure uses Wilson intervals for the stratum proportions and a separate descriptive Tier A point.

## Category co-occurrence audit

The exploratory test uses residual-association families with a non-Other `story_role` represented by at least five families, present in the retained matrix. It excludes zero-variance family columns within each evaluated genome set. Pairwise binary-family correlations are averaged within and between roles; the statistic is the within-role mean minus the between-role mean. The null permutes role labels 2,000 times per genome set using a single random-number stream initialized with seed 0 (the full-panel test runs first). The one-sided Monte Carlo p value is `(count(null statistic >= observed statistic) + 1) / 2001`.

For all 426 genomes the difference is approximately -0.002 (p = 0.613). For the 112 nifHDK-positive genomes it is approximately -0.008 (p = 0.817, rounded to 0.82). This test does not detect greater within-role co-occurrence under this grouping. It neither validates the category scheme nor rules out biological interactions or modules under other definitions. Pairwise entries are not independent observations; the role-label permutation is the stated exploratory reference, not a confirmatory test of mechanism.

## Reproduction

The [annotation replay](../annotations/README.md) regenerates the diagnostics and compares them with frozen expected outputs. `scripts/verify_release.py` also recounts the archived census, checks the current tiers and verifies Figure 5's unchanged source tables. This release does not add a new KEGG/eggNOG annotation, infer orthology from a product name, or promote exploratory BCAA analyses to validated findings.
