# Candidate inventories and data fields

The manuscript uses descriptive names for the two inventories. Machine-readable identifiers remain unchanged so that the published inputs, workbook and reproduction scripts continue to agree.

| Manuscript term | Data identifier | Definition |
|---|---|---|
| Residual-association inventory | `candidate_set = Model-Supported`; source inventory `Tier 1` | 473 families with an estimable positive carrier/non-carrier contrast after purity filtering, plus three ubiquitous controls, for 476 entries. The contrast is computed across the full panel and is not a population-structure-adjusted significance test. |
| High-purity inventory | `candidate_set = Highly Pure`; source inventory `Tier 2` | 981 families with nifHDK-positive carrier fraction at least 0.90. Carrier purity does not measure coverage of diazotroph genomes. |
| Tier A | `priority_tier` beginning `Tier A:` | Eligible residual-association families with adjusted score at least 6.0; 22 families. |
| Tier B | `priority_tier` beginning `Tier B:` | Remaining eligible residual-association families with adjusted score at least 4.5; 140 families. |
| High-purity tier | `priority_tier` beginning `Tier HP:` | Eligible high-purity families with adjusted score at least 4.5; five families. |

The three ubiquitous controls are GF_01297, GF_01899 and GF_01945. They have no non-carrier group and therefore no estimable carrier/non-carrier association. They remain in complete-inventory diagnostics but are excluded from shortlist eligibility. The experimental shortlist has 167 families, of which 65 have at least one active-phase-up proteomics annotation.

## Prediction, inventory ranking and prioritization

**Full-panel inventory ranking** means the ordering attached to the inventories selected using all 426 genomes. In the original classification workflow, feature filtering preceded splitting; per-model importance tables summarize importance across fold fits. Percentiles of absolute importance are averaged across model summaries. This term does not mean that each model was fitted once on all 426 genomes. Importance orders the residual inventory but does not impose an inclusion threshold.

**Training-fold prediction evaluation** fits purity and genus-breadth filters using each training fold. Held-out labels cannot select features. Its importance summaries are supplied separately from the inventory ranking. The family vocabulary and the initial 40-carrier prevalence filter remain fixed across the panel.

**Prioritization score** combines fixed inventory-membership weights with proteomics, condensate and genus-based morphotype annotations, offsets and penalties. It does not use model importance. Tier labels identify score bands, not independent biological validation. The filenames `corrected_AllComposite.csv` and `corrected_TopFamilies.csv` contain the current scores and selections. Earlier score columns are retained only where needed to trace the calculation.

## Annotation columns

| Column or label | Meaning |
|---|---|
| `module_bin` | Coarse product-keyword category. This field also helps define the housekeeping/lineage penalty and must not be changed merely to improve figure labels. |
| `story_role` | Separate broad grouping used for descriptive functional context and the category co-occurrence analysis. |
| FOX-context labels | Interpretive labels concerning fixation in the presence of oxygen. They are distinct from probabilities from the prior gene-level FOX model and add no composite-score term. |
| `consensus_rank_pct_mean` | Mean model-importance percentile used for the inventory ordering. |
| `scope_adjusted_story_score` | Numeric prioritization score used for current tier assignment. |
| `housekeeping_lineage_flag` | Product/category-based scoring penalty and shortlist exclusion indicator. |
| `ubiquitous_control_flag` | Marks the three controls excluded from shortlist eligibility. |

FOX is defined in [Young, Gu and Zhou, Scientific Reports 2026](https://doi.org/10.1038/s41598-026-41873-w). The [functional annotation note](FUNCTIONAL_ANNOTATION.md) specifies category limitations and figure display rules. Numeric results, row membership, source identifiers and scoring rules are unchanged by this terminology mapping.
