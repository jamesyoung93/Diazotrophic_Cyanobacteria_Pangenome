# October 2026 release notes

The manuscript now distinguishes family-content prediction from exploratory annotation. Detailed cohort and annotation audits are available in the repository so the main text can state the central results concisely.

## Scientific scope

- Preserve the 426-genome matrix, 476/981 inventories, saved predictions and reference outputs.
- Supply the corrected 167-family shortlist: remove unverified cross-atlas identifier evidence and its imported HGT penalty, retain the Highly Pure inventory offset, and exclude three ubiquitous controls from prioritization. Full scoring details and migrations are in the annotation bundle.
- Retain all 22 Tier A rows in Figure 3. Label Figure 5C as product-name-pattern matches across study-coverage strata, not curated functional composition.
- Archive the legacy Table 2 keyword census and category co-occurrence result with their limitations. No new ontology annotation or pathway-membership claim is introduced.
- Document assembly-status, sample-profile, official-genus and training-prevalence sensitivities. These fixed-matrix analyses do not rebuild the inventory or establish a new independently collected validation cohort.

## Presentation changes

The main text condenses the cohort discussion and moves the legacy keyword-category census and detailed co-occurrence audit to linked notes. Named candidate descriptions remain. Removing the displayed census does not alter `module_bin`, `story_role`, the housekeeping flag or their contribution to the frozen score. Figure 5C's labels change; its pattern, data, denominators and uncertainty intervals do not.

## Verification

Automated checks cover release hashes, local documentation links, inventory and tier counts, the archived category census, Figure 5 source values, unchanged scientific inputs/reference targets, full annotation replay and saved-model metrics. CI does not refit all models or validate biological mechanisms.

The model bundle records both successful original-runtime replay and a strict clean-backend probability mismatch. The latter remains a failure at the original tolerance, despite unchanged reported metrics and class decisions. Read the [model README](models/README.md) before interpreting numerical reproducibility claims.
