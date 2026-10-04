Portable saved-output rank audit (no model fitting)

Run from this folder with Python3.10.16, numpy2.2.6, pandas2.3.3,
SciPy loaded-module version1.13.1:
    python run_rank_audit.py --out fresh_results --verify

Runtime provenance caveat: this existing Anaconda environment reports
scipy.__version__=1.13.1 but importlib.metadata.version('scipy')=1.13.0.
Both observations are retained; it is not a verified clean installation of
either distribution. Exact saved-table agreement was checked in this existing
runtime. A fresh installation/cross-platform equivalence was not tested.

This recomputes rank_comparison_summary.csv, rank_comparison_family_detail.csv
and all_frozen476_repeated_LR_stability.csv from the supplied model outputs.
The input manifest records staged-file and upstream-source SHA256 hashes.
Coefficient inputs retain seed/fold/family/coefficient for the50 exact-boundary
training-only Logistic Regression fits (10 dependent partitions, seeds42–51).
Global replay inputs retain model/fold/family/importance, without other arms.
Original and corrected seed42 importance exports and frozen476 inventory are
copied byte-for-byte. This bundle does not reproduce model fitting; use the
separate portable model audit for that purpose.

The frozen476 resource is preserved:473 estimable positive univariate
carrier/noncarrier contrasts plus3 ubiquitous empty-comparison edge cases.
All476 have rows in all four importance exports; zero values were admitted.
Mean within-model percentile ranking orders the resource, without imposing
an additional importance threshold. Corrections never silently replace it.

Rank correlations and topN overlap are descriptive. TopN ties are broken by
gene-family identifier. In repeated LR fits, unselected features contribute
zero to the mean absolute coefficient; signed-coefficient stability is counted
only when selected. Signs are conditional on correlated predictors and do not
establish biological direction or family-level validation.
