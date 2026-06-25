Composite proteomics-prioritized family evidence table

This folder joins the current 426-genome protein-family candidate inventory to:
- accession-mapped external literature proteomics studies filtered to
  nitrogen-fixing-condition response evidence;
- related/collaborator protein-family atlas concordance;
- condensate-driver ranking;
- morphotype-breadth flags.

The score is a prioritization/down-filter, not a formal validation statistic.
External literature proteomics contributes to the score only when the mapped
source row increases in the organism-appropriate nitrogen-fixation-active
phase/cell. Broad N-fix-condition response, active-phase-down, transition, and
non-target phase rows are retained as metadata but are not rewarded. Local
Cyanothece proteomics and prior FOX probability layers are excluded from both
scoring and visible composite output. Core nif,
lineage/housekeeping, and HGT-passenger-like features are retained but penalized
in accessory_story_score so they do not dominate the accessory-biology narrative.
