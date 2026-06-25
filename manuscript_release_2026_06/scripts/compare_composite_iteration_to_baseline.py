from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COMPOSITE = ROOT / "evidence_tables" / "composite_proteomics_priority"
BASELINE = COMPOSITE / "baseline_2026-06-23" / "family_composite_all.csv"
CURRENT = COMPOSITE / "family_composite_all.csv"


def main() -> None:
    base = pd.read_csv(BASELINE)
    current = pd.read_csv(CURRENT)

    keep = [
        "gene_family",
        "candidate_set",
        "product",
        "story_role",
        "priority_tier",
        "composite_evidence_score",
        "accessory_story_score",
        "scope_adjusted_story_score",
        "n_literature_studies",
        "n_literature_nfix_response_studies",
        "evidence_alignment_summary",
    ]
    base = base[[c for c in keep if c in base.columns]].add_prefix("baseline_")
    current = current[[c for c in keep if c in current.columns]].add_prefix("current_")

    delta = current.merge(
        base,
        left_on="current_gene_family",
        right_on="baseline_gene_family",
        how="outer",
    )
    delta["gene_family"] = delta["current_gene_family"].fillna(delta["baseline_gene_family"])
    delta["tier_changed"] = delta["current_priority_tier"] != delta["baseline_priority_tier"]
    delta["score_delta"] = (
        pd.to_numeric(delta["current_accessory_story_score"], errors="coerce")
        - pd.to_numeric(delta["baseline_accessory_story_score"], errors="coerce")
    )
    current_scope = (
        pd.to_numeric(delta["current_scope_adjusted_story_score"], errors="coerce")
        if "current_scope_adjusted_story_score" in delta.columns
        else pd.to_numeric(delta["current_accessory_story_score"], errors="coerce")
    )
    baseline_scope = (
        pd.to_numeric(delta["baseline_scope_adjusted_story_score"], errors="coerce")
        if "baseline_scope_adjusted_story_score" in delta.columns
        else pd.to_numeric(delta["baseline_accessory_story_score"], errors="coerce")
    )
    delta["scope_adjusted_score_delta"] = current_scope - baseline_scope
    delta["score_changed"] = delta["score_delta"].abs().fillna(0) > 1e-9
    delta["change_class"] = "unchanged"
    delta.loc[delta["baseline_gene_family"].isna(), "change_class"] = "added"
    delta.loc[delta["current_gene_family"].isna(), "change_class"] = "removed"
    delta.loc[
        (delta["change_class"] == "unchanged") & delta["tier_changed"],
        "change_class",
    ] = "tier_changed"
    delta.loc[
        (delta["change_class"] == "unchanged") & delta["score_changed"],
        "change_class",
    ] = "score_changed_only"

    ordered = [
        "gene_family",
        "change_class",
        "score_delta",
        "current_priority_tier",
        "baseline_priority_tier",
        "current_accessory_story_score",
        "baseline_accessory_story_score",
        "current_scope_adjusted_story_score",
        "baseline_scope_adjusted_story_score",
        "current_candidate_set",
        "current_product",
        "current_story_role",
        "baseline_story_role",
        "current_evidence_alignment_summary",
        "baseline_evidence_alignment_summary",
    ]
    delta = delta[[c for c in ordered if c in delta.columns] + [c for c in delta.columns if c not in ordered]]
    delta = delta.sort_values(
        ["change_class", "score_delta", "gene_family"],
        ascending=[True, False, True],
        na_position="last",
    )

    summary = (
        delta["change_class"]
        .value_counts(dropna=False)
        .rename_axis("change_class")
        .reset_index(name="n_families")
    )
    tier_summary = (
        delta[delta["change_class"].isin(["tier_changed", "score_changed_only"])]
        .groupby(["baseline_priority_tier", "current_priority_tier"], dropna=False)
        .size()
        .reset_index(name="n_families")
    )

    delta.to_csv(COMPOSITE / "iteration_delta_vs_baseline.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    summary.to_csv(COMPOSITE / "iteration_delta_summary.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    tier_summary.to_csv(COMPOSITE / "iteration_tier_transition_summary.csv", index=False, quoting=csv.QUOTE_MINIMAL)

    print(summary.to_string(index=False))
    if not tier_summary.empty:
        print()
        print(tier_summary.to_string(index=False))


if __name__ == "__main__":
    main()
