#!/usr/bin/env python3
"""Map condensate-driver ranking output onto the diazotrophy family inventory.

The input ranking is a cyano-only protein-level list with a driver score and
rank. This script maps WP accessions in either ``uniprot_id`` or ``cluster_id``
onto the current pangenome ``GF_`` namespace, aggregates to one row per family,
and tests whether candidate sets are enriched among high condensate-driver
families.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
RANKING = ROOT / "data" / "full_ranking_conditionC_cyano_seed42_rescued_direct_plus_homologs.csv"
PANGENOME_MAP = ROOT / "_pangenome_data" / "genome_protein_family_map.tsv"
INVENTORY = ROOT / "evidence_tables" / "morphotype_bridge_family_inventory.tsv"
OUT_DIR = ROOT / "evidence_tables"


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def label_join(values: pd.Series, limit: int = 12) -> str:
    vals = sorted({str(v) for v in values.dropna() if str(v) and str(v) != "nan"})
    return "; ".join(vals[:limit])


def main() -> None:
    ranking = pd.read_csv(RANKING)
    ranking["in_training_positive"] = as_bool(ranking["in_training_positive"])
    ranking["withheld_family"] = as_bool(ranking["withheld_family"])

    pg = pd.read_csv(PANGENOME_MAP, sep="\t", dtype=str)
    protein_to_family = pg.drop_duplicates("protein_accession").set_index("protein_accession")[
        "gene_family"
    ]

    ranking["family_from_uniprot_id"] = ranking["uniprot_id"].map(protein_to_family)
    ranking["family_from_cluster_id"] = ranking["cluster_id"].map(protein_to_family)
    ranking["mapped_gene_family"] = ranking["family_from_uniprot_id"].fillna(
        ranking["family_from_cluster_id"]
    )
    ranking["mapping_method"] = "unmapped"
    ranking.loc[ranking["family_from_cluster_id"].notna(), "mapping_method"] = "cluster_id_to_GF"
    ranking.loc[ranking["family_from_uniprot_id"].notna(), "mapping_method"] = "uniprot_id_to_GF"
    ranking["mapping_conflict"] = (
        ranking["family_from_uniprot_id"].notna()
        & ranking["family_from_cluster_id"].notna()
        & (ranking["family_from_uniprot_id"] != ranking["family_from_cluster_id"])
    )

    mapped = ranking[ranking["mapped_gene_family"].notna()].copy()
    mapped["driver_rank_cyano"] = pd.to_numeric(mapped["driver_rank_cyano"], errors="coerce")
    mapped["driver_score"] = pd.to_numeric(mapped["driver_score"], errors="coerce")

    family = (
        mapped.groupby("mapped_gene_family")
        .agg(
            condensate_mapped_rows=("uniprot_id", "size"),
            condensate_best_score=("driver_score", "max"),
            condensate_mean_score=("driver_score", "mean"),
            condensate_median_score=("driver_score", "median"),
            condensate_best_rank_cyano=("driver_rank_cyano", "min"),
            condensate_training_positive_rows=("in_training_positive", "sum"),
            condensate_withheld_rows=("withheld_family", "sum"),
            condensate_label_names=("label_name", label_join),
            condensate_label_sources=("label_source", label_join),
            condensate_best_products=("protein_name", label_join),
            condensate_best_organisms=("organism", label_join),
            condensate_mapping_methods=("mapping_method", label_join),
            condensate_mapping_conflict_rows=("mapping_conflict", "sum"),
        )
        .reset_index()
        .rename(columns={"mapped_gene_family": "gene_family"})
    )

    family = family.sort_values(["condensate_best_score", "condensate_best_rank_cyano"], ascending=[False, True])
    family["condensate_family_rank"] = range(1, len(family) + 1)
    family["condensate_family_percentile"] = family["condensate_family_rank"] / len(family)

    inv = pd.read_csv(INVENTORY, sep="\t", dtype=str)
    family = family.merge(
        inv[
            [
                "gene_family",
                "candidate_set",
                "product",
                "module_bin",
                "fox_status",
                "consensus_rank_pct_mean",
                "diazotroph_pct_mean",
                "primary_bridge_broad_ge1u_ge1f",
                "strict_unicellular_breadth_ge3u_ge1f",
            ]
        ],
        on="gene_family",
        how="left",
    )
    family["candidate_set"] = family["candidate_set"].fillna("Not in candidate inventory")
    family["is_model_supported"] = family["candidate_set"].eq("Model-Supported")
    family["is_highly_pure"] = family["candidate_set"].eq("Highly Pure")
    family["is_candidate_inventory"] = family["is_model_supported"] | family["is_highly_pure"]

    thresholds = [
        ("top_1pct", 0.01),
        ("top_5pct", 0.05),
        ("top_10pct", 0.10),
        ("top_20pct", 0.20),
    ]
    for name, frac in thresholds:
        family[name] = family["condensate_family_percentile"].le(frac)

    rows: list[dict[str, object]] = []
    groups = [
        ("Model-Supported", family["is_model_supported"]),
        ("Highly Pure", family["is_highly_pure"]),
        ("Any candidate inventory", family["is_candidate_inventory"]),
    ]
    for group_name, group_mask in groups:
        scores_group = family.loc[group_mask, "condensate_best_score"]
        scores_other = family.loc[~group_mask, "condensate_best_score"]
        if len(scores_group) and len(scores_other):
            mw = mannwhitneyu(scores_group, scores_other, alternative="greater")
            rows.append(
                {
                    "test": "mannwhitney_best_score_greater",
                    "group": group_name,
                    "threshold": "",
                    "group_n": int(group_mask.sum()),
                    "other_n": int((~group_mask).sum()),
                    "group_hits": "",
                    "other_hits": "",
                    "odds_ratio": "",
                    "p_value": mw.pvalue,
                    "group_median_score": scores_group.median(),
                    "other_median_score": scores_other.median(),
                    "group_mean_score": scores_group.mean(),
                    "other_mean_score": scores_other.mean(),
                }
            )
        for thresh_name, _ in thresholds:
            hit = family[thresh_name]
            table = [
                [int((group_mask & hit).sum()), int((group_mask & ~hit).sum())],
                [int((~group_mask & hit).sum()), int((~group_mask & ~hit).sum())],
            ]
            or_val, p_val = fisher_exact(table, alternative="greater")
            rows.append(
                {
                    "test": "fisher_top_family_percentile_greater",
                    "group": group_name,
                    "threshold": thresh_name,
                    "group_n": int(group_mask.sum()),
                    "other_n": int((~group_mask).sum()),
                    "group_hits": table[0][0],
                    "other_hits": table[1][0],
                    "odds_ratio": or_val,
                    "p_value": p_val,
                    "group_median_score": family.loc[group_mask, "condensate_best_score"].median(),
                    "other_median_score": family.loc[~group_mask, "condensate_best_score"].median(),
                    "group_mean_score": family.loc[group_mask, "condensate_best_score"].mean(),
                    "other_mean_score": family.loc[~group_mask, "condensate_best_score"].mean(),
                }
            )

    summary = pd.DataFrame(rows)

    mapping_summary = pd.DataFrame(
        [
            {"metric": "input_rows", "value": len(ranking)},
            {"metric": "mapped_rows", "value": len(mapped)},
            {"metric": "mapped_unique_families", "value": family["gene_family"].nunique()},
            {"metric": "mapping_conflict_rows", "value": int(ranking["mapping_conflict"].sum())},
            {
                "metric": "uniprot_id_direct_mapped_rows",
                "value": int(ranking["family_from_uniprot_id"].notna().sum()),
            },
            {
                "metric": "cluster_id_mapped_rows",
                "value": int(ranking["family_from_cluster_id"].notna().sum()),
            },
            {
                "metric": "model_supported_mapped_families",
                "value": int(family["is_model_supported"].sum()),
            },
            {"metric": "highly_pure_mapped_families", "value": int(family["is_highly_pure"].sum())},
        ]
    )

    family.to_csv(OUT_DIR / "condensate_ranking_family_overlay.tsv", sep="\t", index=False)
    summary.to_csv(OUT_DIR / "condensate_ranking_enrichment_summary.tsv", sep="\t", index=False)
    mapping_summary.to_csv(OUT_DIR / "condensate_ranking_mapping_summary.tsv", sep="\t", index=False)
    family[family["is_model_supported"]].sort_values("condensate_family_rank").head(100).to_csv(
        OUT_DIR / "condensate_ranking_top_model_supported.tsv", sep="\t", index=False
    )

    print(mapping_summary.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
