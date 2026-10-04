#!/usr/bin/env python3
"""Build compact manuscript evidence tables from inspected result artifacts."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

PUBLIC_RUN = ROOT / (
    "attached_extracted/Diazotrophic_Cyanobacteria_Pangenome-main (1)_extracted/"
    "Diazotrophic_Cyanobacteria_Pangenome-main/unified_pipeline_clean/"
    "unified_pipeline_run_public"
)
RECENT = ROOT / (
    "attached_extracted/DiazPengenome_extracted/DiazPengenome/"
    "Recent_Pangenome_results_ours_alone_extracted/Recent_Pangenome_results"
)
GENUS_SPLIT = ROOT / "attached_extracted/Re_ Genus CV splits (3)_extracted"
COLLAB = ROOT / "attached_extracted/collab2_vs_ours_rerun_outputs (1)_extracted"
OUT = ROOT / "results/tables"


MODULE_OVERRIDES = {
    # Manual audit 2026-06-07. These high-ranked enzymes were carried as
    # "Other / unassigned" in inherited keyword bins, but their product names
    # support the existing central-carbon module used elsewhere in the tables.
    "GF_02141": "Carbon metabolism & NAD(P)H supply",  # transketolase
    "GF_01834": "Carbon metabolism & NAD(P)H supply",  # type I GAPDH
    "GF_00856": "Carbon metabolism & NAD(P)H supply",  # ArsJ-associated GAPDH
}


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def count_rows(path: Path, delimiter: str = ",") -> int:
    return len(read_csv(path, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fieldnames} for row in rows)


def split_gene_families(value: str) -> set[str]:
    families: set[str] = set()
    for part in str(value or "").replace(",", ";").split(";"):
        token = part.strip().split()
        if token and token[0].startswith("GF_"):
            families.add(token[0])
    return families


def corrected_module(row: dict[str, str]) -> str:
    return MODULE_OVERRIDES.get(row.get("gene_family", ""), row.get("module_bin", ""))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    filter_summary_path = RECENT / "results/assemblies/filter_summary.json"
    classification_path = PUBLIC_RUN / "classification_summary.csv"
    matrix_path = PUBLIC_RUN / "gene_family_matrix.csv"
    purity_path = RECENT / "results/modeling/gene_family_purity_stats.csv"
    tier1_path = RECENT / "fox_report/tier1_positive_model_selected.tsv"
    tier2_path = RECENT / "fox_report/tier2_pure_positive_heldout.tsv"
    top100_path = RECENT / "results_tables/Table_Tier1_Top100.tsv"
    tier1_modules_path = RECENT / "results_tables/tier1_module_bin_counts.tsv"
    tier2_modules_path = RECENT / "results_tables/tier2_module_bin_counts.tsv"
    tier1_fox_path = RECENT / "results_tables/tier1_fox_status_counts.tsv"
    tier2_fox_path = RECENT / "results_tables/tier2_fox_status_counts.tsv"
    fox_map_path = ROOT / "enrichment_results/fox_to_tier_mapping.csv"
    decile_path = ROOT / "enrichment_results/decile_enrichment_results.csv"
    mannwhitney_path = ROOT / "enrichment_results/mannwhitney_results.csv"
    confounder_path = ROOT / "enrichment_results/confounder_analysis_results.csv"
    diazo_share_path = (
        GENUS_SPLIT
        / "top100_protein_vs_gene_comparison_package_extracted"
        / "top100_protein_vs_gene_diazotrophy_share.csv"
    )
    functional_comp_path = (
        GENUS_SPLIT
        / "top100_protein_vs_gene_comparison_package_extracted"
        / "top100_protein_vs_gene_functional_composition.csv"
    )
    gene_to_protein_agreement_path = GENUS_SPLIT / "gene_to_protein_top100_agreement.csv"
    cluster_three_way_path = (
        GENUS_SPLIT
        / "cluster_conservation_outputs_extracted"
        / "three_way_conserved_clusters_overlap_ge2.tsv"
    )
    cluster_pairwise_path = (
        GENUS_SPLIT
        / "cluster_conservation_outputs_extracted"
        / "pairwise_cluster_conservation_overlap_ge2.tsv"
    )
    shared_families_path = (
        GENUS_SPLIT
        / "cluster_conservation_outputs_extracted"
        / "gene_families_shared_across_all3_cluster_locations.tsv"
    )
    collab_summary_path = COLLAB / "summary_metrics.json"

    with filter_summary_path.open(encoding="utf-8") as handle:
        filter_summary = json.load(handle)

    class_rows = read_csv(classification_path)
    matrix_rows = read_csv(matrix_path)
    matrix_genomes = len(matrix_rows)
    matrix_families = max(0, len(matrix_rows[0]) - 1) if matrix_rows else 0

    purity_rows = read_csv(purity_path)
    pure_count = sum(1 for row in purity_rows if row.get("is_pure", "").lower() == "true")
    informative_count = len(purity_rows) - pure_count

    tier1_count = count_rows(tier1_path, delimiter="\t")
    tier2_count = count_rows(tier2_path, delimiter="\t")

    key_rows: list[dict[str, object]] = [
        {
            "claim": "Assemblies screened before filtering",
            "value": filter_summary["total"],
            "source": rel(filter_summary_path),
            "note": "Filter summary from inherited recent pangenome run.",
        },
        {
            "claim": "Complete RefSeq genomes retained",
            "value": filter_summary["kept"],
            "source": rel(filter_summary_path),
            "note": "Requires RefSeq GCF and complete-genome status.",
        },
        {
            "claim": "Assemblies dropped by RefSeq filter",
            "value": filter_summary["dropped_by_refseq"],
            "source": rel(filter_summary_path),
            "note": "No assemblies were dropped by completeness in this summary.",
        },
        {
            "claim": "Genomes in gene-family matrix",
            "value": matrix_genomes,
            "source": rel(matrix_path),
            "note": "Rows in presence/absence matrix.",
        },
        {
            "claim": "Gene families in matrix",
            "value": matrix_families,
            "source": rel(matrix_path),
            "note": "Columns excluding genome/accession column.",
        },
        {
            "claim": "Highly pure families removed from model feature space",
            "value": pure_count,
            "source": rel(purity_path),
            "note": "Purity threshold recorded by pipeline as >=0.90 or <=0.10.",
        },
        {
            "claim": "Informative families retained after purity filtering",
            "value": informative_count,
            "source": rel(purity_path),
            "note": "Before minimum-genus filtering in Step 4.",
        },
        {
            "claim": "Tier 1 model-supported families",
            "value": tier1_count,
            "source": rel(tier1_path),
            "note": "Positive, model-selected families in recent result bundle.",
        },
        {
            "claim": "Tier 2 high-purity held-out families",
            "value": tier2_count,
            "source": rel(tier2_path),
            "note": "High-purity positives held out from model feature ranking.",
        },
    ]

    for row in class_rows:
        key_rows.append(
            {
                "claim": f"{row['Model']} cross-validated ROC AUC",
                "value": row["roc_auc"],
                "source": rel(classification_path),
                "note": "Mean +/- SD across genus-blocked folds.",
            }
        )

    write_tsv(
        OUT / "manuscript_key_results.tsv",
        key_rows,
        ["claim", "value", "source", "note"],
    )

    top_rows = read_csv(top100_path, delimiter="\t")[:25]
    for i, row in enumerate(top_rows, start=1):
        row["rank"] = i
        row["module_bin"] = corrected_module(row)
    keep = [
        "rank",
        "gene_family",
        "product",
        "effect_size_mean",
        "consensus_rank_pct_mean",
        "n_models_with_importance",
        "diazotroph_pct_mean",
        "module_bin",
        "fox_status",
    ]
    write_tsv(
        OUT / "tier1_top25_for_manuscript.tsv",
        [{k: row.get(k, "") for k in keep} for row in top_rows],
        keep,
    )

    pg_rows = []
    for row in read_csv(diazo_share_path):
        pg_rows.append(
            {
                "comparison": row["metric"],
                "set_type": row["set_type"],
                "n": row["n"],
                "total": row["total"],
                "pct": row["pct"],
                "source": rel(diazo_share_path),
            }
        )
    write_tsv(
        OUT / "protein_vs_gene_top100_diazotrophy_share.tsv",
        pg_rows,
        ["comparison", "set_type", "n", "total", "pct", "source"],
    )

    func_rows = read_csv(functional_comp_path)
    for row in func_rows:
        row["source"] = rel(functional_comp_path)
    write_tsv(
        OUT / "protein_vs_gene_top100_functional_composition.tsv",
        func_rows,
        ["functional_bin", "n", "pct", "set_type", "source"],
    )

    gene_agreement_rows = read_csv(gene_to_protein_agreement_path)
    inventory_hits = [
        row
        for row in gene_agreement_rows
        if row.get("in_ours_tier", "").strip().lower() in {"tier1", "tier2"}
    ]
    mapped_any_gf = [
        row for row in gene_agreement_rows if row.get("mapped_gene_family", "").strip()
    ]

    three_way_cluster_rows_for_overlap = read_csv(cluster_three_way_path, delimiter="\t")
    pairwise_cluster_rows_for_overlap = read_csv(cluster_pairwise_path, delimiter="\t")
    three_way_cluster_families: set[str] = set()
    pairwise_cluster_families: set[str] = set()
    for row in three_way_cluster_rows_for_overlap:
        three_way_cluster_families.update(split_gene_families(row.get("shared_gene_families", "")))
    for row in pairwise_cluster_rows_for_overlap:
        pairwise_cluster_families.update(split_gene_families(row.get("shared_gene_families", "")))

    inventory_hit_families = {
        row.get("mapped_gene_family", "").strip()
        for row in inventory_hits
        if row.get("mapped_gene_family", "").strip()
    }
    inventory_summary_rows = [
        {
            "metric": "Gene-driven top-100 features",
            "value": len(gene_agreement_rows),
            "source": rel(gene_to_protein_agreement_path),
            "note": "Top 100 collaborator gene-driven random-forest features.",
        },
        {
            "metric": "Gene-driven features mapped to any protein family",
            "value": len(mapped_any_gf),
            "source": rel(gene_to_protein_agreement_path),
            "note": "Mapping used all genes in the gene-family JSON where available.",
        },
        {
            "metric": "Gene-driven features overlapping candidate inventories",
            "value": len(inventory_hits),
            "source": rel(gene_to_protein_agreement_path),
            "note": "Mapped to either original Tier 1 / Model-Supported or original Tier 2 / Highly Pure.",
        },
        {
            "metric": "Overlaps with Model-Supported set",
            "value": sum(1 for row in inventory_hits if row.get("in_ours_tier", "").lower() == "tier1"),
            "source": rel(gene_to_protein_agreement_path),
            "note": "Original Tier 1 source-file label.",
        },
        {
            "metric": "Overlaps with Highly Pure set",
            "value": sum(1 for row in inventory_hits if row.get("in_ours_tier", "").lower() == "tier2"),
            "source": rel(gene_to_protein_agreement_path),
            "note": "Original Tier 2 source-file label.",
        },
        {
            "metric": "Inventory-overlap gene hits in strict three-way conserved clusters",
            "value": len(inventory_hit_families & three_way_cluster_families),
            "source": rel(cluster_three_way_path),
            "note": "Overlap against shared families from three-way conserved cluster modules.",
        },
        {
            "metric": "Inventory-overlap gene hits in pairwise conserved clusters",
            "value": len(inventory_hit_families & pairwise_cluster_families),
            "source": rel(cluster_pairwise_path),
            "note": "Overlap against shared families from pairwise conserved cluster modules with at least two shared families.",
        },
    ]
    write_tsv(
        OUT / "gene_to_protein_top100_inventory_agreement_summary.tsv",
        inventory_summary_rows,
        ["metric", "value", "source", "note"],
    )

    inventory_hit_rows = []
    for row in sorted(inventory_hits, key=lambda r: int(r.get("collab_rank") or 9999)):
        source_label = row.get("in_ours_tier", "").lower()
        inventory_hit_rows.append(
            {
                "collab_rank": row.get("collab_rank", ""),
                "collab_feature": row.get("collab_feature", ""),
                "collab_rep_gene": row.get("collab_rep_gene", ""),
                "mapped_gene_family": row.get("mapped_gene_family", ""),
                "candidate_inventory": (
                    "Model-Supported" if source_label == "tier1" else "Highly Pure"
                ),
                "our_product": row.get("our_product", ""),
                "n_members_in_json": row.get("n_members_in_json", ""),
                "n_members_mapped_to_GF": row.get("n_members_mapped_to_GF", ""),
                "source": rel(gene_to_protein_agreement_path),
            }
        )
    write_tsv(
        OUT / "gene_to_protein_top100_inventory_hits.tsv",
        inventory_hit_rows,
        [
            "collab_rank",
            "collab_feature",
            "collab_rep_gene",
            "mapped_gene_family",
            "candidate_inventory",
            "our_product",
            "n_members_in_json",
            "n_members_mapped_to_GF",
            "source",
        ],
    )

    cluster_rows = read_csv(cluster_three_way_path, delimiter="\t")
    for row in cluster_rows:
        row["source"] = rel(cluster_three_way_path)
    write_tsv(
        OUT / "three_way_conserved_cluster_summary.tsv",
        cluster_rows,
        [
            "cy_cluster",
            "tri_cluster",
            "ana_cluster",
            "shared_count",
            "shared_gene_families",
            "shared_family_products",
            "source",
        ],
    )

    shared_three = read_csv(shared_families_path, delimiter="\t")
    write_tsv(
        OUT / "three_way_shared_gene_families.tsv",
        shared_three,
        [
            "gene_family",
            "product_label",
            "Cyanothece_51142_clusters",
            "Trichodesmium_IMS101_clusters",
            "Anabaena_PCC7120_clusters",
        ],
    )

    with collab_summary_path.open(encoding="utf-8") as handle:
        collab = json.load(handle)
    mapping = collab["mapping_summary"]
    collab_rows = [
        {
            "metric": key,
            "value": value,
            "source": rel(collab_summary_path),
        }
        for key, value in mapping.items()
    ]
    write_tsv(
        OUT / "collaborator_mapping_summary.tsv",
        collab_rows,
        ["metric", "value", "source"],
    )

    tier_module_rows = []
    for tier, path in [
        ("Tier 1", RECENT / "results_tables/tier1_ranked_annotated.tsv"),
        ("Tier 2", RECENT / "results_tables/tier2_pure_positive_annotated.tsv"),
    ]:
        annotated_rows = read_csv(path, delimiter="\t")
        counts = Counter(corrected_module(row) for row in annotated_rows)
        total = sum(counts.values())
        preferred_order = [
            "Other / unassigned",
            "Housekeeping / lineage marker",
            "Stress / membrane / uncharacterized",
            "Transport",
            "Regulation",
            "Respiration & bioenergetics",
            "O2 protection & redox homeostasis",
            "Metallocluster & Fe-S biogenesis",
            "Carbon metabolism & NAD(P)H supply",
            "Nitrogenase / nif machinery",
        ]
        for module in preferred_order:
            if module not in counts:
                continue
            percent = ""
            if tier == "Tier 2" and total:
                percent = f"{counts[module] / total * 100:.1f}"
            tier_module_rows.append(
                {
                    "tier": tier,
                    "module_bin": module,
                    "n_families": counts[module],
                    "percent": percent,
                    "source": rel(path),
                }
            )
    write_tsv(
        OUT / "tier_module_counts_for_manuscript.tsv",
        tier_module_rows,
        ["tier", "module_bin", "n_families", "percent", "source"],
    )

    tier_fox_rows = []
    for tier, path in [("Tier 1", tier1_fox_path), ("Tier 2", tier2_fox_path)]:
        for row in read_csv(path, delimiter="\t"):
            tier_fox_rows.append(
                {
                    "tier": tier,
                    "fox_status": row.get("fox_status", ""),
                    "n_families": row.get("n_families", ""),
                    "percent": row.get("percent", ""),
                    "source": rel(path),
                }
            )
    write_tsv(
        OUT / "tier_fox_status_for_manuscript.tsv",
        tier_fox_rows,
        ["tier", "fox_status", "n_families", "percent", "source"],
    )

    fox_rows = read_csv(fox_map_path)
    fox_summary = [
        {
            "metric": "FOX unknown/proteomics rows",
            "value": len(fox_rows),
            "source": rel(fox_map_path),
            "note": "Rows in FOX_unknown_app_table mapping output.",
        },
        {
            "metric": "Rows captured in tier by direct WP mapping",
            "value": sum(1 for row in fox_rows if row.get("is_tier_wp", "").lower() == "true"),
            "source": rel(fox_map_path),
            "note": "Direct WP-to-family tier mapping.",
        },
        {
            "metric": "Tier 1 rows by direct WP mapping",
            "value": sum(1 for row in fox_rows if row.get("is_tier1_wp", "").lower() == "true"),
            "source": rel(fox_map_path),
            "note": "Direct WP-to-family tier mapping.",
        },
        {
            "metric": "Tier 2 rows by direct WP mapping",
            "value": sum(1 for row in fox_rows if row.get("is_tier2_wp", "").lower() == "true"),
            "source": rel(fox_map_path),
            "note": "Direct WP-to-family tier mapping.",
        },
        {
            "metric": "Rows captured in tier by combined WP/product mapping",
            "value": sum(1 for row in fox_rows if row.get("is_tier_combined", "").lower() == "true"),
            "source": rel(fox_map_path),
            "note": "Combined direct WP and product-name mapping.",
        },
    ]

    for row in read_csv(decile_path):
        if row.get("mapping_strategy") == "WP_only" and row.get("probability_column") == "prob_ensemble_with_position":
            fox_summary.extend(
                [
                    {
                        "metric": "WP-only FOX ensemble top-decile tier fraction",
                        "value": row.get("top_decile_frac", ""),
                        "source": rel(decile_path),
                        "note": "Top FOX-score decile is not enriched for tier membership.",
                    },
                    {
                        "metric": "WP-only FOX ensemble bottom-decile tier fraction",
                        "value": row.get("bottom_decile_frac", ""),
                        "source": rel(decile_path),
                        "note": "Bottom FOX-score decile has higher tier fraction in this output.",
                    },
                    {
                        "metric": "WP-only FOX ensemble enrichment ratio",
                        "value": row.get("enrichment_ratio", ""),
                        "source": rel(decile_path),
                        "note": "Ratio below 1 argues against simple positive FOX-score enrichment.",
                    },
                ]
            )
            break

    for row in read_csv(mannwhitney_path):
        if row.get("mapping_strategy") == "WP_only" and row.get("probability_column") == "prob_ensemble_with_position":
            fox_summary.extend(
                [
                    {
                        "metric": "WP-only tier FOX ensemble median",
                        "value": row.get("tier_median", ""),
                        "source": rel(mannwhitney_path),
                        "note": "Tier rows have lower median score than non-tier rows in this test.",
                    },
                    {
                        "metric": "WP-only non-tier FOX ensemble median",
                        "value": row.get("nontier_median", ""),
                        "source": rel(mannwhitney_path),
                        "note": "Comparison median from Mann-Whitney output.",
                    },
                    {
                        "metric": "WP-only Mann-Whitney p",
                        "value": row.get("p_twosided", ""),
                        "source": rel(mannwhitney_path),
                        "note": "Significant in the opposite direction for simple enrichment.",
                    },
                ]
            )
            break

    for row in read_csv(confounder_path):
        if row.get("confounder") == "diazo_hit_frac" and row.get("stratum") in {
            "diazo_hit_frac_leq_0.5",
            "diazo_hit_frac_leq_0.25",
        }:
            fox_summary.append(
                {
                    "metric": f"FOX score tier-vs-nontier p after {row.get('stratum')}",
                    "value": row.get("p_twosided", ""),
                    "source": rel(confounder_path),
                    "note": "Difference is not significant after diazotrophy-hit-fraction restriction.",
                }
            )

    write_tsv(
        OUT / "fox_proteomics_overlay_summary.tsv",
        fox_summary,
        ["metric", "value", "source", "note"],
    )

    print(f"Wrote manuscript evidence tables to {OUT}")


if __name__ == "__main__":
    main()
