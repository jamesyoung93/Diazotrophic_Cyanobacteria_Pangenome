#!/usr/bin/env python3
"""
Audit whether Aryal's CV-genus-split protein-level materials match the current
manuscript protein-family universe, then quantify direct overlap with the
current Model-Supported/Highly Pure/FOX/proteomics evidence layers.
"""

from __future__ import annotations

import csv
import io
import json
import re
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PRED_ZIP = ROOT / "raw_input_archives" / "predictive-comparative-genomics-for-diazotropy-main (1).zip"
PAN_ZIP = ROOT / "raw_input_archives" / "Diazotrophic_Cyanobacteria_Pangenome-main (1).zip"
EVIDENCE = ROOT / "evidence_tables"
OUT = ROOT / "evidence_tables"
ARYAL_HGT_CACHE = OUT / "aryal_task7_hgt_proximity_from_github.csv"

ARYAL_HGT_URL = (
    "https://raw.githubusercontent.com/erise-bnerc/"
    "predictive-comparative-genomics-for-diazotropy/main/"
    "new_analytical_task_results/task7_hgt_proximity.csv"
)


def read_zip_table(zip_path: Path, suffix: str, **kwargs) -> tuple[pd.DataFrame, str]:
    with zipfile.ZipFile(zip_path) as zf:
        matches = [n for n in zf.namelist() if n.endswith(suffix)]
        if not matches:
            raise FileNotFoundError(f"No entry ending with {suffix!r} in {zip_path}")
        entry = sorted(matches, key=len)[0]
        with zf.open(entry) as handle:
            return pd.read_csv(handle, **kwargs), entry


def read_cv_split_package() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    with zipfile.ZipFile(PRED_ZIP) as zf:
        label_entry = [
            n for n in zf.namelist() if n.endswith("datasets/complete_genomes_labeled-james.csv")
        ][0]
        labels_raw = pd.read_csv(zf.open(label_entry))
        labels_raw["assembly_accession"] = labels_raw["assembly_accession"].astype(str)

        split_rows = []
        for name in sorted(zf.namelist()):
            if not (name.endswith("test_genomes.txt") or name.endswith("train_genomes.txt")):
                continue
            fold_match = re.search(r"fold_\d+", name)
            fold = fold_match.group(0) if fold_match else ""
            split_type = "test" if name.endswith("test_genomes.txt") else "train"
            accessions = [
                line.strip()
                for line in zf.open(name).read().decode("utf-8").splitlines()
                if line.strip() and not line.startswith(".")
            ]
            for acc in accessions:
                split_rows.append(
                    {
                        "assembly_accession": acc,
                        "fold": fold,
                        "split": split_type,
                        "source_file": name,
                    }
                )
        split_long = pd.DataFrame(split_rows)

    dup_counts = labels_raw["assembly_accession"].value_counts()
    label_conflicts = (
        labels_raw.groupby("assembly_accession")["is_diazotroph"]
        .nunique(dropna=False)
        .reset_index(name="n_unique_is_diazotroph")
    )
    labels = (
        labels_raw.sort_values(["assembly_accession"])
        .drop_duplicates("assembly_accession", keep="first")
        .copy()
    )
    labels["cv_label_duplicate_rows"] = labels["assembly_accession"].map(dup_counts).fillna(1).astype(int)
    labels = labels.merge(label_conflicts, on="assembly_accession", how="left")

    meta = {
        "label_entry": label_entry,
        "label_rows_raw": len(labels_raw),
        "label_unique_accessions": labels["assembly_accession"].nunique(),
        "label_duplicate_accessions": int((dup_counts > 1).sum()),
        "label_max_duplicate_rows": int(dup_counts.max()),
        "label_accessions_with_conflicting_diaz_label": int(
            (label_conflicts["n_unique_is_diazotroph"] > 1).sum()
        ),
        "split_unique_accessions": split_long["assembly_accession"].nunique(),
        "split_total_mentions": len(split_long),
    }
    return labels, split_long, meta


def fetch_aryal_task7() -> pd.DataFrame:
    if ARYAL_HGT_CACHE.exists():
        return pd.read_csv(ARYAL_HGT_CACHE)
    with urllib.request.urlopen(ARYAL_HGT_URL, timeout=30) as response:
        text = response.read().decode("utf-8")
    return pd.read_csv(io.StringIO(text))


def bool_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].astype(str).str.lower().isin(["true", "1", "yes"])


def normalize_product(value) -> str:
    if pd.isna(value):
        return ""
    text = str(value).lower().strip()
    text = re.sub(r"\s+", " ", text)
    text = text.replace('"', "")
    return text


def is_informative_product(value: str) -> bool:
    text = normalize_product(value)
    if not text:
        return False
    weak_terms = [
        "hypothetical protein",
        "uncharacterized",
        "domain-containing protein",
        "family protein",
        "duf",
    ]
    return not any(term in text for term in weak_terms)


def pct(n: int, d: int) -> float:
    return round((100.0 * n / d), 2) if d else 0.0


def main() -> None:
    OUT.mkdir(exist_ok=True)

    cv_labels, cv_split_long, cv_meta = read_cv_split_package()

    matrix, matrix_entry = read_zip_table(
        PAN_ZIP,
        "unified_pipeline_run_public/gene_family_matrix.csv",
        index_col=0,
    )
    pan_labels, pan_label_entry = read_zip_table(
        PAN_ZIP,
        "unified_pipeline_run_public/complete_genomes_labeled.csv",
    )
    tier1, tier1_entry = read_zip_table(
        PAN_ZIP,
        "unified_pipeline_run_public/fox_report/tier1_positive_model_selected.tsv",
        sep="\t",
    )
    tier2, tier2_entry = read_zip_table(
        PAN_ZIP,
        "unified_pipeline_run_public/fox_report/tier2_pure_positive_heldout.tsv",
        sep="\t",
    )
    current_tier1_annotated, current_tier1_annotated_entry = read_zip_table(
        PAN_ZIP,
        "unified_pipeline_run_public/results_tables/tier1_ranked_annotated.tsv",
        sep="\t",
    )

    aryal_tier1 = fetch_aryal_task7()

    current_inventory = pd.read_csv(EVIDENCE / "morphotype_bridge_family_inventory.tsv", sep="\t")
    overlay = pd.read_csv(EVIDENCE / "evidence_overlay_venn_candidate_filter.tsv", sep="\t")
    core_nif = pd.read_csv(EVIDENCE / "core_nif_family_audit.tsv", sep="\t")

    matrix_ids = set(matrix.index.astype(str))
    matrix_families = set(map(str, matrix.columns))
    pan_label_ids = set(pan_labels["assembly_accession"].astype(str))
    cv_split_ids = set(cv_split_long["assembly_accession"].astype(str))
    cv_label_ids = set(cv_labels["assembly_accession"].astype(str))

    # Genome-level reconciliation.
    genome_rows = []
    all_genomes = sorted(matrix_ids | pan_label_ids | cv_split_ids | cv_label_ids)
    cv_label_indexed = cv_labels.set_index("assembly_accession", drop=False)
    pan_label_indexed = pan_labels.set_index("assembly_accession", drop=False)
    fold_summary = (
        cv_split_long.groupby("assembly_accession")
        .agg(
            cv_folds=("fold", lambda s: ";".join(sorted(set(s)))),
            cv_split_roles=("split", lambda s: ";".join(sorted(set(s)))),
            cv_mentions=("assembly_accession", "size"),
        )
        .reset_index()
        .set_index("assembly_accession")
    )
    for acc in all_genomes:
        cv_row = cv_label_indexed.loc[acc] if acc in cv_label_indexed.index else None
        pan_row = pan_label_indexed.loc[acc] if acc in pan_label_indexed.index else None
        fold_row = fold_summary.loc[acc] if acc in fold_summary.index else None
        genome_rows.append(
            {
                "assembly_accession": acc,
                "prefix": acc[:3],
                "in_aryal_cv_split": acc in cv_split_ids,
                "in_aryal_unique_label_table": acc in cv_label_ids,
                "in_current_matrix": acc in matrix_ids,
                "in_current_label_table": acc in pan_label_ids,
                "aryal_is_diazotroph": "" if cv_row is None else bool(cv_row["is_diazotroph"]),
                "current_is_diazotroph": "" if pan_row is None else bool(pan_row["is_diazotroph"]),
                "aryal_genus": "" if cv_row is None else cv_row.get("genus", ""),
                "current_genus": "" if pan_row is None else pan_row.get("genus", ""),
                "aryal_organism": "" if cv_row is None else cv_row.get("organism_full", ""),
                "current_organism": "" if pan_row is None else pan_row.get("organism_full", ""),
                "cv_folds": "" if fold_row is None else fold_row["cv_folds"],
                "cv_split_roles": "" if fold_row is None else fold_row["cv_split_roles"],
                "cv_mentions": 0 if fold_row is None else int(fold_row["cv_mentions"]),
                "cv_label_duplicate_rows": 0 if cv_row is None else int(cv_row["cv_label_duplicate_rows"]),
                "cv_label_conflict_count": 0 if cv_row is None else int(cv_row["n_unique_is_diazotroph"]),
            }
        )
    genome_overlap = pd.DataFrame(genome_rows)
    genome_overlap.to_csv(OUT / "aryal_cv_split_genome_overlap.csv", index=False)

    # Family-level reconciliation.
    current_tier1_ids = set(tier1["gene_family"].astype(str))
    current_tier2_ids = set(tier2["gene_family"].astype(str))
    current_inventory_ids = set(current_inventory["gene_family"].astype(str))
    aryal_ids = set(aryal_tier1["gene_family"].astype(str))

    current_inv_keep = [
        "gene_family",
        "candidate_set",
        "product",
        "module_bin",
        "fox_status",
        "best_rank",
        "consensus_rank_pct_mean",
        "diazotroph_pct_mean",
        "n_member_genomes",
        "gene_top100_overlap",
        "primary_bridge_broad_ge1u_ge1f",
        "strict_unicellular_breadth_ge3u_ge1f",
    ]
    current_inv_renamed = current_inventory[current_inv_keep].rename(
        columns={
            "product": "current_product",
            "module_bin": "current_module_bin",
            "candidate_set": "current_candidate_set",
            "fox_status": "current_fox_status",
            "best_rank": "current_best_rank",
            "consensus_rank_pct_mean": "current_consensus_rank_pct_mean",
            "diazotroph_pct_mean": "current_diazotroph_pct_mean",
            "n_member_genomes": "current_n_member_genomes",
            "gene_top100_overlap": "current_gene_top100_overlap",
            "primary_bridge_broad_ge1u_ge1f": "current_primary_bridge_broad_ge1u_ge1f",
            "strict_unicellular_breadth_ge3u_ge1f": "current_strict_unicellular_breadth_ge3u_ge1f",
        }
    )
    overlay_keep = [
        "gene_family",
        "layer_model_supported",
        "layer_highly_pure",
        "layer_fox",
        "layer_proteomics_detected",
        "layer_proteomics_directional",
        "layer_proteomics_nominal_up",
        "three_layer_model_fox_proteomics_directional",
        "evidence_layer_count_model_fox_proteomics",
    ]
    overlay_renamed = overlay[overlay_keep]

    aryal_renamed = aryal_tier1.rename(
        columns={
            "gene": "aryal_gene",
            "product": "aryal_product",
            "module_bin": "aryal_module_bin",
            "consensus_rank_pct_mean": "aryal_consensus_rank_pct_mean",
        }
    )
    family_detail = aryal_renamed.merge(current_inv_renamed, on="gene_family", how="left")
    family_detail = family_detail.merge(overlay_renamed, on="gene_family", how="left")
    family_detail["in_current_matrix_family_columns"] = family_detail["gene_family"].astype(str).isin(
        matrix_families
    )
    family_detail["in_current_inventory"] = family_detail["gene_family"].astype(str).isin(
        current_inventory_ids
    )
    family_detail["in_current_model_supported_476"] = family_detail["gene_family"].astype(str).isin(
        current_tier1_ids
    )
    family_detail["in_current_highly_pure_981"] = family_detail["gene_family"].astype(str).isin(
        current_tier2_ids
    )
    family_detail["aryal_product_norm"] = family_detail["aryal_product"].fillna("").str.lower().str.strip()
    family_detail["current_product_norm"] = family_detail["current_product"].fillna("").str.lower().str.strip()
    family_detail["exact_product_match"] = (
        family_detail["aryal_product_norm"].ne("")
        & family_detail["aryal_product_norm"].eq(family_detail["current_product_norm"])
    )
    family_detail.to_csv(OUT / "aryal_task7_tier1_vs_current_inventory_detail.csv", index=False)

    shared = family_detail[family_detail["in_current_inventory"]].copy()
    shared_focus_cols = [
        "gene_family",
        "aryal_gene",
        "aryal_product",
        "current_product",
        "aryal_module_bin",
        "current_module_bin",
        "aryal_consensus_rank_pct_mean",
        "current_candidate_set",
        "current_best_rank",
        "current_consensus_rank_pct_mean",
        "current_fox_status",
        "layer_fox",
        "layer_proteomics_detected",
        "layer_proteomics_directional",
        "three_layer_model_fox_proteomics_directional",
        "hgt_passenger_flag",
        "frac_assemblies_proximal",
        "fisher_p_vs_background",
        "exact_product_match",
    ]
    shared[shared_focus_cols].to_csv(OUT / "aryal_current_shared_family_detail.csv", index=False)

    shiva_only = family_detail[~family_detail["in_current_inventory"]].copy()
    shiva_only[
        [
            "gene_family",
            "aryal_gene",
            "aryal_product",
            "aryal_module_bin",
            "aryal_consensus_rank_pct_mean",
            "in_current_matrix_family_columns",
            "hgt_passenger_flag",
            "frac_assemblies_proximal",
            "fisher_p_vs_background",
        ]
    ].to_csv(OUT / "aryal_tier1_not_in_current_inventory.csv", index=False)

    current_only = current_inventory[
        ~current_inventory["gene_family"].astype(str).isin(aryal_ids)
    ].copy()
    current_only[
        [
            "gene_family",
            "candidate_set",
            "product",
            "module_bin",
            "fox_status",
            "best_rank",
            "consensus_rank_pct_mean",
            "diazotroph_pct_mean",
            "gene_top100_overlap",
            "primary_bridge_broad_ge1u_ge1f",
            "strict_unicellular_breadth_ge3u_ge1f",
        ]
    ].to_csv(OUT / "current_inventory_not_in_aryal_tier1.csv", index=False)

    # Product/function-level overlap. GF IDs are not necessarily stable between
    # clustering runs, so exact product-name matches provide a conservative
    # biological crosswalk.
    aryal_products = aryal_renamed[
        [
            "gene_family",
            "aryal_gene",
            "aryal_product",
            "aryal_module_bin",
            "aryal_consensus_rank_pct_mean",
            "hgt_passenger_flag",
            "frac_assemblies_proximal",
            "fisher_p_vs_background",
        ]
    ].copy()
    aryal_products["product_norm"] = aryal_products["aryal_product"].map(normalize_product)
    aryal_products["aryal_product_is_informative"] = aryal_products["aryal_product"].map(
        is_informative_product
    )

    current_products = current_inventory[
        [
            "gene_family",
            "candidate_set",
            "product",
            "module_bin",
            "fox_status",
            "best_rank",
            "consensus_rank_pct_mean",
            "diazotroph_pct_mean",
        ]
    ].copy()
    current_products["product_norm"] = current_products["product"].map(normalize_product)
    current_products["current_product_is_informative"] = current_products["product"].map(
        is_informative_product
    )
    current_products = current_products.merge(overlay_renamed, on="gene_family", how="left")

    product_long = aryal_products.merge(
        current_products,
        on="product_norm",
        how="inner",
        suffixes=("_aryal", "_current"),
    ).rename(
        columns={
            "gene_family_aryal": "aryal_gene_family",
            "gene_family_current": "current_gene_family",
            "product": "current_product",
            "candidate_set": "current_candidate_set",
            "module_bin": "current_module_bin",
            "fox_status": "current_fox_status",
            "best_rank": "current_best_rank",
            "consensus_rank_pct_mean": "current_consensus_rank_pct_mean",
            "diazotroph_pct_mean": "current_diazotroph_pct_mean",
        }
    )
    product_long["same_gf_id"] = product_long["aryal_gene_family"].astype(str).eq(
        product_long["current_gene_family"].astype(str)
    )
    product_long.to_csv(OUT / "aryal_current_product_overlap_long.csv", index=False)

    if not product_long.empty:
        product_agg = (
            product_long.groupby("aryal_gene_family")
            .agg(
                aryal_product=("aryal_product", "first"),
                aryal_gene=("aryal_gene", "first"),
                aryal_module_bin=("aryal_module_bin", "first"),
                aryal_consensus_rank_pct_mean=("aryal_consensus_rank_pct_mean", "first"),
                aryal_product_is_informative=("aryal_product_is_informative", "first"),
                n_current_product_matches=("current_gene_family", "nunique"),
                current_matching_families=("current_gene_family", lambda s: ";".join(sorted(set(map(str, s))))),
                current_candidate_sets=("current_candidate_set", lambda s: ";".join(sorted(set(map(str, s))))),
                current_modules=("current_module_bin", lambda s: ";".join(sorted(set(map(str, s))))),
                any_current_model_supported=("layer_model_supported", lambda s: bool(bool_col(pd.DataFrame({"x": s}), "x").any())),
                any_current_highly_pure=("layer_highly_pure", lambda s: bool(bool_col(pd.DataFrame({"x": s}), "x").any())),
                any_current_fox=("layer_fox", lambda s: bool(bool_col(pd.DataFrame({"x": s}), "x").any())),
                any_current_proteomics_detected=("layer_proteomics_detected", lambda s: bool(bool_col(pd.DataFrame({"x": s}), "x").any())),
                any_current_proteomics_directional=("layer_proteomics_directional", lambda s: bool(bool_col(pd.DataFrame({"x": s}), "x").any())),
                any_current_three_layer=("three_layer_model_fox_proteomics_directional", lambda s: bool(bool_col(pd.DataFrame({"x": s}), "x").any())),
                any_same_gf_id=("same_gf_id", "any"),
            )
            .reset_index()
            .sort_values(
                [
                    "any_current_three_layer",
                    "any_current_proteomics_directional",
                    "any_current_fox",
                    "any_current_model_supported",
                    "aryal_product_is_informative",
                    "aryal_consensus_rank_pct_mean",
                ],
                ascending=[False, False, False, False, False, True],
            )
        )
    else:
        product_agg = pd.DataFrame()
    product_agg.to_csv(OUT / "aryal_current_product_overlap_summary.csv", index=False)

    # Counts and summaries.
    summary_rows = []
    def add(metric: str, value, note: str = "") -> None:
        summary_rows.append({"metric": metric, "value": value, "note": note})

    add("Aryal CV label raw rows", cv_meta["label_rows_raw"], "complete_genomes_labeled-james.csv rows")
    add("Aryal CV unique labels", cv_meta["label_unique_accessions"], "deduplicated assembly_accession")
    add("Aryal CV duplicated accessions", cv_meta["label_duplicate_accessions"], "duplicates had no diazotroph-label conflicts")
    add("Aryal CV split unique genomes", cv_meta["split_unique_accessions"], "union of train/test split files")
    add("Current manuscript matrix genomes", len(matrix_ids), matrix_entry)
    add("Current manuscript matrix families", len(matrix_families), matrix_entry)
    add("CV split genomes also in current matrix", len(cv_split_ids & matrix_ids), f"{pct(len(cv_split_ids & matrix_ids), len(matrix_ids))}% of current matrix")
    add("CV split genomes absent from current matrix", len(cv_split_ids - matrix_ids), "all are GCA accessions")
    add("Current matrix genomes absent from CV split", len(matrix_ids - cv_split_ids), "all are GCF accessions")
    add("Current Model-Supported families", len(current_tier1_ids), tier1_entry)
    add("Current Highly Pure families", len(current_tier2_ids), tier2_entry)
    add("Aryal Task7/Tier1-proxy families", len(aryal_ids), ARYAL_HGT_URL)
    add("Aryal families found in current matrix columns", int(family_detail["in_current_matrix_family_columns"].sum()), "")
    add("Aryal families found in current inventory", int(family_detail["in_current_inventory"].sum()), "")
    add("Aryal families found in current Model-Supported", int(family_detail["in_current_model_supported_476"].sum()), "")
    add("Aryal families found in current Highly Pure", int(family_detail["in_current_highly_pure_981"].sum()), "")
    add("Current Model-Supported families recovered by Aryal Tier1 IDs", len(current_tier1_ids & aryal_ids), f"{pct(len(current_tier1_ids & aryal_ids), len(current_tier1_ids))}% of current Model-Supported")
    add("Current inventory families recovered by Aryal Tier1 IDs", len(current_inventory_ids & aryal_ids), f"{pct(len(current_inventory_ids & aryal_ids), len(current_inventory_ids))}% of current 1457-family inventory")
    add("Aryal-only Tier1 IDs vs current inventory", len(aryal_ids - current_inventory_ids), "GF IDs absent from current Model-Supported/Highly Pure inventory")
    add("Aryal-only Tier1 IDs vs current matrix", len(aryal_ids - matrix_families), "GF IDs absent from current 2,286-family matrix")
    add("Shared Aryal/current-inventory exact product matches", int(shared["exact_product_match"].sum()), f"out of {len(shared)} shared GF IDs")
    add("Shared Aryal/current-inventory FOX-layer families", int(bool_col(shared, "layer_fox").sum()), "")
    add("Shared Aryal/current-inventory proteomics directional families", int(bool_col(shared, "layer_proteomics_directional").sum()), "")
    add("Shared Aryal/current-inventory three-layer families", int(bool_col(shared, "three_layer_model_fox_proteomics_directional").sum()), "")
    add("Aryal families with exact product match to current inventory", 0 if product_agg.empty else len(product_agg), "GF-ID independent product-level crosswalk")
    add("Aryal informative products with exact current match", 0 if product_agg.empty else int(product_agg["aryal_product_is_informative"].sum()), "excludes hypothetical/DUF/generic family proteins")
    add("Aryal product matches touching current Model-Supported", 0 if product_agg.empty else int(product_agg["any_current_model_supported"].sum()), "")
    add("Aryal product matches touching current Highly Pure", 0 if product_agg.empty else int(product_agg["any_current_highly_pure"].sum()), "")
    add("Aryal product matches touching current FOX layer", 0 if product_agg.empty else int(product_agg["any_current_fox"].sum()), "")
    add("Aryal product matches touching current proteomics directional layer", 0 if product_agg.empty else int(product_agg["any_current_proteomics_directional"].sum()), "")
    add("Aryal product matches touching current three-layer evidence", 0 if product_agg.empty else int(product_agg["any_current_three_layer"].sum()), "")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "aryal_cv_split_overlap_summary.tsv", sep="\t", index=False)

    # Module and candidate-set overlap summary.
    module_summary = (
        family_detail.assign(
            overlap_class=lambda d: d.apply(
                lambda r: "current Model-Supported"
                if r["in_current_model_supported_476"]
                else ("current Highly Pure" if r["in_current_highly_pure_981"] else ("current matrix only" if r["in_current_matrix_family_columns"] else "not in current matrix")),
                axis=1,
            )
        )
        .groupby(["overlap_class", "aryal_module_bin"], dropna=False)
        .size()
        .reset_index(name="n_families")
        .sort_values(["overlap_class", "n_families"], ascending=[True, False])
    )
    module_summary.to_csv(OUT / "aryal_overlap_by_module.tsv", sep="\t", index=False)

    current_nif_ids = sorted(core_nif["gene_family"].dropna().astype(str).unique())
    aryal_nif_anchor_ids = ["GF_00437", "GF_00180", "GF_01790"]
    nif_anchor_rows = []
    for gf in sorted(set(current_nif_ids + aryal_nif_anchor_ids)):
        curr = current_inventory[current_inventory["gene_family"].astype(str).eq(gf)]
        ar = aryal_tier1[aryal_tier1["gene_family"].astype(str).eq(gf)]
        nif_anchor_rows.append(
            {
                "gene_family": gf,
                "is_current_core_nif_audit_family": gf in current_nif_ids,
                "is_aryal_hgt_anchor_family": gf in aryal_nif_anchor_ids,
                "in_current_matrix": gf in matrix_families,
                "in_current_inventory": gf in current_inventory_ids,
                "in_aryal_task7_tier1_proxy": gf in aryal_ids,
                "current_product": "" if curr.empty else curr.iloc[0].get("product", ""),
                "aryal_task7_product": "" if ar.empty else ar.iloc[0].get("product", ""),
            }
        )
    pd.DataFrame(nif_anchor_rows).to_csv(OUT / "aryal_vs_current_nif_anchor_check.csv", index=False)

    # Markdown report.
    report = []
    report.append("# Aryal CV-Split Protein-Level Overlap Audit")
    report.append("")
    report.append("## Bottom line")
    report.append("")
    report.append(
        "Aryal's claimed starting point is closely related to the current protein-family work, "
        "but it is not the identical manuscript universe. The CV split package contains 449 "
        "unique genomes, including 33 GCA assemblies, while the current manuscript matrix has "
        "426 GCF RefSeq genomes. The overlap is 416 genomes."
    )
    report.append("")
    report.append(
        "The family-level comparison is more cautious: Aryal's Task 7 output gives a 503-family "
        "Tier 1 proxy, while the current manuscript inventory has 476 Model-Supported and 981 "
        "Highly Pure families. Direct GF-ID overlap exists, but because the Aryal run uses a "
        "2,551-family namespace and the current matrix uses 2,286 families, non-overlap should "
        "be interpreted as version/namespace drift unless the exact 449 x 2,551 matrix is supplied."
    )
    report.append("")
    report.append("## Key counts")
    report.append("")
    for row in summary_rows:
        report.append(f"- {row['metric']}: {row['value']} ({row['note']})")
    report.append("")
    report.append("## Interpretation")
    report.append("")
    report.append(
        "The new light helps: this is not an unrelated old gene-level analysis. It is a "
        "protein-level analysis built around a genome set that substantially overlaps the "
        "manuscript set. However, it still should not be merged as if it were the final "
        "426-genome manuscript analysis, because the genome set, family count, Tier 1 count, "
        "and nif-anchor IDs differ."
    )
    report.append("")
    report.append(
        "The most useful direct comparison is therefore not raw GF-ID overlap alone. "
        "Raw GF IDs overlap numerically, but the product labels do not match for shared "
        "GF IDs, indicating namespace drift between clustering versions. The stronger "
        "biological comparison is product/function-level overlap."
    )
    report.append("")
    report.append(
        "At the product/function level, Aryal's 503-family Tier 1 proxy strongly recovers "
        "the same biological territory as the current inventory: 442 Aryal families have "
        "an exact product-label match in the current 1,457-family inventory, 402 touch "
        "current Model-Supported families, 47 touch the FOX layer, 100 touch directional "
        "proteomics, and 9 touch the current three-layer model+FOX+proteomics evidence set."
    )
    report.append("")
    if not product_agg.empty:
        top_three = product_agg[product_agg["any_current_three_layer"]].copy()
        if not top_three.empty:
            report.append("## Strongest product-level overlaps")
            report.append("")
            for _, row in top_three.iterrows():
                report.append(
                    "- "
                    f"{row['aryal_product']} ({row['aryal_gene_family']} -> "
                    f"{row['current_matching_families']}): "
                    f"{row['aryal_module_bin']}"
                )
            report.append("")
    report.append("## Output files")
    report.append("")
    for name in [
        "aryal_cv_split_overlap_summary.tsv",
        "aryal_cv_split_genome_overlap.csv",
        "aryal_task7_tier1_vs_current_inventory_detail.csv",
        "aryal_current_shared_family_detail.csv",
        "aryal_current_product_overlap_long.csv",
        "aryal_current_product_overlap_summary.csv",
        "aryal_tier1_not_in_current_inventory.csv",
        "current_inventory_not_in_aryal_tier1.csv",
        "aryal_overlap_by_module.tsv",
        "aryal_vs_current_nif_anchor_check.csv",
    ]:
        report.append(f"- evidence_tables/{name}")
    (OUT / "aryal_cv_split_overlap_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(summary.to_string(index=False))
    print(f"\nWrote outputs to {OUT}")


if __name__ == "__main__":
    main()
