from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "literature_proteomics_supplement_tables"
FEATURES = SUPP / "feature_tables"
OUT = ROOT / "evidence_tables" / "literature_proteomics_family_mapping"

GF_MAP = ROOT / "_pangenome_data" / "genome_protein_family_map.tsv"
INVENTORY = ROOT / "evidence_tables" / "morphotype_bridge_family_inventory.tsv"
CYANO_CROSSWALK = ROOT / "data" / "cyanothece_old_to_wp_crosswalk.tsv"
SOURCE_MAP = SUPP / "source_map.csv"


ASSEMBLIES = {
    "Panda_2025_Crocosphaera_ATCC51142": "GCF_000017845.1",
    "Welkie_2014_Cyanothece_PCC7822": "GCF_000147335.1",
    "Sandh_2014_Nostoc_punctiforme": "GCF_000020025.1",
    "Held_2022_Trichodesmium_IMS101": "GCF_000014265.1",
}


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def first_present(row: pd.Series, names: list[str]) -> str:
    for name in names:
        if name in row and clean_text(row[name]):
            return clean_text(row[name])
    return ""


def parse_old_locus(attributes: object) -> str:
    text = clean_text(attributes)
    match = re.search(r"(?:^|;)old_locus_tag=([^;]+)", text)
    return match.group(1) if match else ""


def parse_quoted_locus(text: object, prefix_re: str) -> str:
    match = re.search(prefix_re, clean_text(text))
    return match.group(1) if match else ""


def parse_ref_accession(text: object) -> str:
    # Supports ref|YP_003899811.1| and bare WP_ / YP_ / ACB-style accessions.
    s = clean_text(text)
    match = re.search(r"(WP_\d+\.\d+|YP_\d+\.\d+|[A-Z]{3}\d+\.\d+|[A-Z0-9]{6,10})", s)
    return match.group(1) if match else ""


def load_current_family_map() -> pd.DataFrame:
    gf = pd.read_csv(GF_MAP, sep="\t", dtype=str).fillna("")
    gf = gf.drop_duplicates(["gene_family", "genome_accession", "protein_accession"])
    return gf


def load_inventory() -> pd.DataFrame:
    inv = pd.read_csv(INVENTORY, sep="\t", dtype=str).fillna("")
    keep = [
        "gene_family",
        "candidate_set",
        "product",
        "module_bin",
        "fox_status",
        "best_rank",
        "consensus_rank_pct_mean",
        "diazotroph_pct_mean",
        "n_member_genomes",
        "n_member_proteins",
        "nif_or_nitrogenase_audit",
    ]
    return inv[[c for c in keep if c in inv.columns]].drop_duplicates("gene_family")


def load_feature_bridge() -> pd.DataFrame:
    rows: list[pd.DataFrame] = []

    # RefSeq feature tables for non-51142 literature organisms.
    for path in FEATURES.glob("GCF_*_feature_table.txt.gz"):
        with gzip.open(path, "rt", errors="replace") as handle:
            df = pd.read_csv(handle, sep="\t", dtype=str).fillna("")
        df.columns = [c.lstrip("# ").strip() for c in df.columns]
        genes = df[df["feature"] == "gene"][
            ["assembly", "locus_tag", "symbol", "attributes"]
        ].copy()
        genes = genes.rename(columns={"locus_tag": "refseq_locus_tag"})
        genes["old_locus_tag"] = genes["attributes"].map(parse_old_locus)
        cds = df[(df["feature"] == "CDS") & (df["product_accession"] != "")][
            ["assembly", "locus_tag", "product_accession", "name"]
        ].copy()
        cds = cds.rename(
            columns={
                "locus_tag": "refseq_locus_tag",
                "product_accession": "wp_accession",
                "name": "refseq_product",
            }
        )
        bridge = genes.merge(cds, on=["assembly", "refseq_locus_tag"], how="inner")
        bridge["mapping_bridge_source"] = path.name
        rows.append(bridge)

    # Existing ATCC 51142 GenBank old locus/accession -> RefSeq WP bridge.
    cce = pd.read_csv(CYANO_CROSSWALK, sep="\t", dtype=str).fillna("")
    cce = cce.rename(
        columns={
            "wp_accession": "wp_accession",
            "refseq_product": "refseq_product",
            "old_accession": "old_accession",
        }
    )
    cce["assembly"] = "GCF_000017845.1"
    cce["symbol"] = ""
    cce["mapping_bridge_source"] = "data/cyanothece_old_to_wp_crosswalk.tsv"
    rows.append(
        cce[
            [
                "assembly",
                "old_locus_tag",
                "refseq_locus_tag",
                "symbol",
                "old_accession",
                "wp_accession",
                "refseq_product",
                "mapping_bridge_source",
            ]
        ]
    )

    bridge = pd.concat(rows, ignore_index=True).fillna("")
    for col in [
        "assembly",
        "old_locus_tag",
        "refseq_locus_tag",
        "symbol",
        "old_accession",
        "wp_accession",
        "refseq_product",
        "mapping_bridge_source",
    ]:
        if col not in bridge:
            bridge[col] = ""
    return bridge.drop_duplicates(
        ["assembly", "old_locus_tag", "refseq_locus_tag", "old_accession", "wp_accession"]
    )


def map_rows_to_families(norm: pd.DataFrame, bridge: pd.DataFrame, gf: pd.DataFrame) -> pd.DataFrame:
    # Build long-form bridge keys so each source protein can match by old locus,
    # RefSeq locus, old accession, or WP accession without product-name matching.
    key_rows: list[dict[str, str]] = []
    for _, row in bridge.iterrows():
        assembly = row["assembly"]
        keys = [
            ("old_locus_tag", row.get("old_locus_tag", "")),
            ("refseq_locus_tag", row.get("refseq_locus_tag", "")),
            ("old_accession", row.get("old_accession", "")),
            ("wp_accession", row.get("wp_accession", "")),
        ]
        for key_type, key_value in keys:
            key_value = clean_text(key_value)
            if key_value:
                key_rows.append(
                    {
                        "assembly_accession": assembly,
                        "lookup_key_type": key_type,
                        "lookup_key": key_value,
                        "wp_accession": row.get("wp_accession", ""),
                        "refseq_locus_tag": row.get("refseq_locus_tag", ""),
                        "refseq_product": row.get("refseq_product", ""),
                        "mapping_bridge_source": row.get("mapping_bridge_source", ""),
                    }
                )
    bridge_keys = pd.DataFrame(key_rows).drop_duplicates()

    gf_for_join = gf.rename(
        columns={
            "genome_accession": "assembly_accession",
            "protein_accession": "wp_accession",
        }
    )

    all_hits: list[pd.DataFrame] = []
    key_priority = [
        ("wp_accession", "wp_accession"),
        ("old_locus_tag", "old_locus_tag"),
        ("refseq_locus_tag", "refseq_locus_tag"),
        ("protein_accession_original", "old_accession"),
    ]
    for source_col, key_type in key_priority:
        tmp = norm[norm[source_col].map(clean_text) != ""].copy()
        if tmp.empty:
            continue
        tmp["lookup_key_type"] = key_type
        tmp["lookup_key"] = tmp[source_col].map(clean_text)
        hit = tmp.merge(
            bridge_keys,
            on=["assembly_accession", "lookup_key_type", "lookup_key"],
            how="left",
            suffixes=("", "_bridge"),
        )
        # For locus/accession keys, the usable WP accession comes from the bridge.
        # Keep a direct source WP accession only when it was already present.
        for col in ["wp_accession", "refseq_locus_tag", "refseq_product", "mapping_bridge_source"]:
            bridge_col = f"{col}_bridge"
            if bridge_col in hit.columns:
                hit[col] = hit[bridge_col].where(hit[bridge_col] != "", hit[col])
        hit = hit.merge(
            gf_for_join[["assembly_accession", "wp_accession", "gene_family"]],
            on=["assembly_accession", "wp_accession"],
            how="left",
        )
        hit["mapping_method"] = key_type
        all_hits.append(hit)

    if not all_hits:
        norm["map_status"] = "no_lookup_keys"
        return norm

    hits = pd.concat(all_hits, ignore_index=True).fillna("")
    hits = hits[hits["wp_accession"] != ""].copy()
    hits["has_gene_family"] = hits["gene_family"] != ""

    # Prefer direct WP hits, then old locus, RefSeq locus, old accession.
    rank = {"wp_accession": 1, "old_locus_tag": 2, "refseq_locus_tag": 3, "old_accession": 4}
    hits["method_rank"] = hits["mapping_method"].map(rank).fillna(99)
    hits = hits.sort_values(["source_record_id", "has_gene_family", "method_rank"], ascending=[True, False, True])
    hits = hits.drop_duplicates(
        [
            "source_record_id",
            "wp_accession",
            "gene_family",
        ],
        keep="first",
    )

    # Attach compact mapping summaries to one row per source record.
    agg = (
        hits.groupby("source_record_id")
        .agg(
            mapped_wp_accessions=("wp_accession", lambda s: ";".join(sorted(set(x for x in s if x)))),
            mapped_refseq_locus_tags=("refseq_locus_tag", lambda s: ";".join(sorted(set(x for x in s if x)))),
            mapped_gene_families=("gene_family", lambda s: ";".join(sorted(set(x for x in s if x)))),
            mapping_methods=("mapping_method", lambda s: ";".join(sorted(set(x for x in s if x)))),
            mapping_bridge_sources=("mapping_bridge_source", lambda s: ";".join(sorted(set(x for x in s if x)))),
            refseq_products=("refseq_product", lambda s: ";".join(sorted(set(x for x in s if x))[:3])),
        )
        .reset_index()
    )
    out = norm.merge(agg, on="source_record_id", how="left").fillna("")
    out["map_status"] = out.apply(
        lambda r: "mapped_to_family"
        if r["mapped_gene_families"]
        else ("mapped_to_wp_not_in_filtered_family_atlas" if r["mapped_wp_accessions"] else "unmapped"),
        axis=1,
    )
    return out


def normalize_panda() -> pd.DataFrame:
    path = SUPP / "Panda_2025_Crocosphaera_ATCC51142_nitrogen_fixation_supp_table.xlsx"
    df = pd.read_excel(path, sheet_name="Table S3", header=1, dtype=str).fillna("")
    rows = []
    for i, row in df.iterrows():
        locus = first_present(row, ["Locus tag"])
        if not locus:
            continue
        metrics = []
        for col in [
            "Two-way ANOVA p value nitrate",
            "Two-way ANOVA p value lightordark",
            "Two-way ANOVA p value Interaction",
            "D-_avgLFQ",
            "D+_avgLFQ",
            "L-_avgLFQ",
            "L+_avgLFQ",
            "D-_D+(logFC)",
            "L-_L+(logFC)",
        ]:
            if col in df.columns and clean_text(row[col]):
                metrics.append(f"{col}={clean_text(row[col])}")
        rows.append(
            {
                "study_id": "Panda_2025_Crocosphaera_ATCC51142",
                "paper_short": "Panda et al. 2025",
                "organism": "Crocosphaera subtropica ATCC 51142",
                "assembly_accession": ASSEMBLIES["Panda_2025_Crocosphaera_ATCC51142"],
                "source_file": path.name,
                "source_sheet": "Table S3",
                "source_record_id": f"Panda2025|{locus}",
                "source_protein_key": locus,
                "old_locus_tag": locus,
                "refseq_locus_tag": "",
                "protein_accession_original": first_present(row, ["Protein IDs"]),
                "wp_accession": "",
                "gene_symbol": first_present(row, ["Gene name"]),
                "protein_name": first_present(row, ["Protein name"]),
                "evidence_metric": "; ".join(metrics),
            }
        )
    return pd.DataFrame(rows)


def normalize_welkie() -> pd.DataFrame:
    path = SUPP / "Welkie_2014_Cyanothece_PCC7822_TableS4_protein.xlsx"
    df = pd.read_excel(path, sheet_name="total_protein values", dtype=str).fillna("")
    rows = []
    for _, row in df.iterrows():
        locus = first_present(row, ["JGI"])
        if not locus:
            continue
        accession = parse_ref_accession(first_present(row, ["PNames", "Protein"]))
        metrics = []
        for col in ["CT_D0", "CT_D3", "CT_L0", "CT_L3"]:
            if col in df.columns and clean_text(row[col]):
                metrics.append(f"{col}={clean_text(row[col])}")
        rows.append(
            {
                "study_id": "Welkie_2014_Cyanothece_PCC7822",
                "paper_short": "Welkie et al. 2014",
                "organism": "Gloeothece verrucosa PCC 7822",
                "assembly_accession": ASSEMBLIES["Welkie_2014_Cyanothece_PCC7822"],
                "source_file": path.name,
                "source_sheet": "total_protein values",
                "source_record_id": f"Welkie2014|{locus}",
                "source_protein_key": locus,
                "old_locus_tag": locus,
                "refseq_locus_tag": "",
                "protein_accession_original": accession,
                "wp_accession": "",
                "gene_symbol": "",
                "protein_name": first_present(row, ["PNames"]),
                "evidence_metric": "; ".join(metrics),
            }
        )
    return pd.DataFrame(rows)


def normalize_sandh() -> pd.DataFrame:
    files = [
        ("TableS1", SUPP / "Sandh_2014_Nostoc_punctiforme_heterocyst_proteome_TableS1.xlsx"),
        ("TableS2", SUPP / "Sandh_2014_Nostoc_punctiforme_heterocyst_proteome_TableS2.xlsx"),
        ("TableS4", SUPP / "Sandh_2014_Nostoc_punctiforme_heterocyst_proteome_TableS4.xlsx"),
    ]
    rows = []
    for label, path in files:
        raw = pd.read_excel(path, sheet_name=0, header=None, dtype=str).fillna("")
        # These sheets have a title row, then the actual header row.
        header = [clean_text(x) for x in raw.iloc[1].tolist()]
        df = raw.iloc[2:].copy()
        df.columns = header
        for _, row in df.iterrows():
            # The Sandh supplemental column named "RefSeq Locus Tag" contains
            # legacy Npun_F/R locus tags. The downloaded RefSeq feature table
            # stores those as old_locus_tag and uses NPUN_RS identifiers as
            # current RefSeq locus tags.
            old_npun_locus = first_present(row, ["RefSeq Locus Tag"])
            patric_locus = first_present(row, ["PATRIC Locus Tag"])
            if not old_npun_locus and not patric_locus:
                continue
            source_key = old_npun_locus or patric_locus
            metrics = []
            for col in [
                "Log2 Ratio (Heterocyst/Filaments)",
                "24h Log2 Ratio (Het/Fil) (Present study)",
                "Variability [%]",
                "Peptide Pair Count",
                "Steady state Log2 ratio [7]",
            ]:
                if col in df.columns and clean_text(row[col]):
                    metrics.append(f"{col}={clean_text(row[col])}")
            rows.append(
                {
                    "study_id": "Sandh_2014_Nostoc_punctiforme",
                    "paper_short": "Sandh et al. 2014",
                    "organism": "Nostoc punctiforme PCC 73102",
                    "assembly_accession": ASSEMBLIES["Sandh_2014_Nostoc_punctiforme"],
                    "source_file": path.name,
                    "source_sheet": label,
                    "source_record_id": f"Sandh2014|{label}|{source_key}",
                    "source_protein_key": source_key,
                    "old_locus_tag": old_npun_locus,
                    "refseq_locus_tag": "",
                    "protein_accession_original": first_present(row, ["GenBank Accession"]),
                    "wp_accession": "",
                    "gene_symbol": first_present(row, ["Gene Name", "GENE"]),
                    "protein_name": first_present(row, ["Name-annotation [27]"]),
                    "evidence_metric": "; ".join(metrics),
                }
            )
    return pd.DataFrame(rows)


def normalize_held() -> pd.DataFrame:
    path = SUPP / "Held_2022_Trichodesmium_dielproteindata_BCODMO_matrix.csv"
    raw = pd.read_csv(path, dtype=str).fillna("")
    hour_cols = [c for c in raw.columns if c not in {"hourspostdawn"}]
    rows = []
    for _, row in raw.iterrows():
        protein = clean_text(row["hourspostdawn"])
        if not protein or protein.lower() == "cnratio":
            continue
        locus = parse_quoted_locus(protein, r"\b(Tery_\d+)\b")
        if not locus:
            continue
        vals = pd.to_numeric(row[hour_cols], errors="coerce")
        finite = vals.dropna()
        if finite.empty:
            metric = ""
        else:
            max_hour = finite.idxmax()
            min_hour = finite.idxmin()
            metric = (
                f"mean_relative_abundance={finite.mean():.4g}; "
                f"max={finite.max():.4g}@{max_hour}h; "
                f"min={finite.min():.4g}@{min_hour}h; "
                f"dynamic_range={finite.max() - finite.min():.4g}"
            )
        rows.append(
            {
                "study_id": "Held_2022_Trichodesmium_IMS101",
                "paper_short": "Held et al. 2022",
                "organism": "Trichodesmium erythraeum IMS101",
                "assembly_accession": ASSEMBLIES["Held_2022_Trichodesmium_IMS101"],
                "source_file": path.name,
                "source_sheet": "dielproteindata",
                "source_record_id": f"Held2022|{locus}",
                "source_protein_key": locus,
                "old_locus_tag": locus,
                "refseq_locus_tag": "",
                "protein_accession_original": "",
                "wp_accession": "",
                "gene_symbol": "",
                "protein_name": protein,
                "evidence_metric": metric,
            }
        )
    return pd.DataFrame(rows)


def build_family_evidence(mapped: pd.DataFrame, inv: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for _, row in mapped.iterrows():
        families = [x for x in clean_text(row.get("mapped_gene_families", "")).split(";") if x]
        for gf in families:
            rows.append(
                {
                    "study_id": row["study_id"],
                    "paper_short": row["paper_short"],
                    "organism": row["organism"],
                    "gene_family": gf,
                    "source_record_id": row["source_record_id"],
                    "source_protein_key": row["source_protein_key"],
                    "mapped_wp_accessions": row["mapped_wp_accessions"],
                    "protein_name": row["protein_name"],
                    "evidence_metric": row["evidence_metric"],
                    "mapping_methods": row["mapping_methods"],
                }
            )
    long = pd.DataFrame(rows)
    if long.empty:
        return long, long
    fam_study = (
        long.groupby(["gene_family", "study_id", "paper_short", "organism"])
        .agg(
            n_unique_source_proteins=("source_record_id", "nunique"),
            source_protein_keys=("source_protein_key", lambda s: ";".join(sorted(set(clean_text(x) for x in s if clean_text(x)))[:20])),
            mapped_wp_accessions=("mapped_wp_accessions", lambda s: ";".join(sorted(set(";".join(s).split(";")))[:20]).strip(";")),
            representative_protein_names=("protein_name", lambda s: " | ".join(sorted(set(clean_text(x) for x in s if clean_text(x)))[:3])),
            evidence_metric_examples=("evidence_metric", lambda s: " | ".join([clean_text(x) for x in s if clean_text(x)][:3])),
            mapping_methods=("mapping_methods", lambda s: ";".join(sorted(set(";".join(s).split(";")))[:10]).strip(";")),
        )
        .reset_index()
    )
    fam_study = fam_study.merge(inv, on="gene_family", how="left")

    fam_summary = (
        fam_study.groupby("gene_family")
        .agg(
            n_literature_studies=("study_id", "nunique"),
            studies=("paper_short", lambda s: "; ".join(sorted(set(s)))),
            organisms=("organism", lambda s: "; ".join(sorted(set(s)))),
            n_study_family_rows=("study_id", "size"),
            n_unique_source_proteins_total=("n_unique_source_proteins", "sum"),
            source_protein_keys=("source_protein_keys", lambda s: ";".join(sorted(set(";".join(s).split(";")))[:30]).strip(";")),
            representative_protein_names=("representative_protein_names", lambda s: " | ".join(sorted(set(" | ".join(s).split(" | ")))[:4]).strip(" |")),
        )
        .reset_index()
        .merge(inv, on="gene_family", how="left")
    )
    fam_summary = fam_summary.sort_values(
        ["n_literature_studies", "candidate_set", "gene_family"],
        ascending=[False, True, True],
    )
    return fam_study, fam_summary


def build_mapping_summary(mapped: pd.DataFrame) -> pd.DataFrame:
    return (
        mapped.groupby(["study_id", "paper_short", "organism", "map_status"])
        .agg(n_source_proteins=("source_record_id", "nunique"))
        .reset_index()
        .sort_values(["study_id", "map_status"])
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    gf = load_current_family_map()
    inv = load_inventory()
    bridge = load_feature_bridge()

    normalized = pd.concat(
        [
            normalize_panda(),
            normalize_welkie(),
            normalize_sandh(),
            normalize_held(),
        ],
        ignore_index=True,
    ).fillna("")
    # One row per source protein key per study/source sheet. This removes accidental duplicate exports.
    normalized = normalized.drop_duplicates(
        ["study_id", "source_sheet", "source_record_id"], keep="first"
    )

    mapped = map_rows_to_families(normalized, bridge, gf)
    fam_study, fam_summary = build_family_evidence(mapped, inv)
    mapping_summary = build_mapping_summary(mapped)
    unmapped = mapped[mapped["map_status"] != "mapped_to_family"].copy()
    source_map = pd.read_csv(SOURCE_MAP, dtype=str).fillna("")

    bridge_out = bridge.merge(
        gf.rename(
            columns={
                "genome_accession": "assembly",
                "protein_accession": "wp_accession",
            }
        )[["assembly", "wp_accession", "gene_family"]],
        on=["assembly", "wp_accession"],
        how="left",
    ).fillna("")

    outputs = {
        "normalized_proteins.csv": mapped,
        "family_evidence_by_study.csv": fam_study,
        "family_summary.csv": fam_summary,
        "mapping_summary.csv": mapping_summary,
        "unmapped_or_unfiltered_proteins.csv": unmapped,
        "accession_bridge_used.csv": bridge_out,
        "source_map.csv": source_map,
    }
    for name, df in outputs.items():
        df.to_csv(OUT / name, index=False)

    readme = OUT / "README.txt"
    readme.write_text(
        "\n".join(
            [
                "Literature proteomics to protein-family mapping",
                "",
                "Mapping is accession/locus based only: old locus tag, RefSeq locus tag, old accession, or WP accession -> NCBI feature/crosswalk bridge -> current genome_protein_family_map.tsv.",
                "Product-name/fuzzy matching is intentionally not used in the mapped family evidence sheets.",
                "Rows are deduplicated to one source protein record before family summarization, and family_evidence_by_study.csv is one row per study + gene_family.",
                "Proteins with a WP bridge but no gene_family are marked mapped_to_wp_not_in_filtered_family_atlas; these are real proteins but absent from the filtered 2,286-family atlas.",
                "Aryal 2011 binary .xls files are recorded in source_map.csv but not parsed here because this runtime lacks a reliable .xls reader.",
            ]
        ),
        encoding="utf-8",
    )

    print(f"wrote {OUT}")
    print("normalized proteins:", len(mapped))
    print("mapped to family:", int((mapped["map_status"] == "mapped_to_family").sum()))
    print("family-study rows:", len(fam_study))
    print("family summary rows:", len(fam_summary))
    print(mapping_summary.to_string(index=False))


if __name__ == "__main__":
    main()
