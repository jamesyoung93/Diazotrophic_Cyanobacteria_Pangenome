from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence_tables"
LIT = EVIDENCE / "literature_proteomics_family_mapping"
OUT = EVIDENCE / "composite_proteomics_priority"


TRUE_VALUES = {"true", "1", "yes", "y"}
FALSE_VALUES = {"false", "0", "no", "n", ""}


STUDY_CONTEXTS = {
    "Panda et al. 2025": {
        "context": "unicellular Crocosphaera/ATCC51142 N-state and light/dark proteomics",
        "context_group": "unicellular",
        "physiology": "nitrogen-fixing physiology with nitrate-replete comparison",
        "organism_group": "unicellular",
    },
    "Welkie et al. 2014": {
        "context": "unicellular Gloeothece/Cyanothece PCC 7822 diel N-free proteomics",
        "context_group": "unicellular",
        "physiology": "diel nitrogen-fixing physiology",
        "organism_group": "unicellular",
    },
    "Sandh et al. 2014": {
        "context": "Nostoc heterocyst-enriched versus vegetative/filament proteomics",
        "context_group": "heterocyst",
        "physiology": "cell-type resolved nitrogen-fixing physiology",
        "organism_group": "heterocystous filament",
    },
    "Held et al. 2022": {
        "context": "Trichodesmium diel daytime proteome with daytime nitrogenase activity",
        "context_group": "nonheterocystous filament",
        "physiology": "daytime nonheterocystous nitrogen fixation",
        "organism_group": "nonheterocystous filament",
    },
}


STORY_KEYWORDS = [
    (
        "Core nitrogenase / nif machinery",
        [
            "nitrogenase",
            "nif",
            "hesa",
            "hesc",
            "homocitrate",
            "ferredoxin ii",
            "fdxn",
        ],
    ),
    (
        "Nitrogen assimilation and regulation",
        [
            "glutamate synthase",
            "gogat",
            "glutamine",
            "ammonium",
            "nitrate",
            "nitrite",
            "urea",
            "ntca",
            "p-ii",
            "pii",
            "nitrogen regulatory",
        ],
    ),
    (
        "Metallocluster and cofactor support",
        [
            "fe-s",
            "iron-sulfur",
            "suf",
            "isc",
            "ferredoxin",
            "flavodoxin",
            "molyb",
            "cobalamin",
            "cob",
            "radical sam",
            "biotin",
            "thiamine",
            "cofactor",
            "co-chaperone",
            "metallochaperone",
        ],
    ),
    (
        "Respiration and bioenergetics",
        [
            "cytochrome",
            "oxidase",
            "respiration",
            "respiratory",
            "nadph dehydrogenase",
            "ndh",
            "atp synthase",
            "quinone",
            "electron transport",
            "photosystem",
            "ferredoxin--nadp",
        ],
    ),
    (
        "Oxygen, redox, and stress protection",
        [
            "superoxide",
            "peroxiredoxin",
            "thioredoxin",
            "glutaredoxin",
            "rubrerythrin",
            "peroxidase",
            "oxidative",
            "redox",
            "stress",
            "detox",
            "ars",
        ],
    ),
    (
        "Carbon flux and reductant supply",
        [
            "transketolase",
            "glyceraldehyde",
            "gapdh",
            "aldolase",
            "fructose",
            "glucose-6-phosphate",
            "isocitrate",
            "pyruvate",
            "ketoacid",
            "ketol-acid",
            "dehydrogenase",
            "glycogen",
            "glycolate",
            "pentose phosphate",
            "ribulose",
            "rubisco",
            "carbon",
        ],
    ),
    (
        "Transport and envelope remodeling",
        [
            "transporter",
            "permease",
            "porin",
            "efflux",
            "abc",
            "membrane",
            "secretion",
            "glycosyltransferase",
            "capsule",
            "polysaccharide",
        ],
    ),
    (
        "Proteostasis, translation, and repair",
        [
            "ribosomal",
            "translation",
            "trna",
            "elongation",
            "chaperone",
            "chaperonin",
            "heat shock",
            "protease",
            "clp",
            "dna repair",
            "reca",
            "uvr",
            "helicase",
            "peptidase",
        ],
    ),
    ("Regulatory / signaling", ["regulator", "response regulator", "histidine kinase", "ggdef", "eal", "sigma", "transcription"]),
]


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def boolish(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    text = str(value).strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    return bool(text)


def num(value, default=0.0) -> float:
    try:
        if value is None or value == "":
            return default
        out = float(value)
        if math.isnan(out):
            return default
        return out
    except (TypeError, ValueError):
        return default


def listify_semicolon(value) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return [part.strip() for part in str(value).split(";") if part.strip()]


def first_present(row: pd.Series, names: list[str], default=""):
    for name in names:
        if name in row and pd.notna(row[name]) and str(row[name]) != "":
            return row[name]
    return default


def metric_value(metric_text: object, label: str, default=float("nan")) -> float:
    text = "" if metric_text is None else str(metric_text)
    match = re.search(re.escape(label) + r"\s*=\s*([^;|]+)", text)
    if not match:
        return default
    return num(match.group(1), default=default)


def metric_hour(metric_text: object, label: str, default=float("nan")) -> float:
    text = "" if metric_text is None else str(metric_text)
    match = re.search(re.escape(label) + r"\s*=\s*[^;|@]+@([0-9.]+)h", text)
    if not match:
        return default
    return num(match.group(1), default=default)


def in_trichodesmium_day_window(hour: float) -> bool:
    return not math.isnan(hour) and 6 <= hour <= 15


def in_trichodesmium_transition_window(hour: float) -> bool:
    return not math.isnan(hour) and hour in {4, 18}


def classify_literature_nfix_response(row: pd.Series) -> dict[str, object]:
    """Classify mapped literature proteomics rows as N-fix-responsive evidence.

    This intentionally excludes simple detection/mapping. The four downloaded
    literature datasets encode different comparisons, so each receives an
    explicit response criterion:
    - Panda: nitrate-state or interaction p <= 0.05 plus >=0.25 log2 N-fix/nitrate contrast.
    - Sandh: abs heterocyst/filament log2 ratio >= 1.0.
    - Welkie: abs diel range across reported N-fixing time points >= 1.0.
    - Held: relative diel range >= 0.5 in Trichodesmium N-fixing cultures.
    """
    study = str(row.get("paper_short", ""))
    metric = row.get("evidence_metric", "")
    response = False
    response_type = ""
    response_direction = ""
    response_strength = float("nan")
    response_criterion = ""
    active_phase_response = False
    active_phase_up = False
    active_phase_down = False
    phase_alignment_class = ""
    intended_active_phase = ""
    active_phase_metric = ""
    phase_alignment_note = ""

    if study == "Panda et al. 2025":
        nitrate_p = metric_value(metric, "Two-way ANOVA p value nitrate")
        interaction_p = metric_value(metric, "Two-way ANOVA p value Interaction")
        dark_logfc = metric_value(metric, "D-_D+(logFC)")
        light_logfc = metric_value(metric, "L-_L+(logFC)")
        sig = (
            (not math.isnan(nitrate_p) and nitrate_p <= 0.05)
            or (not math.isnan(interaction_p) and interaction_p <= 0.05)
        )
        dark_hit = sig and not math.isnan(dark_logfc) and abs(dark_logfc) >= 0.25
        light_hit = sig and not math.isnan(light_logfc) and abs(light_logfc) >= 0.25
        logfcs = [x for x in [dark_logfc, light_logfc] if not math.isnan(x)]
        max_abs = max([abs(x) for x in logfcs], default=float("nan"))
        response = sig and (not math.isnan(max_abs) and max_abs >= 0.25)
        response_type = "nitrogen-state differential"
        response_strength = max_abs
        if response:
            if all(x > 0 for x in logfcs):
                response_direction = "higher under N-fixing minus-nitrate conditions"
            elif all(x < 0 for x in logfcs):
                response_direction = "lower under N-fixing minus-nitrate conditions"
            else:
                response_direction = "light/dark-dependent N-fixing response"
        response_criterion = "p_nitrate<=0.05 or p_interaction<=0.05, with max abs D-/D+ or L-/L+ logFC >= 0.25"
        intended_active_phase = "dark nitrate-free Crocosphaera phase; D- versus D+ is the active-phase contrast"
        active_phase_metric = f"D-_D+(logFC)={dark_logfc:.4g}; L-_L+(logFC)={light_logfc:.4g}"
        if dark_hit:
            active_phase_response = True
            active_phase_up = dark_logfc > 0
            active_phase_down = dark_logfc < 0
            if light_hit and dark_logfc * light_logfc < 0:
                phase_alignment_class = "active_phase_D_minus_contrast_with_light_dark_interaction"
            elif light_hit:
                phase_alignment_class = "active_phase_D_minus_contrast_broad_minusN"
            else:
                phase_alignment_class = "active_phase_D_minus_contrast_only"
        elif light_hit:
            phase_alignment_class = "not_scored_light_only_minusN_contrast"
            phase_alignment_note = "Light-only Crocosphaera response retained as metadata but not rewarded."

    elif study == "Sandh et al. 2014":
        ratio = metric_value(metric, "Log2 Ratio (Heterocyst/Filaments)")
        if math.isnan(ratio):
            ratio = metric_value(metric, "24h Log2 Ratio (Het/Fil) (Present study)")
        response = not math.isnan(ratio) and abs(ratio) >= 1.0
        response_type = "heterocyst/filament differential"
        response_strength = abs(ratio) if not math.isnan(ratio) else float("nan")
        if response:
            response_direction = "heterocyst-enriched" if ratio > 0 else "filament/vegetative-enriched"
            active_phase_response = True
            active_phase_up = ratio > 0
            active_phase_down = ratio < 0
            phase_alignment_class = (
                "active_phase_heterocyst_enriched"
                if ratio > 0
                else "not_scored_heterocyst_depleted_supporting_filament"
            )
        intended_active_phase = "heterocyst cell type in N2-fixing Nostoc filaments"
        active_phase_metric = f"heterocyst_filament_log2_ratio={ratio:.4g}"
        if active_phase_down:
            phase_alignment_note = "Filament/vegetative-enriched rows are retained as compartment metadata but not rewarded."
        response_criterion = "abs heterocyst/filament log2 ratio >= 1.0"

    elif study == "Welkie et al. 2014":
        labels = ["CT_D0", "CT_D3", "CT_L0", "CT_L3"]
        values_by_label = {label: metric_value(metric, label) for label in labels}
        finite_by_label = {k: v for k, v in values_by_label.items() if not math.isnan(v)}
        values = list(finite_by_label.values())
        dynamic_range = max(values) - min(values) if values else float("nan")
        response = not math.isnan(dynamic_range) and dynamic_range >= 1.0
        response_type = "diel dynamic in N-fixing culture"
        response_strength = dynamic_range
        response_direction = "diel dynamic" if response else ""
        response_criterion = "dynamic range across reported CT_D/CT_L N-fixing time points >= 1.0"
        intended_active_phase = "late-light/early-dark CT_D0 and CT_D3 active-window proxy"
        if finite_by_label:
            max_label = max(finite_by_label, key=finite_by_label.get)
            min_label = min(finite_by_label, key=finite_by_label.get)
            active_phase_metric = f"max={max_label}; min={min_label}"
            if response and max_label.startswith("CT_D"):
                active_phase_response = True
                active_phase_up = True
                phase_alignment_class = "active_phase_CT_D_window_high"
            elif response and min_label.startswith("CT_D"):
                active_phase_response = True
                active_phase_down = True
                phase_alignment_class = "not_scored_CT_D_window_low"
                phase_alignment_note = "CT_D-window depleted rows are retained as metadata but not rewarded."
            elif response:
                phase_alignment_class = "not_scored_non_CT_D_window_extreme"

    elif study == "Held et al. 2022":
        mean_abundance = metric_value(metric, "mean_relative_abundance")
        dynamic_range = metric_value(metric, "dynamic_range")
        rel_dynamic = (
            dynamic_range / mean_abundance
            if mean_abundance and not math.isnan(mean_abundance) and not math.isnan(dynamic_range)
            else float("nan")
        )
        response = not math.isnan(rel_dynamic) and rel_dynamic >= 0.5
        response_type = "relative diel dynamic in N-fixing culture"
        response_strength = rel_dynamic
        response_direction = "diel dynamic" if response else ""
        response_criterion = "dynamic_range / mean_relative_abundance >= 0.5"
        max_hour = metric_hour(metric, "max")
        min_hour = metric_hour(metric, "min")
        intended_active_phase = "daylight Trichodesmium N2-fixing interval; 6-15 h post-dawn scored as active phase"
        active_phase_metric = f"max_hour={max_hour:g}; min_hour={min_hour:g}"
        if response and in_trichodesmium_day_window(max_hour):
            active_phase_response = True
            active_phase_up = True
            phase_alignment_class = "active_phase_daytime_high"
        elif response and in_trichodesmium_day_window(min_hour):
            active_phase_response = True
            active_phase_down = True
            phase_alignment_class = "not_scored_daytime_low"
            phase_alignment_note = "Daytime-low rows are retained as metadata but not rewarded."
        elif response and (in_trichodesmium_transition_window(max_hour) or in_trichodesmium_transition_window(min_hour)):
            phase_alignment_class = "not_scored_transition_window_dynamic"
            phase_alignment_note = "Transition-window dynamics are retained as metadata but not rewarded."
        elif response:
            phase_alignment_class = "not_scored_night_only_dynamic"
            phase_alignment_note = "Night-only Trichodesmium dynamics are retained as metadata but not rewarded."

    return {
        "nfix_response_evidence": response,
        "nfix_response_type": response_type,
        "nfix_response_direction": response_direction,
        "nfix_response_strength": response_strength,
        "nfix_response_criterion": response_criterion,
        "nfix_active_phase_response": active_phase_response,
        "nfix_active_phase_up": active_phase_up,
        "nfix_active_phase_down": active_phase_down,
        "nfix_phase_alignment_class": phase_alignment_class,
        "nfix_intended_active_phase": intended_active_phase,
        "nfix_active_phase_metric": active_phase_metric,
        "nfix_phase_alignment_note": phase_alignment_note,
    }


def classify_story(product: str, module_bin: str, nif_flag: bool) -> str:
    text = f"{product} {module_bin}".lower()
    if nif_flag:
        return "Core nitrogenase / nif machinery"
    for label, keywords in STORY_KEYWORDS:
        if any(keyword in text for keyword in keywords):
            return label
    if "other" not in str(module_bin).lower() and str(module_bin).strip():
        return str(module_bin)
    return "Other / unassigned"


def is_housekeeping_or_lineage(row: pd.Series) -> bool:
    text = " ".join(
        str(row.get(col, ""))
        for col in ["product", "module_bin"]
    ).lower()
    markers = [
        "lineage marker",
        "housekeeping",
        "ribosomal protein",
        "30s ribosomal",
        "50s ribosomal",
        "ribosome-binding factor",
        "elongation factor",
        "translation initiation factor",
        "translation elongation",
        "aminoacyl-trna synthetase",
        "trna synthetase",
    ]
    return any(marker in text for marker in markers)


def merge_one(base: pd.DataFrame, right: pd.DataFrame, on="gene_family") -> pd.DataFrame:
    if right.empty:
        return base
    return base.merge(right, on=on, how="left")


def choose_priority(row: pd.Series) -> str:
    score = num(row.get("scope_adjusted_story_score", row.get("accessory_story_score")))
    candidate = str(row.get("candidate_set", ""))
    if boolish(row.get("nif_or_nitrogenase_audit")):
        return "Core/diagnostic control"
    if boolish(row.get("housekeeping_lineage_flag")):
        return "De-emphasize: housekeeping/lineage"
    if candidate == "Model-Supported" and score >= 6.0:
        return "Tier A: story-leading accessory"
    if candidate == "Model-Supported" and score >= 4.5:
        return "Tier B: strong cross-evidence accessory"
    if candidate == "Highly Pure" and score >= 4.5:
        return "Tier HP: evidence-supported high-purity marker"
    if score >= 4.5:
        return "Tier C: cross-evidence diagnostic/support"
    if num(row.get("n_literature_active_phase_up_studies")) >= 2:
        return "Tier D: active-phase-up proteomics background"
    return "Tier E: lower current evidence"


def alignment_summary(row: pd.Series) -> str:
    pieces: list[str] = []
    candidate = str(row.get("candidate_set", ""))
    if candidate:
        pieces.append(candidate)
    n_lit_up = int(num(row.get("n_literature_active_phase_up_studies")))
    if n_lit_up:
        pieces.append(f"{n_lit_up} active-phase-up proteomics context{'s' if n_lit_up != 1 else ''}")
    else:
        n_lit = int(num(row.get("n_literature_nfix_response_studies")))
        if n_lit:
            pieces.append(f"{n_lit} broad N-fix-responsive proteomics context{'s' if n_lit != 1 else ''} (not scored)")
    if boolish(row.get("related_atlas_match")):
        pieces.append("related protein-family atlas concordance")
    if boolish(row.get("condensate_top20")):
        pieces.append("condensate-ranked")
    if boolish(row.get("primary_bridge_broad_ge1u_ge1f")):
        pieces.append("cross-morphotype bridge")
    return "; ".join(pieces)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for stale_name in [
        "external_literature_evidence_by_study.csv",
        "local_cyanothece_proteomics_rows.csv",
    ]:
        stale_path = OUT / stale_name
        if stale_path.exists():
            try:
                stale_path.unlink()
            except PermissionError:
                pass

    inventory = read_table(EVIDENCE / "morphotype_bridge_family_inventory.tsv", sep="\t")
    lit_summary = read_table(LIT / "family_summary.csv")
    lit_normalized = read_table(LIT / "normalized_proteins.csv")
    aryal = read_table(EVIDENCE / "aryal_current_shared_family_detail.csv")
    condensate = read_table(EVIDENCE / "condensate_ranking_family_overlay.tsv", sep="\t")

    base = inventory.copy()
    base = base[base["candidate_set"].isin(["Model-Supported", "Highly Pure"])].copy()

    lit_keep = [
        "gene_family",
        "n_literature_studies",
        "studies",
        "organisms",
        "n_study_family_rows",
        "n_unique_source_proteins_total",
        "source_protein_keys",
        "representative_protein_names",
    ]
    lit_summary = lit_summary[[c for c in lit_keep if c in lit_summary.columns]].copy()
    lit_summary = lit_summary.rename(
        columns={
            "n_literature_studies": "n_literature_mapped_studies",
            "studies": "mapped_literature_studies",
            "organisms": "mapped_literature_organisms",
            "n_study_family_rows": "n_mapped_study_family_rows",
            "n_unique_source_proteins_total": "n_mapped_source_proteins_total",
            "source_protein_keys": "mapped_source_protein_keys",
            "representative_protein_names": "mapped_representative_protein_names",
        }
    )

    if not lit_normalized.empty:
        classified_rows = []
        for _, row in lit_normalized.iterrows():
            row_dict = row.to_dict()
            row_dict.update(classify_literature_nfix_response(row))
            classified_rows.append(row_dict)
        lit_classified = pd.DataFrame(classified_rows)
    else:
        lit_classified = pd.DataFrame()

    if not lit_classified.empty:
        mapped = lit_classified[lit_classified["map_status"].eq("mapped_to_family")].copy()
        response_mapped = mapped[mapped["nfix_response_evidence"].apply(boolish)].copy()
        exploded_rows = []
        for _, row in response_mapped.iterrows():
            families = [x for x in str(row.get("mapped_gene_families", "")).split(";") if x]
            for gf in families:
                item = row.to_dict()
                item["gene_family"] = gf
                exploded_rows.append(item)
        lit_response_proteins = pd.DataFrame(exploded_rows)
    else:
        lit_response_proteins = pd.DataFrame()

    if not lit_response_proteins.empty:
        lit_response_summary = (
            lit_response_proteins.groupby("gene_family")
            .agg(
                n_literature_nfix_response_studies=("paper_short", lambda s: len(set(x for x in s if x))),
                nfix_response_studies=("paper_short", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                nfix_response_organisms=("organism", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                n_nfix_response_study_family_rows=("study_id", "size"),
                n_nfix_response_source_proteins_total=("source_protein_key", lambda s: len(set(x for x in s if x))),
                nfix_response_source_protein_keys=("source_protein_key", lambda s: ";".join(sorted(set(x for x in s if x))[:50])),
                nfix_response_representative_protein_names=("protein_name", lambda s: " | ".join(sorted(set(str(x) for x in s if str(x).strip()))[:8])),
                nfix_response_types=("nfix_response_type", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                nfix_response_directions=("nfix_response_direction", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                max_nfix_response_strength=("nfix_response_strength", "max"),
                nfix_response_criteria=("nfix_response_criterion", lambda s: " | ".join(sorted(set(x for x in s if x)))),
                nfix_response_metric_examples=("evidence_metric", lambda s: " | ".join([str(x) for x in s if str(x).strip()][:3])),
            )
            .reset_index()
        )
        context_rows = []
        for _, row in lit_response_summary.iterrows():
            studies = listify_semicolon(row.get("nfix_response_studies"))
            groups = sorted({STUDY_CONTEXTS.get(study, {}).get("context_group", "other") for study in studies})
            context_rows.append(
                {
                    "gene_family": row["gene_family"],
                    "nfix_response_context_groups": "; ".join(groups),
                    "nfix_response_context_count": len([g for g in groups if g and g != "other"]),
                    "nfix_response_has_unicellular": "unicellular" in groups,
                    "nfix_response_has_heterocyst": "heterocyst" in groups,
                    "nfix_response_has_trichodesmium": "nonheterocystous filament" in groups,
                }
            )
        lit_context = pd.DataFrame(context_rows)
        lit_response_summary = lit_response_summary.merge(lit_context, on="gene_family", how="left")
    else:
        lit_response_proteins = pd.DataFrame()
        lit_response_summary = pd.DataFrame(columns=["gene_family"])

    if not lit_response_proteins.empty and "nfix_active_phase_up" in lit_response_proteins.columns:
        active_up_proteins = lit_response_proteins[
            lit_response_proteins["nfix_active_phase_up"].apply(boolish)
        ].copy()
    else:
        active_up_proteins = pd.DataFrame()

    if not active_up_proteins.empty:
        active_up_summary = (
            active_up_proteins.groupby("gene_family")
            .agg(
                n_literature_active_phase_up_studies=("paper_short", lambda s: len(set(x for x in s if x))),
                active_phase_up_studies=("paper_short", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                active_phase_up_organisms=("organism", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                n_active_phase_up_study_family_rows=("study_id", "size"),
                n_active_phase_up_source_proteins_total=("source_protein_key", lambda s: len(set(x for x in s if x))),
                active_phase_up_source_protein_keys=("source_protein_key", lambda s: ";".join(sorted(set(x for x in s if x))[:50])),
                active_phase_up_representative_protein_names=("protein_name", lambda s: " | ".join(sorted(set(str(x) for x in s if str(x).strip()))[:8])),
                active_phase_up_alignment_classes=("nfix_phase_alignment_class", lambda s: "; ".join(sorted(set(x for x in s if x)))),
                active_phase_up_metric_examples=("nfix_active_phase_metric", lambda s: " | ".join([str(x) for x in s if str(x).strip()][:3])),
            )
            .reset_index()
        )
        active_context_rows = []
        for _, row in active_up_summary.iterrows():
            studies = listify_semicolon(row.get("active_phase_up_studies"))
            groups = sorted({STUDY_CONTEXTS.get(study, {}).get("context_group", "other") for study in studies})
            active_context_rows.append(
                {
                    "gene_family": row["gene_family"],
                    "active_phase_up_context_groups": "; ".join(groups),
                    "active_phase_up_context_count": len([g for g in groups if g and g != "other"]),
                    "active_phase_up_has_unicellular": "unicellular" in groups,
                    "active_phase_up_has_heterocyst": "heterocyst" in groups,
                    "active_phase_up_has_trichodesmium": "nonheterocystous filament" in groups,
                }
            )
        active_context = pd.DataFrame(active_context_rows)
        active_up_summary = active_up_summary.merge(active_context, on="gene_family", how="left")
    else:
        active_up_summary = pd.DataFrame(columns=["gene_family"])

    if not aryal.empty:
        aryal = aryal.copy()
        aryal["related_atlas_match"] = True
        aryal_cols = [
            "gene_family",
            "related_atlas_match",
            "aryal_gene",
            "aryal_product",
            "aryal_module_bin",
            "aryal_consensus_rank_pct_mean",
            "exact_product_match",
            "hgt_passenger_flag",
            "frac_assemblies_proximal",
            "fisher_p_vs_background",
        ]
        aryal = aryal[[c for c in aryal_cols if c in aryal.columns]].copy()

    cond_cols = [
        "gene_family",
        "condensate_mapped_rows",
        "condensate_best_score",
        "condensate_family_rank",
        "condensate_family_percentile",
        "condensate_label_names",
        "condensate_label_sources",
        "top_1pct",
        "top_5pct",
        "top_10pct",
        "top_20pct",
    ]
    condensate = condensate[[c for c in cond_cols if c in condensate.columns]].copy()
    condensate = condensate.rename(
        columns={
            "top_1pct": "condensate_top1",
            "top_5pct": "condensate_top5",
            "top_10pct": "condensate_top10",
            "top_20pct": "condensate_top20",
        }
    )

    composite = base
    composite = merge_one(base, lit_summary)
    for table in [lit_response_summary, active_up_summary, aryal, condensate]:
        composite = merge_one(composite, table)

    for col in ["n_literature_mapped_studies", "n_literature_nfix_response_studies", "n_literature_active_phase_up_studies"]:
        if col in composite.columns:
            composite[col] = composite[col].fillna(0).astype(int)
        else:
            composite[col] = 0
    for col in [
        "nfix_response_context_count",
        "active_phase_up_context_count",
        "max_nfix_response_strength",
        "aryal_consensus_rank_pct_mean",
        "condensate_family_percentile",
    ]:
        if col in composite.columns:
            composite[col] = pd.to_numeric(composite[col], errors="coerce")

    for col in [
        "nfix_response_has_unicellular",
        "nfix_response_has_heterocyst",
        "nfix_response_has_trichodesmium",
        "active_phase_up_has_unicellular",
        "active_phase_up_has_heterocyst",
        "active_phase_up_has_trichodesmium",
        "related_atlas_match",
        "exact_product_match",
        "hgt_passenger_flag",
        "condensate_top1",
        "condensate_top5",
        "condensate_top10",
        "condensate_top20",
    ]:
        if col in composite.columns:
            composite[col] = composite[col].apply(boolish)

    composite["story_role"] = composite.apply(
        lambda row: classify_story(
            str(row.get("product", "")),
            str(row.get("module_bin", "")),
            boolish(row.get("nif_or_nitrogenase_audit")),
        ),
        axis=1,
    )
    composite["housekeeping_lineage_flag"] = composite.apply(is_housekeeping_or_lineage, axis=1)

    composite["model_score"] = composite["candidate_set"].map({"Model-Supported": 2.0, "Highly Pure": 1.0}).fillna(0.0)
    composite["literature_proteomics_score"] = (
        composite["n_literature_active_phase_up_studies"].clip(upper=3) * 1.00
        + composite.get("active_phase_up_context_count", pd.Series(0, index=composite.index)).fillna(0).clip(upper=3).map({0: 0.0, 1: 0.0, 2: 0.5, 3: 1.0}).fillna(1.0)
    )
    composite["related_atlas_score"] = (
        composite.get("related_atlas_match", pd.Series(False, index=composite.index)).astype(float) * 0.50
        + composite.get("exact_product_match", pd.Series(False, index=composite.index)).astype(float) * 0.25
        + (composite.get("aryal_consensus_rank_pct_mean", pd.Series(float("nan"), index=composite.index)) <= 0.25).fillna(False).astype(float) * 0.25
    )
    composite["condensate_score"] = 0.0
    composite.loc[composite.get("condensate_top20", pd.Series(False, index=composite.index)), "condensate_score"] = 0.20
    composite.loc[composite.get("condensate_top10", pd.Series(False, index=composite.index)), "condensate_score"] = 0.30
    composite.loc[composite.get("condensate_top5", pd.Series(False, index=composite.index)), "condensate_score"] = 0.40
    composite.loc[composite.get("condensate_top1", pd.Series(False, index=composite.index)), "condensate_score"] = 0.50
    composite["morphotype_breadth_score"] = (
        composite.get("primary_bridge_broad_ge1u_ge1f", pd.Series(False, index=composite.index)).apply(boolish).astype(float) * 1.25
        + composite.get("strict_unicellular_breadth_ge3u_ge1f", pd.Series(False, index=composite.index)).apply(boolish).astype(float) * 1.50
    )
    component_cols = [
        "model_score",
        "literature_proteomics_score",
        "related_atlas_score",
        "condensate_score",
        "morphotype_breadth_score",
    ]
    composite["composite_evidence_score"] = composite[component_cols].sum(axis=1).round(3)
    composite["story_penalty"] = (
        composite["housekeeping_lineage_flag"].astype(float) * 1.00
        + composite["nif_or_nitrogenase_audit"].apply(boolish).astype(float) * 1.00
        + composite.get("hgt_passenger_flag", pd.Series(False, index=composite.index)).astype(float) * 0.75
    )
    composite["accessory_story_score"] = (composite["composite_evidence_score"] - composite["story_penalty"]).round(3)
    composite["related_atlas_applicability"] = composite["candidate_set"].map(
        {
            "Highly Pure": "not applicable by collaborator-atlas scope",
            "Model-Supported": "applicable",
        }
    ).fillna("applicable")
    composite["related_atlas_scope_adjustment_score"] = 0.0
    hp_mask = composite["candidate_set"].eq("Highly Pure")
    composite.loc[hp_mask, "related_atlas_scope_adjustment_score"] = (
        1.0 - composite.loc[hp_mask, "related_atlas_score"].fillna(0.0)
    ).clip(lower=0.0)
    composite["scope_adjusted_story_score"] = (
        composite["accessory_story_score"] + composite["related_atlas_scope_adjustment_score"]
    ).round(3)
    composite["priority_tier"] = composite.apply(choose_priority, axis=1)
    composite["evidence_alignment_summary"] = composite.apply(alignment_summary, axis=1)
    composite = composite.drop(columns=[c for c in composite.columns if "fox" in c.lower()], errors="ignore")

    preferred_cols = [
        "priority_tier",
        "gene_family",
        "candidate_set",
        "product",
        "module_bin",
        "story_role",
        "composite_evidence_score",
        "accessory_story_score",
        "scope_adjusted_story_score",
        "model_score",
        "literature_proteomics_score",
        "related_atlas_score",
        "related_atlas_scope_adjustment_score",
        "related_atlas_applicability",
        "condensate_score",
        "morphotype_breadth_score",
        "story_penalty",
        "n_literature_active_phase_up_studies",
        "active_phase_up_studies",
        "active_phase_up_context_groups",
        "active_phase_up_organisms",
        "n_active_phase_up_source_proteins_total",
        "active_phase_up_representative_protein_names",
        "active_phase_up_alignment_classes",
        "n_literature_nfix_response_studies",
        "nfix_response_studies",
        "nfix_response_context_groups",
        "nfix_response_organisms",
        "n_nfix_response_source_proteins_total",
        "nfix_response_types",
        "nfix_response_directions",
        "max_nfix_response_strength",
        "n_literature_mapped_studies",
        "mapped_literature_studies",
        "mapped_literature_organisms",
        "related_atlas_match",
        "aryal_product",
        "aryal_consensus_rank_pct_mean",
        "exact_product_match",
        "condensate_family_percentile",
        "condensate_best_score",
        "primary_bridge_broad_ge1u_ge1f",
        "strict_unicellular_breadth_ge3u_ge1f",
        "nif_or_nitrogenase_audit",
        "housekeeping_lineage_flag",
        "hgt_passenger_flag",
        "evidence_alignment_summary",
    ]
    ordered_cols = [col for col in preferred_cols if col in composite.columns] + [
        col for col in composite.columns if col not in preferred_cols
    ]
    composite = composite[ordered_cols].sort_values(
        ["scope_adjusted_story_score", "accessory_story_score", "composite_evidence_score", "n_literature_active_phase_up_studies"],
        ascending=[False, False, False, False],
    )

    story_top = composite[
        composite["priority_tier"].isin(
            [
                "Tier A: story-leading accessory",
                "Tier B: strong cross-evidence accessory",
                "Tier HP: evidence-supported high-purity marker",
                "Tier C: cross-evidence diagnostic/support",
            ]
        )
    ].copy()
    tier_order = {
        "Tier A: story-leading accessory": 1,
        "Tier B: strong cross-evidence accessory": 2,
        "Tier HP: evidence-supported high-purity marker": 3,
        "Tier C: cross-evidence diagnostic/support": 4,
    }
    story_top["_tier_order"] = story_top["priority_tier"].map(tier_order).fillna(9)
    story_top = story_top.sort_values(
        ["_tier_order", "scope_adjusted_story_score", "accessory_story_score", "composite_evidence_score"],
        ascending=[True, False, False, False],
    ).drop(columns=["_tier_order"])
    hp_markers = composite[
        composite["priority_tier"].eq("Tier HP: evidence-supported high-purity marker")
    ].copy()
    hp_markers = hp_markers.sort_values(
        ["scope_adjusted_story_score", "accessory_story_score", "composite_evidence_score", "n_literature_active_phase_up_studies"],
        ascending=[False, False, False, False],
    )

    component_definitions = pd.DataFrame(
        [
            ["model_score", "2.0 for Model-Supported; 1.0 for Highly Pure", "Separates classifier-supported accessory candidates from near-diagnostic carrier-purity features."],
            ["literature_proteomics_score", "1.0 per external active-phase-up proteomics study, capped at 3; +0.5 for two active-up context groups; +1.0 for three", "Rewards only proteins that increase in the organism-appropriate N-fix-active phase/cell. Broad response, active-phase-down, transition, and non-target phase rows are retained as metadata but do not score."],
            ["related_atlas_score", "+0.5 related protein-family atlas match; +0.25 exact product match; +0.25 top-quartile atlas ranking", "Caps collaborator/related-atlas convergence at 1.0 so it supports but does not dominate the biological filter."],
            ["condensate_score", "+0.2 top 20%; +0.3 top 10%; +0.4 top 5%; +0.5 top 1% condensate-driver rank", "Keeps the condensate-associated layer as a minor prioritization signal."],
            ["morphotype_breadth_score", "+1.25 broad unicellular/filament bridge; +1.5 strict unicellular breadth", "Makes cross-morphotype conservation a primary biological support axis."],
            ["story_penalty", "-1.0 housekeeping/lineage; -1.0 core nif/nitrogenase; -0.75 HGT-passenger flag", "Does not delete these features, but prevents them from leading the accessory-biology narrative."],
            ["accessory_story_score", "composite_evidence_score minus story_penalty", "Raw biological-story score before scope adjustment."],
            ["related_atlas_scope_adjustment_score", "For Highly Pure only: 1.0 minus observed related_atlas_score", "Neutralizes missing collaborator-atlas evidence when the related atlas was scoped around Model-Supported families."],
            ["scope_adjusted_story_score", "accessory_story_score plus related_atlas_scope_adjustment_score", "Primary tiering/sort score for top-family tables."],
        ],
        columns=["component", "rule", "interpretation"],
    )

    module_rows = []
    for role, group in composite.groupby("story_role", dropna=False):
        top_examples = group.head(5)
        module_rows.append(
            {
                "story_role": role,
                "n_families": len(group),
                "n_model_supported": int((group["candidate_set"] == "Model-Supported").sum()),
                "n_highly_pure": int((group["candidate_set"] == "Highly Pure").sum()),
                "n_with_external_nfix_response_proteomics": int((group["n_literature_nfix_response_studies"] > 0).sum()),
                "n_with_external_active_phase_up_proteomics": int((group["n_literature_active_phase_up_studies"] > 0).sum()),
                "n_with_related_atlas_match": int(group.get("related_atlas_match", pd.Series(False, index=group.index)).astype(bool).sum()),
                "mean_accessory_story_score": round(group["accessory_story_score"].mean(), 3),
                "mean_scope_adjusted_story_score": round(group["scope_adjusted_story_score"].mean(), 3),
                "top_family_examples": "; ".join(
                    f"{r.gene_family} ({str(r.product)[:60]})"
                    for r in top_examples.itertuples(index=False)
                ),
            }
        )
    module_summary = pd.DataFrame(module_rows).sort_values(
        ["mean_scope_adjusted_story_score", "n_with_external_active_phase_up_proteomics", "n_families"],
        ascending=[False, False, False],
    )

    lit_response_aug = lit_response_proteins.copy()
    if not lit_response_aug.empty and "paper_short" in lit_response_aug.columns:
        for key in ["context", "context_group", "physiology", "organism_group"]:
            lit_response_aug[key] = lit_response_aug["paper_short"].map(
                {study: meta[key] for study, meta in STUDY_CONTEXTS.items()}
            )
    lit_response_aug = lit_response_aug.merge(
        composite[["gene_family", "candidate_set", "product", "module_bin", "story_role", "priority_tier", "accessory_story_score"]],
        on="gene_family",
        how="left",
    )

    classified_exploded_rows = []
    if not lit_classified.empty:
        for _, row in lit_classified.iterrows():
            families = [x for x in str(row.get("mapped_gene_families", "")).split(";") if x]
            if not families:
                item = row.to_dict()
                item["gene_family"] = ""
                classified_exploded_rows.append(item)
            for gf in families:
                item = row.to_dict()
                item["gene_family"] = gf
                classified_exploded_rows.append(item)
    lit_classified_out = pd.DataFrame(classified_exploded_rows)
    if not lit_classified_out.empty and "gene_family" in lit_classified_out.columns:
        lit_classified_out = lit_classified_out.merge(
            composite[["gene_family", "candidate_set", "product", "module_bin", "story_role", "priority_tier", "accessory_story_score"]],
            on="gene_family",
            how="left",
        )

    summary = pd.DataFrame(
        [
            ["all_candidate_families", len(composite)],
            ["model_supported_families", int((composite["candidate_set"] == "Model-Supported").sum())],
            ["highly_pure_families", int((composite["candidate_set"] == "Highly Pure").sum())],
            ["families_with_external_literature_mapping", int((composite["n_literature_mapped_studies"] > 0).sum())],
            ["families_with_external_nfix_response_proteomics", int((composite["n_literature_nfix_response_studies"] > 0).sum())],
            ["families_with_external_active_phase_up_proteomics", int((composite["n_literature_active_phase_up_studies"] > 0).sum())],
            ["families_with_two_or_more_external_active_phase_up_contexts", int((composite["n_literature_active_phase_up_studies"] >= 2).sum())],
            ["families_with_related_protein_atlas_match", int(composite.get("related_atlas_match", pd.Series(False, index=composite.index)).astype(bool).sum())],
            ["tier_a_story_leading_accessory", int((composite["priority_tier"] == "Tier A: story-leading accessory").sum())],
            ["tier_b_strong_accessory", int((composite["priority_tier"] == "Tier B: strong cross-evidence accessory").sum())],
            ["tier_hp_evidence_supported_high_purity_marker", int((composite["priority_tier"] == "Tier HP: evidence-supported high-purity marker").sum())],
            ["core_or_diagnostic_controls", int((composite["priority_tier"] == "Core/diagnostic control").sum())],
            ["housekeeping_lineage_deemphasized", int((composite["priority_tier"] == "De-emphasize: housekeeping/lineage").sum())],
        ],
        columns=["metric", "value"],
    )

    outputs = {
        "family_composite_all.csv": composite,
        "top_story_families.csv": story_top,
        "high_purity_marker_families.csv": hp_markers,
        "module_story_summary.csv": module_summary,
        "component_definitions.csv": component_definitions,
        "external_literature_nfix_response_evidence.csv": lit_response_aug,
        "external_literature_protein_rows_nfix_classified.csv": lit_classified_out,
        "summary_metrics.csv": summary,
    }
    for name, frame in outputs.items():
        frame.to_csv(OUT / name, index=False, quoting=csv.QUOTE_MINIMAL)

    readme = """Composite proteomics-prioritized family evidence table

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
"""
    (OUT / "README.txt").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
