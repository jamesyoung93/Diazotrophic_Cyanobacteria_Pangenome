from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT_WORKBOOK = ROOT / "outputs" / "composite_proteomics_priority_morphotype_weighted_hp" / "proteomics_composite_family_evidence.xlsx"
OUTDIR = ROOT / "evidence_tables" / "literature_deep_dive" / "phase_alignment_audit"


def as_float(value: object, default: float = math.nan) -> float:
    try:
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none"}:
            return default
        return float(text)
    except Exception:
        return default


def metric_value(metric_text: object, label: str, default: float = math.nan) -> float:
    text = "" if metric_text is None else str(metric_text)
    match = re.search(re.escape(label) + r"\s*=\s*([^;|]+)", text)
    if not match:
        return default
    return as_float(match.group(1), default=default)


def metric_hour(metric_text: object, label: str, default: float = math.nan) -> float:
    text = "" if metric_text is None else str(metric_text)
    match = re.search(re.escape(label) + r"\s*=\s*[^;|@]+@([0-9.]+)h", text)
    if not match:
        return default
    return as_float(match.group(1), default=default)


def day_window(hour: float) -> bool:
    return not math.isnan(hour) and 6 <= hour <= 15


def transition_window(hour: float) -> bool:
    return not math.isnan(hour) and hour in {4, 18}


def classify_row(row: pd.Series) -> dict[str, object]:
    study = str(row.get("paper_short", ""))
    metric = row.get("evidence_metric", "")
    out = {
        "intended_nfix_phase": "",
        "phase_alignment_class": "not_classified",
        "active_phase_response": False,
        "active_phase_up": False,
        "active_phase_down": False,
        "transition_or_ambiguous": False,
        "phase_primary_metric": "",
        "phase_secondary_metric": "",
        "phase_reason": "",
    }

    if study == "Panda et al. 2025":
        nitrate_p = metric_value(metric, "Two-way ANOVA p value nitrate")
        interaction_p = metric_value(metric, "Two-way ANOVA p value Interaction")
        dark_logfc = metric_value(metric, "D-_D+(logFC)")
        light_logfc = metric_value(metric, "L-_L+(logFC)")
        sig = (not math.isnan(nitrate_p) and nitrate_p <= 0.05) or (
            not math.isnan(interaction_p) and interaction_p <= 0.05
        )
        dark_hit = sig and not math.isnan(dark_logfc) and abs(dark_logfc) >= 0.25
        light_hit = sig and not math.isnan(light_logfc) and abs(light_logfc) >= 0.25
        out["intended_nfix_phase"] = "dark nitrate-free Crocosphaera phase; D- versus D+ is the primary active-phase contrast"
        out["phase_primary_metric"] = f"D-_D+(logFC)={dark_logfc:.4g}; p_nitrate={nitrate_p:.4g}; p_interaction={interaction_p:.4g}"
        out["phase_secondary_metric"] = f"L-_L+(logFC)={light_logfc:.4g}"
        if dark_hit:
            out["active_phase_response"] = True
            out["active_phase_up"] = dark_logfc > 0
            out["active_phase_down"] = dark_logfc < 0
            if light_hit and dark_logfc * light_logfc < 0:
                out["phase_alignment_class"] = "aligned_dark_active_contrast_with_light_dark_interaction"
                out["phase_reason"] = "The dark N-fixing contrast passes threshold and the light contrast has the opposite sign."
            elif light_hit:
                out["phase_alignment_class"] = "aligned_dark_active_contrast_broad_minusN"
                out["phase_reason"] = "The dark N-fixing contrast passes threshold; the light contrast also changes in the same direction."
            else:
                out["phase_alignment_class"] = "aligned_dark_active_contrast_only"
                out["phase_reason"] = "The row is supported specifically by the dark N-fixing contrast."
        elif light_hit:
            out["phase_alignment_class"] = "not_active_phase_light_only"
            out["phase_reason"] = "The current filter keeps this row because the light minus-nitrate contrast passes, but the dark active-phase contrast does not."
        else:
            out["phase_alignment_class"] = "ambiguous_after_reparse"
            out["transition_or_ambiguous"] = True
            out["phase_reason"] = "The row was retained by the workbook but did not re-pass the parsed active-phase threshold."

    elif study == "Held et al. 2022":
        max_hour = metric_hour(metric, "max")
        min_hour = metric_hour(metric, "min")
        out["intended_nfix_phase"] = "daylight Trichodesmium N2-fixing interval; 6-15 h post-dawn treated as active phase and 4/18 h as transitions"
        out["phase_primary_metric"] = f"max_hour={max_hour:g}; min_hour={min_hour:g}"
        if day_window(max_hour):
            out["phase_alignment_class"] = "aligned_daytime_active_phase_high"
            out["active_phase_response"] = True
            out["active_phase_up"] = True
            out["phase_reason"] = "The diel maximum falls in the daylight interval used as the active Trichodesmium N-fixing window."
        elif day_window(min_hour):
            out["phase_alignment_class"] = "aligned_daytime_active_phase_low"
            out["active_phase_response"] = True
            out["active_phase_down"] = True
            out["phase_reason"] = "The diel minimum falls in the daylight interval, so this is an active-phase depletion rather than an active-phase increase."
        elif transition_window(max_hour) or transition_window(min_hour):
            out["phase_alignment_class"] = "transition_window_dynamic"
            out["transition_or_ambiguous"] = True
            out["phase_reason"] = "The strongest extrema are near the light/dark transition; keep only as transitional or preparatory evidence."
        else:
            out["phase_alignment_class"] = "not_active_phase_night_only_dynamic"
            out["phase_reason"] = "The retained diel dynamic is centered outside the daylight active N-fix interval."

    elif study == "Welkie et al. 2014":
        labels = ["CT_D0", "CT_D3", "CT_L0", "CT_L3"]
        vals = {label: metric_value(metric, label) for label in labels}
        finite = {k: v for k, v in vals.items() if not math.isnan(v)}
        out["intended_nfix_phase"] = "late-light/early-dark CT_D0 and CT_D3 window in a unicellular diazotroph diel program"
        out["phase_primary_metric"] = "; ".join(f"{k}={v:.4g}" for k, v in finite.items())
        if finite:
            max_label = max(finite, key=finite.get)
            min_label = min(finite, key=finite.get)
            out["phase_secondary_metric"] = f"max={max_label}; min={min_label}"
            if max_label.startswith("CT_D"):
                out["phase_alignment_class"] = "aligned_CT_D_active_window_high"
                out["active_phase_response"] = True
                out["active_phase_up"] = True
                out["phase_reason"] = "The highest abundance is in the CT_D0/CT_D3 active-window proxy."
            elif min_label.startswith("CT_D"):
                out["phase_alignment_class"] = "aligned_CT_D_active_window_low"
                out["active_phase_response"] = True
                out["active_phase_down"] = True
                out["phase_reason"] = "The lowest abundance is in the CT_D0/CT_D3 active-window proxy, so this is active-window depletion rather than an increase."
            else:
                out["phase_alignment_class"] = "not_CT_D_active_window_extreme"
                out["phase_reason"] = "The dynamic range does not place CT_D0 or CT_D3 at either abundance extreme."
        else:
            out["phase_alignment_class"] = "ambiguous_after_reparse"
            out["transition_or_ambiguous"] = True
            out["phase_reason"] = "No parsed CT_D/CT_L abundance values were available."

    elif study == "Sandh et al. 2014":
        direction = str(row.get("nfix_response_direction", ""))
        out["intended_nfix_phase"] = "heterocyst cell type in N2-fixing Nostoc filaments"
        out["phase_primary_metric"] = str(row.get("evidence_metric", ""))
        if direction == "heterocyst-enriched":
            out["phase_alignment_class"] = "aligned_heterocyst_active_cell_high"
            out["active_phase_response"] = True
            out["active_phase_up"] = True
            out["phase_reason"] = "The protein is enriched in isolated heterocysts relative to N2-fixing filaments."
        elif direction == "filament/vegetative-enriched":
            out["phase_alignment_class"] = "aligned_heterocyst_active_cell_low_or_supporting_filament"
            out["active_phase_response"] = True
            out["active_phase_down"] = True
            out["phase_reason"] = "The protein is lower in heterocysts than the N2-fixing filament mixture; use as compartment contrast/supporting-vegetative evidence, not as heterocyst-up evidence."
        else:
            out["phase_alignment_class"] = "ambiguous_after_reparse"
            out["transition_or_ambiguous"] = True
            out["phase_reason"] = "The heterocyst/filament direction was not parsed."

    return out


def semicolon_items(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def build_report(response: pd.DataFrame, top_join: pd.DataFrame, summary: pd.DataFrame, narrative_summary: pd.DataFrame) -> str:
    strict_up = int(response["active_phase_up"].sum())
    active = int(response["active_phase_response"].sum())
    total = len(response)
    not_active = total - active - int(response["transition_or_ambiguous"].sum())

    lines = [
        "# Nitrogen-fixation phase alignment audit",
        "",
        f"Input workbook: `{INPUT_WORKBOOK.relative_to(ROOT)}`",
        "",
        "## Bottom line",
        "",
        (
            f"The current literature-proteomics layer contains {total} retained response rows. "
            f"{active} rows are aligned to the intended active nitrogen-fixation phase or cell-type contrast, "
            f"and {strict_up} of those are higher in the active phase/cell. "
            f"{not_active} rows are not active-phase-specific under the stricter audit, with the remainder flagged as transition/ambiguous."
        ),
        "",
        "This means the existing overlay is usable as a broad nitrogen-fixation-condition response layer, but the manuscript should not describe every retained row as an active-phase increase. For the strongest reviewer-proof framing, use two labels: active-phase up and active-phase differential response.",
        "",
        "## Study-level interpretation",
        "",
    ]

    for study in sorted(response["paper_short"].dropna().unique()):
        sub = response[response["paper_short"] == study]
        counts = sub["phase_alignment_class"].value_counts().to_dict()
        lines.append(f"### {study}")
        for cls, count in counts.items():
            lines.append(f"- {cls}: {count}")
        if study == "Panda et al. 2025":
            lines.append("- Recommendation: make D- versus D+ the primary Crocosphaera active-phase contrast. Keep light-only rows only as light/dark interaction or non-active-phase nitrate-response evidence.")
        elif study == "Held et al. 2022":
            lines.append("- Recommendation: for a strict Trichodesmium active-fixation layer, distinguish daytime-high proteins from daytime-low proteins and downgrade night-only dynamics.")
        elif study == "Welkie et al. 2014":
            lines.append("- Recommendation: keep rows where CT_D0 or CT_D3 is an abundance extreme as late-light/early-dark active-window differential evidence; only the CT_D-high subset should be used as 'up in the N-fixing window' examples.")
        elif study == "Sandh et al. 2014":
            lines.append("- Recommendation: call heterocyst-enriched rows direct active-cell evidence; call filament/vegetative-enriched rows supporting compartment evidence.")
        lines.append("")

    lines.extend(
        [
            "## Top-family effect",
            "",
            f"Top-family rows with any broad literature response: {int((top_join['n_literature_nfix_response_studies'].fillna(0) > 0).sum())}",
            f"Top-family rows with active-phase differential evidence after this audit: {int((top_join['n_phase_aligned_response_studies'].fillna(0) > 0).sum())}",
            f"Top-family rows with active-phase-up evidence after this audit: {int((top_join['n_phase_up_response_studies'].fillna(0) > 0).sum())}",
            "",
            "The difference between these numbers is the exact reviewer-risk zone. It is not a mapping failure; it is a wording and scoring distinction.",
            "",
            "## Files written",
            "",
            "- `phase_alignment_response_rows.csv`: row-level audit of every retained literature-proteomics response.",
            "- `phase_alignment_study_summary.csv`: study-level counts by phase-alignment class.",
            "- `phase_alignment_top_family_summary.csv`: family-level phase-aligned study counts joined to the top-family table.",
            "- `phase_alignment_narrative_summary.csv`: narrative-level support after phase audit.",
            "- `nfix_phase_alignment_audit.xlsx`: sortable workbook with the same audit tables.",
            "",
        ]
    )

    if not narrative_summary.empty:
        lines.extend(["## Narrative-level signal", ""])
        preview = narrative_summary.sort_values(
            ["n_phase_aligned_families", "n_phase_up_families", "n_families"],
            ascending=False,
        ).head(12)
        for _, row in preview.iterrows():
            lines.append(
                f"- {row['narrative_title']}: {int(row['n_phase_aligned_families'])} phase-aligned families, "
                f"{int(row['n_phase_up_families'])} active-phase-up families, {int(row['n_families'])} total narrative families."
            )
        lines.append("")

    return "\n".join(lines)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    response = pd.read_excel(INPUT_WORKBOOK, sheet_name="LiteratureNfixResponse")
    top = pd.read_excel(INPUT_WORKBOOK, sheet_name="TopFamilies")

    phase = response.apply(classify_row, axis=1, result_type="expand")
    response_phase = pd.concat([response, phase], axis=1)

    study_summary = (
        response_phase.groupby(["study_id", "paper_short", "phase_alignment_class"], dropna=False)
        .agg(
            n_rows=("gene_family", "size"),
            n_families=("gene_family", "nunique"),
            n_active_phase_response=("active_phase_response", "sum"),
            n_active_phase_up=("active_phase_up", "sum"),
            n_active_phase_down=("active_phase_down", "sum"),
            n_transition_or_ambiguous=("transition_or_ambiguous", "sum"),
        )
        .reset_index()
    )

    family_phase = (
        response_phase.groupby("gene_family", dropna=False)
        .agg(
            n_phase_audited_rows=("gene_family", "size"),
            n_phase_aligned_rows=("active_phase_response", "sum"),
            n_phase_up_rows=("active_phase_up", "sum"),
            n_phase_down_rows=("active_phase_down", "sum"),
            n_phase_transition_or_ambiguous_rows=("transition_or_ambiguous", "sum"),
            n_phase_aligned_response_studies=("paper_short", lambda s: len(set(response_phase.loc[s.index[response_phase.loc[s.index, "active_phase_response"]], "paper_short"]))),
            phase_aligned_studies=("paper_short", lambda s: "; ".join(sorted(set(response_phase.loc[s.index[response_phase.loc[s.index, "active_phase_response"]], "paper_short"])))),
            n_phase_up_response_studies=("paper_short", lambda s: len(set(response_phase.loc[s.index[response_phase.loc[s.index, "active_phase_up"]], "paper_short"]))),
            phase_up_studies=("paper_short", lambda s: "; ".join(sorted(set(response_phase.loc[s.index[response_phase.loc[s.index, "active_phase_up"]], "paper_short"])))),
            phase_alignment_classes=("phase_alignment_class", lambda s: "; ".join(sorted(set(str(x) for x in s if str(x).strip())))),
            phase_alignment_reasons=("phase_reason", lambda s: " | ".join(list(dict.fromkeys(str(x) for x in s if str(x).strip()))[:4])),
        )
        .reset_index()
    )

    top_join = top.merge(family_phase, on="gene_family", how="left")
    fill_zero = [
        "n_phase_audited_rows",
        "n_phase_aligned_rows",
        "n_phase_up_rows",
        "n_phase_down_rows",
        "n_phase_transition_or_ambiguous_rows",
        "n_phase_aligned_response_studies",
        "n_phase_up_response_studies",
    ]
    for col in fill_zero:
        if col in top_join.columns:
            top_join[col] = top_join[col].fillna(0).astype(int)
    for col in ["phase_aligned_studies", "phase_up_studies", "phase_alignment_classes", "phase_alignment_reasons"]:
        if col in top_join.columns:
            top_join[col] = top_join[col].fillna("")

    narrative_path = ROOT / "evidence_tables" / "literature_deep_dive" / "narrative_family_memberships.csv"
    if narrative_path.exists():
        narrative = pd.read_csv(narrative_path)
        narrative_join = narrative.merge(family_phase, on="gene_family", how="left")
        for col in fill_zero:
            if col in narrative_join.columns:
                narrative_join[col] = narrative_join[col].fillna(0).astype(int)
        narrative_summary = (
            narrative_join.groupby(["narrative_id", "narrative_title"], dropna=False)
            .agg(
                n_families=("gene_family", "nunique"),
                n_phase_aligned_families=("n_phase_aligned_response_studies", lambda s: int((s > 0).sum())),
                n_phase_up_families=("n_phase_up_response_studies", lambda s: int((s > 0).sum())),
                n_phase_aligned_rows=("n_phase_aligned_rows", "sum"),
                n_phase_up_rows=("n_phase_up_rows", "sum"),
                phase_aligned_studies=("phase_aligned_studies", lambda s: "; ".join(sorted(set(x for value in s for x in semicolon_items(value))))),
            )
            .reset_index()
        )
    else:
        narrative_summary = pd.DataFrame()

    response_phase.to_csv(OUTDIR / "phase_alignment_response_rows.csv", index=False)
    study_summary.to_csv(OUTDIR / "phase_alignment_study_summary.csv", index=False)
    top_join.to_csv(OUTDIR / "phase_alignment_top_family_summary.csv", index=False)
    if not narrative_summary.empty:
        narrative_summary.to_csv(OUTDIR / "phase_alignment_narrative_summary.csv", index=False)

    workbook_path = OUTDIR / "nfix_phase_alignment_audit.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        study_summary.to_excel(writer, sheet_name="StudySummary", index=False)
        top_join.to_excel(writer, sheet_name="TopFamiliesPhaseAudit", index=False)
        response_phase.to_excel(writer, sheet_name="ResponseRows", index=False)
        if not narrative_summary.empty:
            narrative_summary.to_excel(writer, sheet_name="NarrativeSummary", index=False)

        for worksheet in writer.book.worksheets:
            worksheet.freeze_panes = "A2"
            for col_cells in worksheet.columns:
                header = str(col_cells[0].value or "")
                width = min(max(len(header) + 2, 12), 48)
                worksheet.column_dimensions[col_cells[0].column_letter].width = width

    report = build_report(response_phase, top_join, study_summary, narrative_summary)
    (OUTDIR / "nfix_phase_alignment_report.md").write_text(report, encoding="utf-8")

    print(report)


if __name__ == "__main__":
    main()
