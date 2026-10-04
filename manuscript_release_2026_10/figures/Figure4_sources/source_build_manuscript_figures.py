"""Build manuscript figures from inspected evidence tables.

The figures are deterministic and data-driven. Every value is read from the
evidence tables in ``results/tables/`` so the display items cannot drift from
the committed numbers. Styling targets Genome Biology display quality: a single
sans-serif type family, a colorblind-safe Okabe-Ito palette, editable vector
text in the PDF copies, and value labels where they aid reading. The factual
content and table sources are unchanged from the first-pass figures.
"""

from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"

# Okabe-Ito colorblind-safe qualitative palette.
OKABE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#999999",
}
TWO_GROUP = (OKABE["blue"], OKABE["orange"])

# Module color map (keys must match the module_bin strings in the tables).
PALETTE = {
    "Nitrogenase / nif machinery": OKABE["black"],
    "Carbon metabolism & NAD(P)H supply": OKABE["blue"],
    "Metallocluster & Fe-S biogenesis": OKABE["orange"],
    "O2 protection & redox homeostasis": OKABE["green"],
    "Respiration & bioenergetics": OKABE["vermillion"],
    "Regulation": OKABE["purple"],
    "Transport": OKABE["sky"],
    "Stress / membrane / uncharacterized": OKABE["yellow"],
    "Housekeeping / lineage marker": OKABE["grey"],
    "Other / unassigned": "#cfcfcf",
}

GRID = "#e5e7eb"
INK = "#222222"
SUBTLE = "#4b5563"


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 11,
            "axes.titlesize": 12.5,
            "axes.titleweight": "bold",
            "axes.titlepad": 10,
            "axes.labelsize": 11.5,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": SUBTLE,
            "axes.linewidth": 0.8,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": SUBTLE,
            "ytick.color": SUBTLE,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 300,
            # Embed editable TrueType text in vector output for journal editing.
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLE_DIR / name, sep="\t")


def save(fig: plt.Figure, name: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / name, dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / name.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)


def wrap_label(label: str, width: int = 28) -> str:
    return "\n".join(fill(str(label), width=width).splitlines())


def label_barh(ax, bars, values, fmt: str, pad: float) -> None:
    for bar, value in zip(bars, values):
        width = bar.get_width()
        if width <= 0:
            continue
        ax.text(
            width + pad,
            bar.get_y() + bar.get_height() / 2,
            fmt.format(value),
            va="center",
            ha="left",
            fontsize=8,
            color=SUBTLE,
        )


# Faithful short forms for long product strings, used only for figure axis
# labels. The full product names remain in Table 3 and the evidence tables.
PRODUCT_SHORT = {
    "CobW family GTP-binding protein": "CobW GTP-binding protein",
    "type I glyceraldehyde-3-phosphate dehydrogenase": "GAPDH (type I)",
    "50S ribosomal protein L5": "ribosomal protein L5",
    "ArsJ-associated glyceraldehyde-3-phosphate dehydrogenase": "ArsJ-associated GAPDH",
    "magnesium-protoporphyrin IX monomethyl ester (oxidative) cyclase": "Mg-protoporphyrin cyclase",
    "STAS domain-containing protein": "STAS domain protein",
    "glutamate synthase small subunit": "glutamate synthase (small)",
    "fumarate reductase/succinate dehydrogenase flavoprotein subunit": "fumarate reductase/SDH",
    "manganese catalase family protein": "manganese catalase",
    "methionine adenosyltransferase": "Met adenosyltransferase",
    "amino acid ABC transporter ATP-binding protein": "amino acid ABC transporter",
    "sulfate adenylyltransferase": "sulfate adenylyltransferase",
    "ATP-dependent Clp endopeptidase proteolytic subunit ClpP": "Clp protease ClpP",
    "response regulator transcription factor": "response regulator TF",
    "3-deoxy-7-phosphoheptulonate synthase": "DAHP synthase",
    "photosystem I biogenesis protein BtpA": "PSI biogenesis BtpA",
    "uroporphyrinogen decarboxylase": "uroporphyrinogen decarboxylase",
}


def short_product(product: str, maxlen: int = 26) -> str:
    p = str(product)
    if p in PRODUCT_SHORT:
        return PRODUCT_SHORT[p]
    if len(p) <= maxlen:
        return p
    cut = p[:maxlen]
    if " " in cut:
        cut = cut[: cut.rfind(" ")]
    return cut


def figure1_pipeline() -> None:
    steps = [
        ("489 cyanobacterial\nassemblies screened", "Complete RefSeq\nGCF filter"),
        ("nifH, nifD, nifK\nHMM marker calls", "Diazotrophy\npotential labels"),
        ("Predicted proteins\nclustered (MMseqs2)", "426 x 2286\nbinary matrix"),
        ("Genus-blocked\nsupervised models", "LR, RF, GB,\nXGBoost"),
        ("Two candidate\ninventories", "Model-Supported\nHighly Pure"),
        ("Interpretation\nlayers", "Modules, synteny,\nprotein vs gene"),
    ]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # Two rows of three boxes in a left-to-right then right-to-left snake.
    col_x = [0.18, 0.5, 0.82]
    top_y, bot_y = 0.74, 0.34
    pos = [
        (col_x[0], top_y),
        (col_x[1], top_y),
        (col_x[2], top_y),
        (col_x[2], bot_y),
        (col_x[1], bot_y),
        (col_x[0], bot_y),
    ]
    hw, hh = 0.145, 0.135
    for (title, detail), (x, y) in zip(steps, pos):
        box = FancyBboxPatch(
            (x - hw, y - hh),
            2 * hw,
            2 * hh,
            boxstyle="round,pad=0.008,rounding_size=0.02",
            facecolor="#eef2f7",
            edgecolor="#334155",
            linewidth=1.1,
        )
        ax.add_patch(box)
        ax.add_patch(
            plt.Rectangle((x - hw, y + hh - 0.014), 2 * hw, 0.014, facecolor="#334155", edgecolor="none")
        )
        ax.text(x, y + 0.052, title, ha="center", va="center", fontsize=9.6, fontweight="bold", color="#0f172a")
        ax.text(x, y - 0.062, detail, ha="center", va="center", fontsize=8.9, color="#475569")

    def arrow(p0, p1):
        ax.add_patch(
            FancyArrowPatch(
                p0, p1, arrowstyle="-|>", mutation_scale=15, linewidth=1.4, color="#64748b"
            )
        )

    arrow((col_x[0] + hw + 0.004, top_y), (col_x[1] - hw - 0.004, top_y))
    arrow((col_x[1] + hw + 0.004, top_y), (col_x[2] - hw - 0.004, top_y))
    arrow((col_x[2], top_y - hh - 0.004), (col_x[2], bot_y + hh + 0.004))
    arrow((col_x[2] - hw - 0.004, bot_y), (col_x[1] + hw + 0.004, bot_y))
    arrow((col_x[1] - hw - 0.004, bot_y), (col_x[0] + hw + 0.004, bot_y))

    ax.text(
        0.5,
        0.085,
        "Outputs are ranked protein families, high-purity context families,",
        ha="center",
        fontsize=9.3,
        color=SUBTLE,
    )
    ax.text(
        0.5,
        0.03,
        "conserved neighborhoods, and engineering-prioritization tables",
        ha="center",
        fontsize=9.3,
        color=SUBTLE,
    )
    save(fig, "figure1_pipeline_study_design.png")


def figure2_model_performance() -> None:
    key = read_tsv("manuscript_key_results.tsv")
    rows = key[key["claim"].str.contains("cross-validated ROC AUC", regex=False)].copy()
    parsed = rows["value"].str.extract(r"(?P<mean>[0-9.]+)\s*\+/-\s*(?P<sd>[0-9.]+)").astype(float)
    rows = pd.concat([rows, parsed], axis=1)
    rows["model"] = rows["claim"].str.replace(" cross-validated ROC AUC", "", regex=False)
    order = ["Random Forest", "Gradient Boosting", "Logistic Regression", "XGBoost"]
    rows["model"] = pd.Categorical(rows["model"], order, ordered=True)
    rows = rows.sort_values("model")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    x = range(len(rows))
    ax.bar(
        x,
        rows["mean"],
        yerr=rows["sd"],
        capsize=5,
        color=OKABE["blue"],
        edgecolor=INK,
        linewidth=0.7,
        error_kw={"elinewidth": 1.1, "ecolor": SUBTLE},
        width=0.64,
    )
    for xi, (mean, sd) in enumerate(zip(rows["mean"], rows["sd"])):
        ax.text(xi, mean + sd + 0.02, f"{mean:.3f}", ha="center", va="bottom", fontsize=10, color=INK)
    ax.axhline(0.5, ls="--", lw=1.0, color=OKABE["grey"])
    ax.text(-0.45, 0.515, "chance (0.5)", ha="left", va="bottom", fontsize=9.5, color=SUBTLE)
    ax.axhline(1.0, ls=":", lw=0.9, color=OKABE["grey"])
    ax.text(-0.45, 1.005, "ROC AUC maximum (1.0)", ha="left", va="bottom", fontsize=9, color=SUBTLE)
    ax.set_ylim(0, 1.14)
    ax.set_xlim(-0.6, len(rows) - 0.4)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("Mean ROC AUC across genus-blocked folds")
    ax.set_xticks(list(x))
    ax.set_xticklabels(rows["model"].astype(str), rotation=15, ha="right")
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    save(fig, "figure2_model_performance_auc.png")


def figure3_module_counts() -> None:
    modules = read_tsv("tier_module_counts_for_manuscript.tsv")
    pivot = modules.pivot_table(index="module_bin", columns="tier", values="n_families", aggfunc="sum").fillna(0)
    order = pivot.sum(axis=1).sort_values().index
    pivot = pivot.loc[order]

    fig, ax = plt.subplots(figsize=(6.5, 4.9))
    y = range(len(pivot))
    ms = pivot.get("Tier 1", pd.Series(0, index=pivot.index))
    hp = pivot.get("Tier 2", pd.Series(0, index=pivot.index))
    bars_ms = ax.barh([i - 0.2 for i in y], ms, height=0.38, color=TWO_GROUP[0], label="Model-Supported", edgecolor="white", linewidth=0.5)
    bars_hp = ax.barh([i + 0.2 for i in y], hp, height=0.38, color=TWO_GROUP[1], label="Highly Pure", edgecolor="white", linewidth=0.5)
    pad = pivot.values.max() * 0.01
    label_barh(ax, bars_ms, ms, "{:.0f}", pad)
    label_barh(ax, bars_hp, hp, "{:.0f}", pad)
    ax.set_yticks(list(y))
    ax.set_yticklabels([wrap_label(x, 28) for x in pivot.index], fontsize=9.5)
    ax.set_xlabel("Number of families")
    ax.set_xlim(0, pivot.values.max() * 1.12)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.legend(frameon=False, loc="lower right", fontsize=10)
    ax.set_axisbelow(True)
    save(fig, "figure3_tier_module_counts.png")


def figure4_top_candidates() -> None:
    top = read_tsv("tier1_top25_for_manuscript.tsv").head(15).copy()
    top = top.sort_values("consensus_rank_pct_mean")
    fig, ax = plt.subplots(figsize=(6.5, 5.6))
    colors = [PALETTE.get(m, OKABE["grey"]) for m in top["module_bin"]]
    labels = [f"{row.gene_family}  {short_product(row.product)}" for row in top.itertuples()]
    values = top["consensus_rank_pct_mean"].to_numpy()
    xmin = values.min() - 0.008
    xmax = min(1.0, values.max() + 0.012)
    for i, (value, color) in enumerate(zip(values, colors)):
        ax.plot([xmin, value], [i, i], color="#d3d9e0", linewidth=1.3, zorder=1)
    ax.scatter(values, range(len(top)), c=colors, s=95, edgecolor=INK, linewidth=0.6, zorder=2)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=9.5)
    ax.set_ylim(-0.6, len(top) - 0.4)
    ax.set_xlabel("Mean consensus rank percentile across the four classifiers")
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)

    present = [m for m in PALETTE if m in set(top["module_bin"])]
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=PALETTE[m], markeredgecolor=INK, markersize=8, label=m)
        for m in present
    ]
    ax.legend(
        handles=handles,
        title="Functional module",
        frameon=False,
        fontsize=8.8,
        title_fontsize=9.2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=2,
        handletextpad=0.4,
        columnspacing=1.2,
    )
    save(fig, "figure4_top_tier1_candidates.png")


def figure5_synteny_summary() -> None:
    cross = read_tsv("conserved_cluster_ncbi_product_crosswalk.tsv")
    green = OKABE["green"]
    light = "#cdd4dc"
    dpor_tag = {"GF_00326": "DPOR ATP-binding", "GF_01623": "DPOR subunit N"}

    def rows_for(module: str) -> list[dict]:
        out = []
        for _, r in cross[cross["module"] == module].iterrows():
            recon = [g for g in str(r["actual_pangenome_gene_families_for_product"]).split(";") if g]
            checked = int(r["representatives_checked"])
            exact = int(r["representatives_with_exact_product"])
            full = exact == checked and len(recon) == 1
            name = str(r["cluster_output_product_label"]).split()[-1]
            out.append({"recon": recon, "checked": checked, "exact": exact, "full": full, "name": name})
        return out

    def draw_box(x, y, hw, hh, fill, line1, line2):
        ax.add_patch(
            FancyBboxPatch(
                (x - hw, y - hh), 2 * hw, 2 * hh,
                boxstyle="round,pad=0.006,rounding_size=0.012",
                facecolor=fill, edgecolor=INK, linewidth=0.7,
                linestyle="solid" if fill == green else (0, (3, 2)),
            )
        )
        tcol = "white" if fill == green else "#1f2937"
        ax.text(x, y + 0.028, line1, ha="center", va="center", fontsize=9, fontweight="bold", color=tcol)
        ax.text(x, y - 0.032, line2, ha="center", va="center", fontsize=8.6, color=tcol)

    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # DPOR neighborhood on top, the prioritized candidate.
    ax.add_patch(plt.Rectangle((0.04, 0.64), 0.92, 0.24, facecolor="#eaf4ef", edgecolor="none", zorder=0))
    ax.text(0.06, 0.935, "DPOR conserved neighborhood", ha="left", va="center",
            fontsize=10.5, fontweight="bold", color="#0f172a")
    dpor = rows_for("cy71_tri28_ana67")
    dpor_x = [0.33, 0.67]
    for d, x in zip(dpor, dpor_x):
        gf = d["recon"][0]
        draw_box(x, 0.76, 0.15, 0.075, green, gf, dpor_tag.get(gf, "DPOR"))

    # Ribosomal and translation neighborhood below, the synteny control.
    ax.text(0.06, 0.55, "Ribosomal and translation neighborhood", ha="left", va="center",
            fontsize=10.5, fontweight="bold", color="#0f172a")
    ribo = rows_for("cy59_tri56_ana56")
    n = len(ribo)
    bw, gap = 0.165, 0.022
    start = (1.0 - (n * bw + (n - 1) * gap)) / 2
    for i, d in enumerate(ribo):
        cx = start + bw / 2 + i * (bw + gap)
        if d["full"]:
            draw_box(cx, 0.37, bw / 2, 0.075, green, d["recon"][0], d["name"])
        else:
            tag = f"{d['exact']} of {d['checked']}" if d["exact"] < d["checked"] else "split"
            draw_box(cx, 0.37, bw / 2, 0.075, light, d["name"], tag)

    # Legend.
    ax.add_patch(plt.Rectangle((0.06, 0.135), 0.028, 0.03, facecolor=green, edgecolor=INK, linewidth=0.6))
    ax.text(0.10, 0.15, "Model-Supported, shared across all three", ha="left", va="center", fontsize=8.8, color=INK)
    ax.add_patch(plt.Rectangle((0.55, 0.135), 0.028, 0.03, facecolor=light, edgecolor=INK, linewidth=0.6, linestyle=(0, (3, 2))))
    ax.text(0.59, 0.15, "split or representative-specific", ha="left", va="center", fontsize=8.8, color=INK)

    save(fig, "figure5_conserved_candidate_neighborhoods.png")


def figure6_functional_composition() -> None:
    comp = read_tsv("protein_vs_gene_top100_functional_composition.tsv")
    pivot = comp.pivot_table(index="functional_bin", columns="set_type", values="pct", aggfunc="sum").fillna(0)
    pivot = pivot.loc[pivot.max(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(9.4, 6.8))
    y = range(len(pivot))
    protein = pivot.get("Protein-driven", pd.Series(0, index=pivot.index))
    gene = pivot.get("Gene-driven", pd.Series(0, index=pivot.index))
    bars_p = ax.barh([i - 0.19 for i in y], protein, height=0.36, color=TWO_GROUP[0], label="Protein-family top 100", edgecolor="white", linewidth=0.5)
    bars_g = ax.barh([i + 0.19 for i in y], gene, height=0.36, color=TWO_GROUP[1], label="Gene-level top 100", edgecolor="white", linewidth=0.5)
    label_barh(ax, bars_p, protein, "{:.0f}", 0.4)
    label_barh(ax, bars_g, gene, "{:.0f}", 0.4)
    ax.set_yticks(list(y))
    ax.set_yticklabels([wrap_label(x, 34) for x in pivot.index], fontsize=9.2)
    ax.set_xlabel("Share of top-100 list (%)")
    ax.set_xlim(0, pivot.values.max() * 1.12)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.legend(frameon=False, loc="lower right")
    ax.set_axisbelow(True)
    save(fig, "figure6_protein_vs_gene_functional_composition.png")


def main() -> None:
    apply_style()
    figure1_pipeline()
    figure2_model_performance()
    figure3_module_counts()
    figure4_top_candidates()
    figure5_synteny_summary()
    figure6_functional_composition()


if __name__ == "__main__":
    main()
