#!/usr/bin/env python3
"""
editorial_diagnostics.py
------------------------
Reproduces every quantitative test behind 07_EDITORIAL_ARCHITECTURE.md.
Each test targets a specific claim in the "Beyond nif" manuscript.

Usage:
    python scripts/historical_diagnostics.py [--workbook PATH --matrix PATH --labels PATH]

Requires: pandas, numpy, scipy, statsmodels, openpyxl (matplotlib for --figure)
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, spearmanr
import statsmodels.formula.api as smf
import warnings; warnings.filterwarnings("ignore")

T = lambda s: s.astype(str).str.lower().isin(["true", "1", "1.0"])
ABUNDANT = ("clp|chaperon|dnak|groe|grol|dnaj|grpe|release factor|translation|elongation factor|"
            "typa|trna|nyn|ribosom|isopropylmalate|threonine synthase|ketol-acid|dihydroxy-acid|"
            "aminomutase|amino acid")
CARBON = (r"glucose-6-phosphate|phosphogluconate|transketolase|transaldolase|isocitrate|"
          r"glyceraldehyde-3-phosphate|fructose-bisphosphate|fructose-1,6-bisphosphat|phosphoglycerate|enolase|"
          r"pyruvate|citrate synthase|aconit|malate|fumarate|succinyl|2-oxoglutarate|phosphoenolpyruvate|"
          r"glycogen|ribulose|sedoheptulose|ribose-5-phosphate|triosephosphate|glucokinase|hexokinase|"
          r"phosphofructo|acetyl-coa")
GENERIC = r"hypothetical|DUF\d|uncharacteri|domain-containing protein$|family protein$"

def hdr(s): print("\n" + "=" * 92 + "\n" + s + "\n" + "=" * 92)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    package=Path(__file__).resolve().parents[1]
    ap.add_argument('--workbook',type=Path,default=package/'inputs/released/proteomics_composite_family_evidence.xlsx')
    ap.add_argument('--matrix',type=Path,default=package/'inputs/raw/gene_family_matrix.csv')
    ap.add_argument('--labels',type=Path,default=package/'inputs/raw/complete_genomes_labeled.csv')
    ap.add_argument('--figure')
    a=ap.parse_args()
    ac=pd.read_excel(a.workbook,sheet_name='AllComposite')
    tf=pd.read_excel(a.workbook,sheet_name='TopFamilies')
    print('HISTORICAL RELEASE DIAGNOSTICS: frozen1457 inventory, original31 TierA and174 shortlist. Not current167 results.')
    ms = ac[ac.candidate_set == "Model-Supported"].copy()
    tier = ms.priority_tier.fillna("")

    # ------------------------------------------------------------------ D1
    hdr("D1  Table 2 'Other / unassigned' is not an uncharacterized bin")
    oth = ms[ms.module_bin == "Other / unassigned"]
    gen = oth["product"].str.contains(GENERIC, case=False, na=False)
    print(f"  Model-Supported 'Other / unassigned': {len(oth)}")
    print(f"    specific named function : {(~gen).sum()} ({(~gen).mean():.0%})")
    print(f"    generic / hypothetical  : {gen.sum()} ({gen.mean():.0%})")
    allgen = ms["product"].str.contains(GENERIC, case=False, na=False)
    print(f"  truly generic/hypothetical names across ALL Model-Supported bins: {allgen.sum()} of {len(ms)} ({allgen.mean():.0%})")

    hdr("D2  Carbon metabolism is under-counted ~8x by the keyword binner")
    c = ms["product"].str.contains(CARBON, case=False, na=False)
    print(f"  Model-Supported central carbon-metabolism enzymes by product name: {c.sum()}")
    print(f"  of which binned 'Carbon metabolism & NAD(P)H supply': {(c & (ms.module_bin=='Carbon metabolism & NAD(P)H supply')).sum()}")
    print(f"  landed in: {ms[c].module_bin.value_counts().to_dict()}")

    hdr("D3  Substring bug: 'fur' (ferric uptake regulator) matches inside 'sulfur'")
    reg = ac[ac.module_bin == "Regulation"]
    bug = reg["product"].str.contains("sulfur|desulfur", case=False, na=False) & \
          ~reg["product"].str.contains("regulat|ferric uptake", case=False, na=False)
    print(f"  sulfur-metabolism enzymes filed as 'Regulation': {bug.sum()}")
    for g, p in reg[bug][["gene_family", "product"]].values: print(f"    {g}  {p}")

    # ------------------------------------------------------------------ D4
    hdr("D4  What the composite score is made of (Tier A + B, Model-Supported)")
    top = tf[tf.priority_tier.str.startswith(("Tier A", "Tier B"))]
    tot = top.scope_adjusted_story_score.mean()
    for k in ["model_score", "morphotype_breadth_score", "literature_proteomics_score",
              "related_atlas_score", "condensate_score"]:
        print(f"  {k:<30} mean {top[k].mean():.2f}  ({top[k].mean()/tot:.0%} of score)")
    bare = (top.literature_proteomics_score == 0) & (top.related_atlas_score == 0) & (top.condensate_score == 0)
    tb = top.priority_tier.str.startswith("Tier B")
    print(f"  Tier B 'strong cross-evidence' families with ZERO external evidence: {(bare & tb).sum()} of {tb.sum()}")
    ta = tf[tf.priority_tier.str.startswith("Tier A")]
    print(f"  Tier A families passing both morphotype flags: {(ta.morphotype_breadth_score==2.75).sum()} of {len(ta)} (saturated gate)")
    r, p = spearmanr(ms.consensus_rank_pct_mean, ms.scope_adjusted_story_score, nan_policy="omit")
    print(f"  Spearman(model consensus rank, composite score): rho={r:+.2f}  -> score is nearly blind to model rank")
    print(f"  median model rank pct: Tier A {ms.loc[tier.str.startswith('Tier A'),'consensus_rank_pct_mean'].median():.2f}"
          f" vs other MS {ms.loc[~tier.str.startswith(('Tier A','Tier B')),'consensus_rank_pct_mean'].median():.2f}")

    # ------------------------------------------------------------------ D5
    hdr("D5  Proteomics annotation association conditional on detection opportunities")
    mm = ms[ms.n_literature_mapped_studies >= 1].copy()
    mm["up"] = (mm.n_literature_active_phase_up_studies >= 1).astype(int)
    mm["strict"] = T(mm.strict_unicellular_breadth_ge3u_ge1f).astype(int)
    mm["nmap"] = mm.n_literature_mapped_studies; mm["rank"] = mm.consensus_rank_pct_mean
    f = smf.logit("up ~ strict + nmap + rank", data=mm.dropna(subset=["rank"])).fit(disp=0)
    for k in ["strict", "nmap", "rank"]:
        print(f"  logit up ~ ... : {k:<7} OR={np.exp(f.params[k]):.2f}  p={f.pvalues[k]:.3g}")
    mm["per_opp"] = mm.n_literature_active_phase_up_studies / mm.nmap
    print(f"  mean family-level active-up/mapped-study ratio: strict {mm[mm.strict==1].per_opp.mean():.3f}  vs  not strict {mm[mm.strict==0].per_opp.mean():.3f}")
    print(f"  overall mean family-level ratio: {mm.per_opp.mean():.3f}  (descriptive family-weighted ratio; not a per-study probability)")
    base = (ac.n_literature_active_phase_up_studies >= 1).mean()
    print(f"  => '71 of 174' ({(tf.n_literature_active_phase_up_studies>=1).mean():.0%}) vs {base:.1%} baseline is not evidence of convergence")

    # ------------------------------------------------------------------ D6
    hdr("D6  Condensate 'selectivity' for Model-Supported is confounded with prevalence")
    cc = ac[ac.condensate_mapped_rows.fillna(0) >= 1].copy()
    cc["ms"] = (cc.candidate_set == "Model-Supported").astype(int); cc["top10"] = T(cc.condensate_top10).astype(int)
    cc["nrows"] = cc.condensate_mapped_rows
    f1 = smf.logit("top10 ~ ms + np.log(nrows)", data=cc).fit(disp=0)
    f2 = smf.logit("top10 ~ ms + np.log(nrows) + np.log(n_member_genomes)", data=cc).fit(disp=0)
    print(f"  Model-Supported OR, controlling #mapped proteins only : {np.exp(f1.params['ms']):.2f} (p={f1.pvalues['ms']:.2g})")
    print(f"  Model-Supported OR, also controlling carrier genomes   : {np.exp(f2.params['ms']):.2f} (p={f2.pvalues['ms']:.2g})")
    print(f"  carrier-genome prevalence OR                           : {np.exp(f2.params['np.log(n_member_genomes)']):.2f} (p={f2.pvalues['np.log(n_member_genomes)']:.2g})")
    one = cc[cc.nrows == 1].groupby("ms").top10.mean()
    print(f"  (descriptive single-mapped-row stratum:, MS {one[1]:.1%} vs HP {one[0]:.1%})")
    med = cc.groupby("candidate_set").n_member_genomes.median()
    print(f"  median carrier genomes: Model-Supported {med['Model-Supported']:.0f}, Highly Pure {med['Highly Pure']:.0f}"
          f" (HP max {cc[cc.ms==0].n_member_genomes.max():.0f}) -> inventories barely overlap in prevalence")
    print("  within matched prevalence strata (one-sided Fisher, uncorrected):")
    for lo, hi in [(40, 60), (60, 80), (80, 130)]:
        s_ = cc[(cc.n_member_genomes >= lo) & (cc.n_member_genomes < hi)]
        x, y = s_[s_.ms == 1].top10.astype(bool), s_[s_.ms == 0].top10.astype(bool)
        OR, pv = fisher_exact([[x.sum(), (~x).sum()], [y.sum(), (~y).sum()]], alternative="greater")
        print(f"    {lo:>3}-{hi:<3} carriers: MS {x.mean():5.1%} ({x.sum()}/{len(x)})  HP {y.mean():5.1%} ({y.sum()}/{len(y)})  OR={OR:.2f} p={pv:.3g}")
    print("  => contrast is confounded with prevalence and inconsistent within matched strata;")
    print("     the manuscript's 'selectivity is biologically informative' claim needs a prevalence-matched test")

    # ------------------------------------------------------------------ D7
    hdr("D7  Functional composition shifts with proteomic detectability")
    ms["abund"] = ms["product"].str.lower().fillna("").str.contains(ABUNDANT)
    ms["nmap"] = ms.n_literature_mapped_studies.fillna(0).astype(int)
    for n, g in ms.groupby("nmap"):
        print(f"  detected in {n} proteomes: {g.abund.mean():4.0%} proteostasis/translation/amino-acid  (n={len(g)})")
    print(f"  Tier A: {ms.loc[tier.str.startswith('Tier A'),'abund'].mean():.0%}")

    # ------------------------------------------------------------------ D8
    hdr("D8  Within- versus between-role co-occurrence correlation (exploratory permutation)")
    M = pd.read_csv(a.matrix, index_col=0)
    lab = pd.read_csv(a.labels)
    mz = ms[(ms.story_role != "Other / unassigned") & ms.gene_family.isin(M.columns)].drop_duplicates("gene_family")
    k = mz.story_role.value_counts(); mz = mz[mz.story_role.isin(k[k >= 5].index)]
    role = dict(zip(mz.gene_family, mz.story_role.astype(str))); rng = np.random.default_rng(0)
    dz = set(lab.loc[lab.is_diazotroph == True, "assembly_accession"])
    for label, X in [("all 426 genomes", M), ("112 diazotrophs", M.loc[[i for i in M.index if i in dz]])]:
        fam = [g for g in mz.gene_family if X[g].std() > 0]
        C = np.corrcoef(X[fam].values.astype(float).T); L = np.array([role[g] for g in fam], dtype=object)
        iu = np.triu_indices(len(fam), 1); v = C[iu]
        d = lambda l: (lambda s: v[s].mean() - v[~s].mean())((l[:, None] == l[None, :])[iu])
        obs = d(L); null = np.array([d(rng.permutation(L)) for _ in range(2000)])
        print(f"  {label:<16} within-minus-between r = {obs:+.3f}   permutation p = {(np.sum(null>=obs)+1)/2001:.3f}")

    if a.figure: make_figure(ac, ms, a.figure)

def make_figure(ac, ms, out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
    S1, S2 = "#2a78d6", "#eb6834"   # validated categorical slots 1-2
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.5), facecolor=SURF)
    for x in ax: x.set_facecolor(SURF); x.yaxis.grid(True, color=GRID, lw=0.8); x.set_axisbelow(True)

    # A: per-stratum active-phase-up rate by breadth
    mm = ms[ms.n_literature_mapped_studies >= 1].copy()
    mm["up"] = mm.n_literature_active_phase_up_studies >= 1
    mm["strict"] = T(mm.strict_unicellular_breadth_ge3u_ge1f)
    for flag, col, name in [(True, S1, "Broad (strict cross-morphotype)"), (False, S2, "Not broad")]:
        g = mm[mm.strict == flag].groupby("n_literature_mapped_studies").up.agg(["mean", "size"])
        g = g[g["size"] >= 10]
        ax[0].plot(g.index, 100 * g["mean"], "-o", color=col, lw=2, ms=6, mec=SURF, mew=1.5, label=name)
    ax[0].set_xlabel("Proteomes the family is detected in"); ax[0].set_ylabel("Families active-phase-up (%)")
    ax[0].set_title("A  Proteomics support tracks detection,\n    not breadth", loc="left", fontsize=10, color=INK)
    ax[0].set_xticks([1, 2, 3, 4]); ax[0].set_ylim(0, 80); ax[0].legend(frameon=False, fontsize=8, loc="upper left")

    # B: condensate top-10% by prevalence, per inventory
    cc = ac[ac.condensate_mapped_rows.fillna(0) >= 1].copy(); cc["top10"] = T(cc.condensate_top10)
    edges = [40, 60, 80, 110, 426]; cc["pbin"] = pd.cut(cc.n_member_genomes, edges, include_lowest=True)
    for cs, col, name in [("Model-Supported", S1, "Model-Supported"), ("Highly Pure", S2, "Highly Pure")]:
        g = cc[cc.candidate_set == cs].groupby("pbin", observed=True).top10.agg(["mean", "size"])
        g = g[g["size"] >= 15]; xs = [edges.index(int(iv.left) if iv.left >= 40 else 40) for iv in g.index]
        ax[1].plot(xs, 100 * g["mean"], "-o", color=col, lw=2, ms=6, mec=SURF, mew=1.5, label=name)
    ax[1].set_xticks(range(4)); ax[1].set_xticklabels(["40–60", "60–80", "80–110", "110+"])
    ax[1].set_xlabel("Carrier genomes per family"); ax[1].set_ylabel("In top 10% condensate rank (%)")
    ax[1].set_title("B  The inventory contrast is\n    a prevalence contrast", loc="left", fontsize=10, color=INK)
    ax[1].text(3, 37, "median carriers:\nHighly Pure 55\nModel-Supported 117", ha="right", va="top", fontsize=7.5, color=INK2)
    ax[1].legend(frameon=False, fontsize=8, loc="upper left"); ax[1].set_ylim(0, 40)

    # C: detectability gradient in functional composition
    ms2 = ms.copy(); ms2["abund"] = ms2["product"].str.lower().fillna("").str.contains(ABUNDANT)
    g = ms2.groupby(ms2.n_literature_mapped_studies.fillna(0).astype(int)).abund.mean()
    ax[2].plot(g.index, 100 * g.values, "-o", color=S1, lw=2, ms=6, mec=SURF, mew=1.5)
    ta = 100 * ms2.loc[ms2.priority_tier.fillna("").str.startswith("Tier A"), "abund"].mean()
    ax[2].axhline(ta, color=INK2, lw=1, ls=(0, (4, 3)))
    ax[2].text(0, ta + 2, f"Tier A: {ta:.0f}%", color=INK2, fontsize=8)
    ax[2].set_xticks(range(5)); ax[2].set_ylim(0, 70)
    ax[2].set_xlabel("Proteomes the family is detected in"); ax[2].set_ylabel("Proteostasis / translation /\namino-acid families (%)")
    ax[2].set_title("C  Functional mix shifts\n    with detectability", loc="left", fontsize=10, color=INK)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(out, dpi=300, facecolor=SURF); print(f"\nfigure written: {out}")

if __name__ == "__main__":
    sys.exit(main())
