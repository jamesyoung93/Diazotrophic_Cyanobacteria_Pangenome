#!/usr/bin/env python3
import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, List

import pandas as pd


def find_scripts_dir(run_dir: Path) -> Path:
    """
    Locate nif_downstream_code containing 11_build_fox_gene_report.py and 12_make_results_tables.py.
    Works whether you run from unified_pipeline_clean or elsewhere.
    """
    candidates = [
        Path.cwd() / "nif_downstream_code",
        Path.cwd().parent / "nif_downstream_code",
        run_dir.parent / "nif_downstream_code",
        run_dir.parent.parent / "nif_downstream_code",
    ]
    for d in candidates:
        if (d / "11_build_fox_gene_report.py").is_file() and (d / "12_make_results_tables.py").is_file():
            return d.resolve()
    raise FileNotFoundError(
        "Could not locate nif_downstream_code with 11_build_fox_gene_report.py and 12_make_results_tables.py.\n"
        f"Tried: {[str(c) for c in candidates]}"
    )


def stage_modeling_artifacts(run_dir: Path, overwrite: bool = False) -> Path:
    """
    Ensure run_dir/results/modeling contains the artifacts Step 11 commonly scans.
    Copies from run_dir root if present.
    """
    src_root = run_dir
    dest = run_dir / "results" / "modeling"
    dest.mkdir(parents=True, exist_ok=True)

    # Copy these if they exist at run root
    explicit = [
        "gene_family_purity_stats.csv",
        "classification_summary.csv",
        "roc_curves.png",
        "feature_importance_plot.png",
    ]

    patterns = [
        "feature_importance_*.csv",
        "fold_metrics_*.csv",
        "feature_directionality_*.csv",
    ]

    def copy_one(src: Path, dst: Path):
        if not src.exists():
            return
        if dst.exists() and not overwrite:
            return
        shutil.copy2(src, dst)

    # Explicit files
    for name in explicit:
        copy_one(src_root / name, dest / name)

    # Patterned files
    for pat in patterns:
        for src in src_root.glob(pat):
            copy_one(src, dest / src.name)

    return dest


def run_step11_12(
    scripts_dir: Path,
    run_dir: Path,
    purity_threshold: float,
    max_members: int,
    top_n: int,
) -> None:
    step11 = scripts_dir / "11_build_fox_gene_report.py"
    step12 = scripts_dir / "12_make_results_tables.py"

    fox_report = run_dir / "fox_report"
    results_tables = run_dir / "results_tables"

    # Run Step 11
    cmd11 = [
        sys.executable, str(step11),
        "--root", str(run_dir),
        "--outdir", str(fox_report),
        "--purity_threshold", str(purity_threshold),
        "--max_members", str(max_members),
    ]
    print("\n[RUN] Step 11:", " ".join(cmd11))
    subprocess.run(cmd11, check=True)

    # Run Step 12
    cmd12 = [
        sys.executable, str(step12),
        "--report_dir", str(fox_report),
        "--outdir", str(results_tables),
        "--top_n", str(top_n),
    ]
    print("\n[RUN] Step 12:", " ".join(cmd12))
    subprocess.run(cmd12, check=True)


def file_nlines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return sum(1 for _ in f)


def write_tier2_counts(results_tables: Path) -> List[Path]:
    """
    Generate Tier2 counts from tier2_pure_positive_annotated.tsv if present.
    Writes:
      - tier2_module_bin_counts.tsv
      - tier2_fox_status_counts.tsv
    """
    tier2_annot = results_tables / "tier2_pure_positive_annotated.tsv"
    if not tier2_annot.exists():
        return []

    df = pd.read_csv(tier2_annot, sep="\t")

    out_paths: List[Path] = []

    # module_bin counts
    if "module_bin" in df.columns:
        counts = df["module_bin"].fillna("Other / unassigned").value_counts()
        out = counts.reset_index()
        out.columns = ["module_bin", "n_families"]
        out["percent"] = (out["n_families"] / out["n_families"].sum() * 100).round(1)
        p = results_tables / "tier2_module_bin_counts.tsv"
        out.to_csv(p, sep="\t", index=False)
        out_paths.append(p)
    else:
        print("[WARN] tier2_pure_positive_annotated.tsv has no 'module_bin' column; skipping module bin counts.")

    # fox_status counts
    if "fox_status" in df.columns:
        counts = df["fox_status"].fillna("Novel/uncertain candidate").value_counts()
        out = counts.reset_index()
        out.columns = ["fox_status", "n_families"]
        out["percent"] = (out["n_families"] / out["n_families"].sum() * 100).round(1)
        p = results_tables / "tier2_fox_status_counts.tsv"
        out.to_csv(p, sep="\t", index=False)
        out_paths.append(p)
    else:
        print("[WARN] tier2_pure_positive_annotated.tsv has no 'fox_status' column; skipping FOX status counts.")

    return out_paths


def clean_outputs(run_dir: Path) -> None:
    """
    Remove fox_report and results_tables to avoid stale files.
    """
    for d in [run_dir / "fox_report", run_dir / "results_tables"]:
        if d.exists():
            shutil.rmtree(d)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Rebuild Tier2 tables and counts by staging modeling artifacts then rerunning Step11/12."
    )
    ap.add_argument("--run-dir", default="unified_pipeline_run", help="Path to unified_pipeline_run directory")
    ap.add_argument("--purity-threshold", type=float, default=0.90, help="Purity threshold for Tier2 (default 0.90)")
    ap.add_argument("--max-members", type=int, default=10, help="Max example members per family (default 10)")
    ap.add_argument("--top-n", type=int, default=100, help="Top N for Tier1 table (default 100)")
    ap.add_argument("--clean", action="store_true", help="Delete fox_report/ and results_tables/ before regenerating")
    ap.add_argument("--overwrite-staged", action="store_true", help="Overwrite files in results/modeling if present")
    ap.add_argument("--only-counts", action="store_true", help="Do not rerun Step11/12; only compute Tier2 counts if Tier2 annotated exists")
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"ERROR: run dir not found: {run_dir}")

    results_tables = run_dir / "results_tables"
    fox_report = run_dir / "fox_report"

    if args.clean:
        print(f"[CLEAN] Removing {fox_report} and {results_tables}")
        clean_outputs(run_dir)

    if not args.only_counts:
        staged = stage_modeling_artifacts(run_dir, overwrite=args.overwrite_staged)
        print(f"[OK] Ensured staged modeling artifacts under: {staged}")

        scripts_dir = find_scripts_dir(run_dir)
        print(f"[OK] Using scripts dir: {scripts_dir}")

        run_step11_12(
            scripts_dir=scripts_dir,
            run_dir=run_dir,
            purity_threshold=args.purity_threshold,
            max_members=args.max_members,
            top_n=args.top_n,
        )

    # Diagnose Tier2 existence
    tier2_src = fox_report / "tier2_pure_positive_heldout.tsv"
    tier2_annot = results_tables / "tier2_pure_positive_annotated.tsv"

    src_lines = file_nlines(tier2_src)
    annot_lines = file_nlines(tier2_annot)

    print("\n[CHECK] Tier2 source:", tier2_src)
    print(f"[CHECK] Lines: {src_lines} (1 means header-only / empty)")

    print("[CHECK] Tier2 annotated:", tier2_annot)
    print(f"[CHECK] Lines: {annot_lines} (1 means header-only / empty)")

    # Write Tier2 counts if possible
    if tier2_annot.exists() and annot_lines > 1:
        out_paths = write_tier2_counts(results_tables)
        if out_paths:
            print("\n[OK] Wrote Tier2 count tables:")
            for p in out_paths:
                print(" -", p)
        else:
            print("\n[WARN] Tier2 annotated exists but expected columns missing; no count tables produced.")
    else:
        print(
            "\n[INFO] Tier2 annotated table is missing or empty.\n"
            "This can happen if no families meet the Tier2 criteria for this run, or if purity/model artifacts were not discoverable.\n"
            "You can inspect the Tier2 source file above and gene_family_purity_stats.csv under results/modeling."
        )


if __name__ == "__main__":
    main()
