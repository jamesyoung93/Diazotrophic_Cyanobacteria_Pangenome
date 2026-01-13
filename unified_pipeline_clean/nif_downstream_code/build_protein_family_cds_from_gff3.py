#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple, Set, Optional
import zipfile

import pandas as pd


def run(cmd, cwd: Optional[Path] = None) -> None:
    print("[CMD]", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), cwd=str(cwd) if cwd else None, check=True)


def ensure_gff3_cache(
    run_dir: Path,
    datasets_cmd: str,
    include: str = "gff3",
    force: bool = False,
) -> Path:
    """
    Ensures a rehydrated NCBI datasets directory exists with GFF3 files.

    Returns the path to the NCBI datasets data directory:
      <run_dir>/downloads/genomes_pkg/ncbi_genomes/ncbi_dataset/data
    """
    genomes_file = run_dir / "genome_accessions.txt"
    if not genomes_file.exists():
        raise FileNotFoundError(f"Missing {genomes_file}. Run the upstream pipeline first.")

    pkg_dir = run_dir / "downloads" / "genomes_pkg"
    zip_path = pkg_dir / "ncbi_genomes.zip"
    extract_dir = pkg_dir / "ncbi_genomes"
    data_dir = extract_dir / "ncbi_dataset" / "data"

    # If already present and contains GFF3, reuse unless forced
    if data_dir.exists() and not force:
        gffs = list(data_dir.rglob("*.gff")) + list(data_dir.rglob("*.gff3"))
        if len(gffs) > 0:
            print(f"[OK] Using existing GFF3 cache: {data_dir} ({len(gffs)} gff files)")
            return data_dir

    pkg_dir.mkdir(parents=True, exist_ok=True)

    # Download dehydrated package with gff3 only
    run([
        datasets_cmd, "download", "genome", "accession",
        "--inputfile", str(genomes_file),
        "--dehydrated",
        "--include", include,
        "--filename", str(zip_path),
    ], cwd=run_dir)

    # Extract zip using Python to avoid external unzip dependency
    if extract_dir.exists():
        # Keep it simple: wipe old extract if forcing or stale
        # (You can remove this if you prefer incremental behavior.)
        import shutil
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Extracting {zip_path} -> {extract_dir}")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    # Rehydrate
    run([datasets_cmd, "rehydrate", "--directory", str(extract_dir)], cwd=run_dir)

    # Verify
    if not data_dir.exists():
        raise RuntimeError(f"Expected NCBI datasets data dir not found at {data_dir}")

    gffs = list(data_dir.rglob("*.gff")) + list(data_dir.rglob("*.gff3"))
    if len(gffs) == 0:
        raise RuntimeError(f"No GFF3 files found under {data_dir}. Download/include may have failed.")
    print(f"[OK] GFF3 cache ready: {data_dir} ({len(gffs)} gff files)")
    return data_dir


def parse_gff3_for_products(gff_root: Path, needed_proteins: Set[str]) -> Dict[str, Tuple[str, str, str]]:
    """
    Parse GFF3 to map protein_accession -> (product, gene, locus_tag).
    Only retains proteins in needed_proteins for performance.

    Returns dict: protein_accession -> (product, gene, locus_tag)
    """
    def parse_attrs(attr: str) -> Dict[str, str]:
        d = {}
        for part in attr.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                d[k] = v
        return d

    ann: Dict[str, Tuple[str, str, str]] = {}
    gffs = list(gff_root.rglob("*.gff")) + list(gff_root.rglob("*.gff3"))
    for gff in gffs:
        with gff.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue
                ftype = parts[2]
                if ftype not in ("CDS", "polypeptide"):
                    continue
                attrs = parse_attrs(parts[8])

                pid = attrs.get("protein_id") or attrs.get("proteinId") or attrs.get("Name")
                if not pid:
                    continue
                if pid not in needed_proteins:
                    continue
                if pid in ann:
                    continue

                product = attrs.get("product", "") or ""
                gene = attrs.get("gene", "") or ""
                locus = attrs.get("locus_tag", "") or attrs.get("locusTag", "") or ""
                ann[pid] = (product, gene, locus)

    return ann


def main() -> None:
    ap = argparse.ArgumentParser(description="Build protein_family_cds.tsv from NCBI GFF3 annotations (reproducible).")
    ap.add_argument("--run-dir", default="unified_pipeline_run", help="Path to unified_pipeline_run directory")
    ap.add_argument("--datasets-cmd", default="datasets", help="NCBI datasets CLI command (default: datasets)")
    ap.add_argument("--force-download", action="store_true", help="Redownload/rehydrate GFF3 even if cache exists")
    ap.add_argument("--no-download", action="store_true", help="Do not download GFF3; require existing cache")
    ap.add_argument("--gff-root", default=None, help="Optional override for GFF3 root directory")
    ap.add_argument("--map-file", default=None, help="Optional override for genome_protein_family_map.tsv path")
    ap.add_argument("--out", default=None, help="Output path for protein_family_cds.tsv")
    args = ap.parse_args()

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"ERROR: run_dir not found: {run_dir}")

    map_path = Path(args.map_file).resolve() if args.map_file else (run_dir / "genome_protein_family_map.tsv")
    if not map_path.exists():
        raise SystemExit(f"ERROR: missing {map_path}. Generate genome_protein_family_map.tsv first.")

    out_path = Path(args.out).resolve() if args.out else (run_dir / "protein_family_cds.tsv")

    # Locate or build GFF3 cache
    if args.gff_root:
        gff_root = Path(args.gff_root).resolve()
        if not gff_root.exists():
            raise SystemExit(f"ERROR: gff-root not found: {gff_root}")
    else:
        default_root = run_dir / "downloads" / "genomes_pkg" / "ncbi_genomes" / "ncbi_dataset" / "data"
        if args.no_download:
            gff_root = default_root
            if not gff_root.exists():
                raise SystemExit(f"ERROR: no-download set but cache not found at {gff_root}")
        else:
            gff_root = ensure_gff3_cache(
                run_dir=run_dir,
                datasets_cmd=args.datasets_cmd,
                include="gff3",
                force=args.force_download,
            )

    # Load mapping
    m = pd.read_csv(map_path, sep="\t", dtype=str)
    required_cols = {"gene_family", "protein_accession"}
    if not required_cols.issubset(set(m.columns)):
        raise SystemExit(f"ERROR: {map_path} missing required columns {required_cols}. Found: {list(m.columns)}")

    m["gene_family"] = m["gene_family"].fillna("").astype(str).str.strip()
    m["protein_accession"] = m["protein_accession"].fillna("").astype(str).str.strip()
    m = m[(m["gene_family"] != "") & (m["protein_accession"] != "")].copy()

    needed = set(m["protein_accession"].unique())
    print(f"[INFO] Need annotations for {len(needed)} unique proteins")

    ann = parse_gff3_for_products(gff_root, needed)
    print(f"[INFO] Found annotations for {len(ann)} proteins in GFF3")

    # Build annotation dataframe
    a = pd.DataFrame(
        [{"protein_accession": k, "product": v[0], "gene": v[1], "locus_tag": v[2]} for k, v in ann.items()],
        dtype=str,
    )

    out = m.merge(a, on="protein_accession", how="left")
    out["product"] = out["product"].fillna("")
    out["locus_tag"] = out["locus_tag"].fillna("")
    out["go_ids"] = ""  # UniProt step can populate later if available

    out2 = out[["gene_family", "protein_accession", "product", "go_ids", "locus_tag"]].copy()
    out2.to_csv(out_path, sep="\t", index=False)

    nonempty = (out2["product"].astype(str).str.strip() != "").sum()
    frac = nonempty / len(out2) if len(out2) else 0.0

    print(f"[OK] wrote {out_path}")
    print(f"[CHECK] Non-empty product rows: {nonempty} of {len(out2)} ({frac:.3f})")

    # Fail loudly if products are essentially missing
    if frac < 0.20:
        raise SystemExit(
            "ERROR: product annotation coverage is too low. "
            "Module binning will collapse to Other / unassigned.\n"
            "Verify GFF3 cache exists and that protein_accession IDs match GFF3 protein_id fields."
        )


if __name__ == "__main__":
    main()
