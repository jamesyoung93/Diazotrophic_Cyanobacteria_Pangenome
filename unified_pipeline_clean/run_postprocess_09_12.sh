#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ROOT="${ROOT}/unified_pipeline_run"
SCRIPTS="${ROOT}/nif_downstream_code"

# Set to 1 to skip UniProt/GO enrichment (useful if compute nodes have no internet)
SKIP_UNIPROT_GO="${SKIP_UNIPROT_GO:-0}"

# --- sanity checks ---
[[ -d "${RUN_ROOT}" ]] || { echo "ERROR: missing run dir: ${RUN_ROOT}"; exit 1; }
[[ -f "${SCRIPTS}/09_top_features_narratives.py" ]] || { echo "ERROR: missing ${SCRIPTS}/09_top_features_narratives.py"; exit 1; }
[[ -f "${SCRIPTS}/10_make_narrative_and_viz.py" ]] || { echo "ERROR: missing ${SCRIPTS}/10_make_narrative_and_viz.py"; exit 1; }
[[ -f "${SCRIPTS}/11_build_fox_gene_report.py" ]] || { echo "ERROR: missing ${SCRIPTS}/11_build_fox_gene_report.py"; exit 1; }
[[ -f "${SCRIPTS}/12_make_results_tables.py" ]] || { echo "ERROR: missing ${SCRIPTS}/12_make_results_tables.py"; exit 1; }

# --- find where the modeling outputs landed (feature_importance_*.csv) ---
IMP_FILE="$(find "${RUN_ROOT}" -type f -name 'feature_importance_*.csv' | head -n 1 || true)"
[[ -n "${IMP_FILE}" ]] || { echo "ERROR: no feature_importance_*.csv found under ${RUN_ROOT}"; exit 1; }
MODEL_DIR="$(cd "$(dirname "${IMP_FILE}")" && pwd)"
echo "[INFO] Using model directory: ${MODEL_DIR}"

# --- ensure gene_family_info.csv exists (for representative proteins) ---
GF_INFO="${RUN_ROOT}/gene_family_info.csv"
[[ -f "${GF_INFO}" ]] || { echo "ERROR: missing ${GF_INFO}"; exit 1; }

# --- Build protein->family maps (used by later reports as optional enrichment) ---
# Requires gene_families_clusters.tsv created by MMseqs
CLUSTERS="${RUN_ROOT}/gene_families_clusters.tsv"
if [[ -f "${CLUSTERS}" ]]; then
  echo "[STEP] Build protein_family_map.tsv + genome_protein_family_map.tsv"
  python "${SCRIPTS}/make_protein_family_map.py" \
    --clusters "${CLUSTERS}" \
    --family-info "${GF_INFO}" \
    --out "${RUN_ROOT}/protein_family_map.tsv" \
    --keep-genome-prefix

  # Expand into 3-column membership map for FOX report script (genome + protein + gene_family)
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{print "gene_family","genome_accession","protein_accession"; next}
    {
      split($1,a,"|");
      gf=$2; g=a[1]; p=a[2];
      if (g=="" || p=="") next;
      print gf,g,p
    }' "${RUN_ROOT}/protein_family_map.tsv" > "${RUN_ROOT}/genome_protein_family_map.tsv"
else
  echo "[WARN] ${CLUSTERS} not found; skipping protein_family_map generation."
fi

# --- Ensure protein_family_cds.tsv exists for Step 09 narrative table ---
# 09 requires protein_family_cds.tsv (it will error if missing) 
PF_CDS="$(find "${RUN_ROOT}" -maxdepth 3 -name 'protein_family_cds.tsv' | head -n 1 || true)"
if [[ -z "${PF_CDS}" ]]; then
  echo "[STEP] Create minimal protein_family_cds.tsv (from all_proteins.faa + protein_family_map.tsv)"
  [[ -f "${RUN_ROOT}/all_proteins.faa" ]] || { echo "ERROR: missing ${RUN_ROOT}/all_proteins.faa"; exit 1; }
  [[ -f "${RUN_ROOT}/protein_family_map.tsv" ]] || { echo "ERROR: missing ${RUN_ROOT}/protein_family_map.tsv (needed to build protein_family_cds.tsv)"; exit 1; }

  python - <<'PY'
from pathlib import Path
import pandas as pd

run = Path("unified_pipeline_run")
pf_map = run / "protein_family_map.tsv"
faa = run / "all_proteins.faa"
out = run / "protein_family_cds.tsv"

# Map protein_accession -> product (description) from FASTA headers
prod = {}
with faa.open() as f:
    for line in f:
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        toks = h.split()
        if not toks:
            continue
        pid = toks[0]                 # genome|protein or protein
        protein = pid.split("|",1)[1] if "|" in pid else pid
        desc = " ".join(toks[1:]).strip()
        prod[protein] = desc

m = pd.read_csv(pf_map, sep="\t")
# expected columns: protein_id, gene_family (from make_protein_family_map.py) 
m["protein_accession"] = m["protein_id"].astype(str).str.split("|", n=1, expand=True)[1].fillna(m["protein_id"].astype(str))
m["product"] = m["protein_accession"].map(prod).fillna("")
m["go_ids"] = ""
m["locus_tag"] = ""

out_df = m[["gene_family","protein_accession","product","go_ids","locus_tag"]].copy()
out_df.to_csv(out, sep="\t", index=False)
print(f"[OK] wrote {out}")
PY

  PF_CDS="${RUN_ROOT}/protein_family_cds.tsv"
fi
echo "[INFO] protein_family_cds: ${PF_CDS}"

# --- Optional: enrich with UniProt + GO via UniProt REST API ---
# add_uniprot_go.py uses UniProt ID mapping endpoints /idmapping/run|status|stream and UniProtKB stream. 
# It writes refseq_to_uniprot.tsv and uniprot_to_go.tsv and a merged TSV. 
PF_CDS_USED="${PF_CDS}"
if [[ "${SKIP_UNIPROT_GO}" -eq 0 && -f "${SCRIPTS}/add_uniprot_go.py" ]]; then
  echo "[STEP] UniProt/GO enrichment (internet required)"
  python "${SCRIPTS}/add_uniprot_go.py" \
    --in "${PF_CDS}" \
    --out "${RUN_ROOT}/protein_family_cds_uniprot_go.tsv" \
    --collapse-by-protein
  PF_CDS_USED="${RUN_ROOT}/protein_family_cds_uniprot_go.tsv"
else
  echo "[INFO] Skipping UniProt/GO enrichment (set SKIP_UNIPROT_GO=0 to enable)."
fi

# --- Step 09: top feature narrative table + module bins + outline ---
# Produces: top_features_narrative.tsv, module_summary.tsv, narrative_outline.md 
echo "[STEP] 09_top_features_narratives.py"
python "${SCRIPTS}/09_top_features_narratives.py" \
  --run-dir "${MODEL_DIR}" \
  --top-k 60 \
  --protein-family-cds "${PF_CDS_USED}" \
  --gene-family-info "${GF_INFO}"

# --- Step 10: figure pack + auto narrative markdown ---
# Uses top_features_narrative.tsv as input and writes narrative_viz outputs 
echo "[STEP] 10_make_narrative_and_viz.py"
python "${SCRIPTS}/10_make_narrative_and_viz.py" \
  --table "${MODEL_DIR}/top_features_narrative.tsv" \
  --top-k 60 \
  --outdir "${MODEL_DIR}/narrative_viz"

# --- Step 11: FOX report across results* directories under RUN_ROOT ---
# Writes tier tables + master table + plot 
echo "[STEP] 11_build_fox_gene_report.py"
python "${SCRIPTS}/11_build_fox_gene_report.py" \
  --root "${RUN_ROOT}" \
  --outdir "${RUN_ROOT}/fox_report" \
  --purity_threshold 0.90 \
  --max_members 10

# --- Step 12: manuscript-facing tables + bin counts ---
# Reads tier files from Step 11 and writes results_tables outputs 
echo "[STEP] 12_make_results_tables.py"
python "${SCRIPTS}/12_make_results_tables.py" \
  --report_dir "${RUN_ROOT}/fox_report" \
  --outdir "${RUN_ROOT}/results_tables" \
  --top_n 100

echo ""
echo "[DONE] Key outputs:"
echo " - ${MODEL_DIR}/top_features_narrative.tsv"
echo " - ${MODEL_DIR}/narrative_viz/ (panels + auto_narrative.md)"
echo " - ${RUN_ROOT}/fox_report/ (tier tables + master table)"
echo " - ${RUN_ROOT}/results_tables/ (manuscript tables)"
