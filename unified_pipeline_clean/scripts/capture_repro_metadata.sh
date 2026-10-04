#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="reproducibility"
mkdir -p "${OUT_DIR}"

stamp="$(date +%Y%m%d_%H%M%S)"
run_dir="${OUT_DIR}/${stamp}"
mkdir -p "${run_dir}"

echo "[INFO] Writing reproducibility metadata to: ${run_dir}"

# Module list (Lmod)
if command -v module >/dev/null 2>&1; then
  # module list writes to stderr on many systems
  module list 2> "${run_dir}/module_list.txt" || true
fi

# Tool versions
( python -V ) > "${run_dir}/python_version.txt" 2>&1 || true
( datasets --version ) > "${run_dir}/datasets_version.txt" 2>&1 || true
( mmseqs version ) > "${run_dir}/mmseqs_version.txt" 2>&1 || true
( hmmsearch -h | head -n 5 ) > "${run_dir}/hmmer_hmmsearch_version.txt" 2>&1 || true

# Python packages
( python -m pip freeze --all ) > "${run_dir}/pip_freeze_all.txt" 2>&1 || true

# Conda environment export if available
if command -v conda >/dev/null 2>&1; then
  conda env export > "${run_dir}/conda_env_export.yml" 2>/dev/null || true
fi

echo "[INFO] Done."
