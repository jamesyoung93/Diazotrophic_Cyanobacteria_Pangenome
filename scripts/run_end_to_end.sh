#!/usr/bin/env bash
set -euo pipefail

# End to end helper.
# Assumes you already created and activated the conda environment and loaded required modules.

export ENTREZ_EMAIL="${ENTREZ_EMAIL:-}"
if [[ -z "${ENTREZ_EMAIL}" ]]; then
  echo "ERROR: ENTREZ_EMAIL is not set. Example:"
  echo "  export ENTREZ_EMAIL="your.email@institution.edu""
  exit 1
fi

cd unified_pipeline_clean

chmod u+x run_unified_pipeline.sh run_postprocess_09_12.sh

./run_unified_pipeline.sh

# Default to skipping UniProt on HPC unless explicitly enabled
SKIP_UNIPROT_GO="${SKIP_UNIPROT_GO:-1}" ./run_postprocess_09_12.sh
