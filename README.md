# FOX pangenome unified pipeline

This repository contains a reproducible, publication ready code package for the FOX pangenome analysis workflow.
The workflow integrates

- Upstream detection of nifH, nifD, nifK marker hits using profile HMM search
- Assembly quality metadata enrichment and filtering to RefSeq accessions and complete genomes
- Downstream pangenome construction using MMseqs2 gene family clustering
- Supervised classification and feature ranking to identify gene families associated with oxic nitrogen fixation
- Postprocessing scripts that generate tiered result tables, narrative summaries, and manuscript facing figures

The runtime workflow downloads public reference data at execution time and therefore does not include large data payloads in this code archive.

## Directory layout

- `unified_pipeline_clean/`  
  Primary entrypoint for end to end reproducible runs.
- `unified_pipeline_clean/nif_hdk_scan_release_clean/`  
  Upstream marker scan, hit summarization, and assembly metadata enrichment.
- `unified_pipeline_clean/nif_downstream_code/`  
  Downstream pangenome build, modeling, and postprocessing scripts.
- `tests/`  
  Lightweight unit tests for filtering and mode selection logic.
- `docs/`  
  Reproducibility and archival guidance.
- `scripts/`  
  Helper scripts for capturing environment metadata.

## Quick start for HPC environments

These steps reflect the environment used during pipeline bring up on an Lmod based HPC cluster.
If your site uses different module names, adapt accordingly.

### 1. Create and activate the conda environment

NCBI Datasets CLI is required to fetch genome packages and metadata. Installation options are documented by NCBI at https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/.


```bash
# load your site conda or mamba module first if needed
conda env create -f environment.yml
conda activate pangenome_fox
```

Alternatively, install `ncbi-datasets-cli` from conda-forge, or install the official binary from NCBI.


### 2. Load external tool modules

The pipeline expects the following command line tools to be available on PATH.

```bash
module load hmmer/3.4
module load mmseqs2/15-6f452
```

Confirm availability.

```bash
which hmmsearch
which mmseqs
datasets --version
python -V
```

### 3. Run the unified pipeline

Set the NCBI Entrez email required by the download utilities.

```bash
export ENTREZ_EMAIL="your.email@institution.edu"
```

Run from the unified pipeline directory.

```bash
cd unified_pipeline_clean
chmod u+x run_unified_pipeline.sh run_postprocess_09_12.sh
./run_unified_pipeline.sh
```

### 4. Run postprocessing and manuscript tables

Some compute environments block outbound requests to UniProt.
If UniProt is blocked, run with `SKIP_UNIPROT_GO=1` to bypass UniProt GO enrichment.

```bash
python nif_downstream_code/build_protein_family_cds_from_gff3.py --run-dir unified_pipeline_run

SKIP_UNIPROT_GO=1 ./run_postprocess_09_12.sh
```

Outputs are written under `unified_pipeline_clean/unified_pipeline_run/`.

## Citation and archival

- GitHub citation files: https://docs.github.com/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files

- Zenodo GitHub release archiving: https://help.zenodo.org/docs/github/archive-software/github-upload/
- Zenodo software metadata guidance: https://help.zenodo.org/docs/github/describe-software/


See `docs/ZENODO_GITHUB.md` for a step by step release workflow.

## Reproducibility checklist

Before creating an archival release, capture the environment used for the run.

```bash
bash scripts/capture_repro_metadata.sh
```

This writes `reproducibility/` artifacts such as module lists and package inventories.

