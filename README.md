# Protein-family pangenome analysis of cyanobacterial diazotrophy

This repository contains the reproducible code and manuscript-facing result tables for a cyanobacterial pangenome analysis of diazotrophy beyond canonical `nif` genes.

The original unified workflow performs:

- upstream detection of `nifH`, `nifD`, and `nifK` marker hits using profile HMM search;
- assembly metadata enrichment and filtering to complete RefSeq (`GCF`) cyanobacterial genomes;
- MMseqs2 protein-family clustering and presence/absence matrix construction;
- supervised classification and feature ranking of protein families associated with diazotrophy;
- postprocessing to separate model-selected accessory families from near-diagnostic high-purity `nif`/context families.

## Current manuscript release

The current manuscript-facing June 2026 release is in:

`manuscript_release_2026_06/`

That folder contains the current tables and scripts used to align the core 426-genome protein-family model with biological interpretation layers:

- `476` Model-Supported accessory candidate families;
- `981` Highly Pure `nif`/near-core context families;
- related protein-family atlas product/function concordance;
- external literature proteomics evidence restricted to nitrogen-fixation-relevant active-phase-up responses;
- morphotype-breadth proxy and annotation review;
- HGT-proximity and alternative-nitrogenase audit summaries where available;
- condensate-driver ranking overlays used only as indirect prioritization evidence.

The local Cyanothece proteomics screen and FOX ensemble probability are retained only as historical/contextual analyses in older materials; they are not part of the current composite scoring used in the June 2026 manuscript release.

## Directory layout

- `manuscript_release_2026_06/`  
  Current manuscript-facing result tables, release notes, figure assets, and overlay scripts.
- `unified_pipeline_clean/`  
  Primary entrypoint for the core pangenome/modeling workflow.
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

These steps reflect the environment used during pipeline bring up on an Lmod based HPC cluster. If your site uses different module names, adapt accordingly.

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
grep '^GCF_' unified_pipeline_run/genome_accessions.txt > unified_pipeline_run/genome_accessions.gcf.txt
mv unified_pipeline_run/genome_accessions.gcf.txt unified_pipeline_run/genome_accessions.txt

./run_unified_pipeline.sh
```

### 4. Run postprocessing and manuscript tables

Some compute environments block outbound requests to UniProt. If UniProt is blocked, run with `SKIP_UNIPROT_GO=1` to bypass UniProt GO enrichment.

```bash
python nif_downstream_code/build_protein_family_cds_from_gff3.py --run-dir unified_pipeline_run

SKIP_UNIPROT_GO=1 ./run_postprocess_09_12.sh
python regen_tier2_tables.py --run-dir unified_pipeline_run --clean
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