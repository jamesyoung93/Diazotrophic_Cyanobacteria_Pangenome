# Reproducing the unified pipeline run

These instructions assume you are starting from a clean clone, with your working directory set to `unified_pipeline_clean`.

## 1. Prepare the environment

1. Create and activate the conda environment

```bash
mamba env create -f environment.yml
mamba activate pangenome_fox
```

2. Load external tool modules

```bash
module load hmmer/3.4
module load mmseqs2/15-6f452
```

3. Set your Entrez contact email

```bash
export ENTREZ_EMAIL="your.email@institution.edu"
```

## 2. Run the unified pipeline

```bash
chmod u+x run_unified_pipeline.sh run_postprocess_09_12.sh
./run_unified_pipeline.sh
```

This produces a run directory at `unified_pipeline_clean/unified_pipeline_run/`.

## 3. Run postprocessing

If you do not have outbound access to UniProt from compute nodes, skip UniProt GO mapping.

```bash
SKIP_UNIPROT_GO=1 ./run_postprocess_09_12.sh
```

If UniProt access is available, you can enable GO mapping by setting `SKIP_UNIPROT_GO=0`.
Downloaded protein payloads remain local to the run.

## 4. Capture provenance

```bash
bash scripts/capture_repro_metadata.sh
```

This writes a `reproducibility/` folder suitable for inclusion in an archival release.
