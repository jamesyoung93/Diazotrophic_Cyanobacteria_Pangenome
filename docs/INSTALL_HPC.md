# HPC installation notes

This workflow was executed on an Lmod based HPC environment. The exact module names vary by site.

## External tool dependencies

The downstream steps require these executables on PATH.

- hmmsearch from HMMER
- mmseqs from MMseqs2
- datasets from NCBI Datasets CLI

If your cluster provides these as modules, load them before running.

```bash
module load hmmer/3.4
module load mmseqs2/15-6f452
```

If `datasets` is not available as a module, install it in your conda environment.
NCBI install docs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ 

## Conda environment

The provided `environment.yml` creates a single environment that includes both the Python dependencies and the NCBI Datasets CLI.

```bash
mamba env create -f environment.yml
mamba activate pangenome_fox
```

If your site does not allow conda, see the NCBI documentation for downloading prebuilt binaries and adding them to PATH. 
