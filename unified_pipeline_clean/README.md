# Sequence acquisition and family construction

This folder supplies acquisition, HMM scanning, assembly filtering and MMseqs2 clustering code underlying the fixed manuscript matrix. The manuscript's prediction and candidate-score results are reproduced by the [release workflow](../manuscript_release_2026_10/README.md).

## Run on an HPC or Linux environment

From this directory:

```sh
conda env create -f environment.yml
conda activate pangenome_fox
export ENTREZ_EMAIL="your.email@institution.edu"
./run_unified_pipeline.sh
```

HMMER, MMseqs2 and NCBI Datasets must be available on PATH. See [installation](docs/INSTALL_HPC.md), [workflow details](docs/REPRODUCE_PIPELINE.md) and [troubleshooting](docs/TROUBLESHOOTING.md).

The driver can download current assemblies and write a new working directory. Such a run changes the data collection and does not replace the manuscript's frozen 426-assembly panel. Its downstream discovery exports supply provenance; use the release's model commands for training-fold feature filtering and reported prediction sensitivities.

Generated downloads, scan logs and earlier result directories are excluded from the current tree. Exact source revisions for the fixed analysis are recorded in [the model input manifest](../manuscript_release_2026_10/models/input_manifest.json).
