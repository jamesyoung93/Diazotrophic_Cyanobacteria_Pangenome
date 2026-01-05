# Troubleshooting

## Permission denied when running a shell script

If you see `Permission denied`, make the script executable.

```bash
chmod u+x unified_pipeline_clean/run_unified_pipeline.sh
chmod u+x unified_pipeline_clean/run_postprocess_09_12.sh
```

## mmseqs not found

If you see `FileNotFoundError: ... 'mmseqs'`, load your MMseqs2 module or install MMseqs2 and ensure it is on PATH.

```bash
module avail mmseqs
module load mmseqs2/15-6f452
which mmseqs
```

## datasets not found

If you see errors about missing NCBI Datasets CLI, install it using conda.
NCBI install docs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ 

```bash
conda install -c conda-forge ncbi-datasets-cli
datasets --version
```

## UniProt GO mapping fails on HPC

Some institutional networks block outbound HTTPS from compute nodes or require a proxy.
If UniProt mapping fails, run postprocessing with `SKIP_UNIPROT_GO=1` and rely on the NCBI GFF annotations for product strings.
