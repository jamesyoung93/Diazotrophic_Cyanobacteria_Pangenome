These are full-file replacements.

1) Replace:
   unified_pipeline_clean/run_unified_pipeline.sh
   unified_pipeline_clean/nif_downstream_code/pangenome_pipeline_consolidated2.py

2) Ensure executables:
   chmod u+x unified_pipeline_clean/run_unified_pipeline.sh

Key fixes:
- Downstream script path now points to unified_pipeline_clean/nif_downstream_code/... (not repo_root/nif_downstream_code)
- Robust script discovery that works after chdir and in extracted zip trees.
- Uses sys.executable for all internal step-script invocations.
