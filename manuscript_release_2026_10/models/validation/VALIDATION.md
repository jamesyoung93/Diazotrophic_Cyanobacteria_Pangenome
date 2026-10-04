# Model validation

Fresh separate copies of the bundle regenerated all five analysis modes in the
existing Windows scientific runtime. All **38,706** OOF probabilities exactly
matched their saved references (maximum absolute difference **0**), with exact
feature masks and assignments and **2,784 independent metric checks**. The runs
performed **460 model fits**, each restricted to one CPU/BLAS thread.

The first fresh copy ran primary, repeated, sample-identity/group and current-status
modes. A second fresh copy tested the newly added official-taxonomy mode. Numeric
primary/repeated function ASTs and all sensitivity analysis source files in the
delivered bundle match the versions exercised in these copies; later dispatch
additions did not alter their numeric code. Small per-mode JSONs preserve the
evidence. `run_audit.py --mode verify --out reference` rechecks saved outputs without
fitting. Comparison tolerances remain rtol1e-10 and atol1e-12; no relaxation occurred.

The reference environment used Python3.10.16, numpy2.2.6, pandas2.3.3,
scikit-learn1.6.1, xgboost3.1.1 and threadpoolctl3.5.0. It reports SciPy module1.13.1
but distribution metadata1.13.0, which is explicitly disclosed. Same-runtime
fresh-copy verification is not a clean installation or cross-platform guarantee.

The installation specification now selects threadpoolctl3.6.0 to address a Windows
long-DLL-path failure encountered during a separately created isolated-environment
test. Its documented upstream fix is https://github.com/joblib/threadpoolctl/pull/189.
Clean-install validation, if completed, is reported separately and does not replace
the original reference-environment provenance.

## Separate isolated-environment result

All20 primary fits completed after upgrading threadpoolctl to3.6.0, which resolves
the observed native-crash mechanism consistent with the documented Windows
long-path bug. The **strict probability comparison remains FAIL**:140 of1704 LR
probabilities exceed the original rtol1e-10/atol1e-12, with maximum difference
4.66904e-10. No tolerance was changed. Exact masks, predicted classes, all reported
AUC/AP values and audit importance ranks agree. One nearly tied same-negative-class
LR pair reverses order; it does not change metrics. Reported numerical values agree
to four decimals. The reference SciPy uses MKL; the isolated wheel uses OpenBLAS.

See [clean-environment report](clean_environment/REPORT.txt), saved clean outputs,
runtime/native-library provenance and the independent
comparison tables. Only the primary run was fitted in the isolated environment;
the460-fit exact fresh-copy validation above used the original existing runtime.
The diagnostic compare_clean_backend.py can regenerate the comparison, and its
--assert-strict flag reproduces the preserved nonzero-exit assertion.

## Separate additional training-prevalence sensitivity

The optional run_training_prevalence.py command was run from a third fresh copy.
Its20 fits reproduce all1704 reference probabilities exactly (maximum difference0),
with144 independent metric checks, five nested masks and426 reconstructed LR
probabilities. This evidence is saved separately and is not added to the main
five-mode460-fit contract above. The optional arm applies absolute training
carrier prevalence>=40 and is not a replacement primary evaluation.
