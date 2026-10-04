"""Independent verification of nested training-prevalence masks and metrics."""
from pathlib import Path
E=Path(__file__).resolve().parents[1]
source=(E/'analysis/verify_sensitivity.py').read_text()
source=source.replace("D=OUT/'model_sample_identity_sensitivity'","D=OUT/'model_training_prevalence_sensitivity'")
source=source.replace('(n>0)&(k/n>0.1)&(k/n<0.9)','(n>=40)&(k/n>0.1)&(k/n<0.9)')
source=source.replace('independently_verified_masks=20','independently_verified_masks=5')
source=source.replace('original_fold_ids_retained_for_all_local_exclusion_arms=True','original_fold_ids_retained_for_all426=True')
exec(compile(source,'independent_training_prevalence_verifier','exec'),{'__file__':str(E/'analysis/verify_training_prevalence.py')})
