#!/usr/bin/env python
"""Portable sample/group/current-status sensitivities; no canonical mutations."""
from pathlib import Path
import argparse,os,sys,subprocess,json,hashlib,importlib.metadata
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parent
def verify_reference(out,scope):
    actual=pd.read_csv(out/scope/'oof_predictions.csv');expected=pd.read_csv(BASE/'reference'/scope/'oof_predictions.csv')
    keys=['scenario','model','fold','assembly_accession']
    merged=actual.merge(expected,on=keys,suffixes=('_new','_ref'),validate='one_to_one')
    assert len(merged)==len(actual)==len(expected)
    np.testing.assert_allclose(merged.y_prob_new,merged.y_prob_ref,rtol=1e-10,atol=1e-12)
    np.testing.assert_array_equal(merged.y_true_new,merged.y_true_ref)
    np.testing.assert_array_equal(merged.y_pred_new,merged.y_pred_ref)
    for name in ['fold_assignments.csv','fold_family_features.csv']:
        pd.testing.assert_frame_equal(pd.read_csv(out/scope/name),pd.read_csv(BASE/'reference'/scope/name))
    if scope=='current_status':
        status=pd.read_csv(out/scope/'frozen_current_status_with_labels.csv').set_index('assembly_accession')
        independent=pd.read_csv(BASE/'inputs/public_status/cohort_status_reconciled.csv').set_index('assembly_accession')
        assert status.status.sort_index().equals(independent.public_status.sort_index())
        assert status.status.value_counts().to_dict()=={'current':401,'suppressed':22,'previous':3}
        assign=pd.read_csv(out/scope/'fold_assignments.csv')
        for scenario,keep in [('exclude_currently_suppressed22',status.status!='suppressed'),('current_assemblies_only401',status.status=='current')]:
            rows=assign[assign.scenario==scenario]
            assert set(rows.assembly_accession)==set(status.index[keep]) and rows.y_true.sum()==106
    result=dict(status='PASS',reference_prediction_rows=len(actual),maximum_reference_probability_difference=float((merged.y_prob_new-merged.y_prob_ref).abs().max()),exact_fold_assignments_and_feature_masks_match=True)
    (out/scope/'portable_reference_verification.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(scope,result,flush=True)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['sample_identity','current_status','official_taxonomy','all','verify'],default='all')
    parser.add_argument('--out',type=Path,default=BASE/'rerun_outputs')
    args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((BASE/'input_manifest.json').read_text())
    for name,digest in manifest['sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest,name
    environment=os.environ.copy();environment['MODEL_AUDIT_OUTPUT']=str(out)
    for scope in ['sample_identity','current_status','official_taxonomy']:
        if args.mode not in [scope,'all','verify']:continue
        if args.mode!='verify':
            with (out/(scope+'_training.log')).open('w',encoding='utf8') as log:
                subprocess.run([sys.executable,str(BASE/'analysis'/f'{scope}.py')],env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
        elif not (out/scope/'oof_predictions.csv').exists():continue
        environment['MODEL_AUDIT_SCOPE']=scope
        with (out/(scope+'_verification.log')).open('w',encoding='utf8') as log:
            subprocess.run([sys.executable,str(BASE/'analysis'/('verify_official_taxonomy.py' if scope=='official_taxonomy' else 'verify_sensitivity.py'))],env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
        verify_reference(out,scope)
    runtime=dict(python=sys.version,distribution_metadata={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','scikit-learn','xgboost','threadpoolctl']})
    (out/'sensitivity_runtime.json').write_text(json.dumps(runtime,indent=2),encoding='utf8')
if __name__=='__main__':main()
