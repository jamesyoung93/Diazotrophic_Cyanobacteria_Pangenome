#!/usr/bin/env python
"""Separate additional20-fit training-carrier-prevalence>=40 sensitivity."""
from pathlib import Path
import argparse,os,sys,json,hashlib,subprocess
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out',type=Path,default=BASE/'extra_prevalence_outputs')
parser.add_argument('--verify-only',action='store_true')
args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
for name,digest in json.loads((BASE/'input_manifest.json').read_text())['sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
environment=os.environ.copy();environment['MODEL_AUDIT_OUTPUT']=str(out);environment['MODEL_AUDIT_SCOPE']='training_prevalence'
if not args.verify_only:
    with (out/'training_prevalence_training.log').open('w',encoding='utf8') as log:
        subprocess.run([sys.executable,str(BASE/'analysis/training_prevalence.py')],env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
with (out/'training_prevalence_verification.log').open('w',encoding='utf8') as log:
    subprocess.run([sys.executable,str(BASE/'analysis/verify_training_prevalence.py')],env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
D=out/'training_prevalence';reference=BASE/'reference/training_prevalence'
p=pd.read_csv(D/'oof_predictions.csv');r=pd.read_csv(reference/'oof_predictions.csv')
keys=['scenario','model','fold','assembly_accession'];j=p.merge(r,on=keys,suffixes=('_new','_ref'),validate='one_to_one');assert len(j)==len(p)==len(r)==1704
np.testing.assert_allclose(j.y_prob_new,j.y_prob_ref,rtol=1e-10,atol=1e-12)
np.testing.assert_array_equal(j.y_pred_new,j.y_pred_ref);np.testing.assert_array_equal(j.y_true_new,j.y_true_ref)
for name in ['fold_assignments.csv','fold_family_features.csv']:pd.testing.assert_frame_equal(pd.read_csv(D/name),pd.read_csv(reference/name),check_exact=True)
result=dict(status='PASS',scope='Separate additional20-fit sensitivity; main all-mode460-fit contract unchanged.',reference_prediction_rows=len(j),maximum_reference_probability_difference=float((j.y_prob_new-j.y_prob_ref).abs().max()),exact_assignments_and_masks_match=True)
(D/'portable_reference_verification.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result,indent=2))
