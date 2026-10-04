"""Training-fold carrier prevalence>=40 sensitivity, original folds and defaults.

Read-only scientific sensitivity. Original matrix, labels and inventory stay frozen.
Raw public API snapshots are inputs; this script makes no network requests.
"""
from pathlib import Path
import os, ast, json, hashlib, sys, time
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name]='1'
import numpy as np
import pandas as pd
import scipy, sklearn, xgboost
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score
from threadpoolctl import threadpool_limits

E=Path(__file__).resolve().parents[1]
D=Path(os.environ['MODEL_AUDIT_OUTPUT'])/'training_prevalence'; D.mkdir(parents=True,exist_ok=True)
R=E
P=E/'inputs'
PUBLIC=E/'inputs/public_status'
X=pd.read_csv(P/'gene_family_matrix.csv',index_col=0)
M=pd.read_csv(P/'complete_genomes_with_proteins.csv').drop_duplicates('assembly_accession').set_index('assembly_accession')
ids=sorted(set(X.index)&set(M.index)); X=X.loc[ids]; M=M.loc[ids]; y=M.is_diazotroph.astype(int)
original=pd.read_csv(E/'inputs/model_fold_assignments.csv').set_index('assembly_accession')
assert len(X)==426 and X.shape[1]==2286

# Extract model factories and scoring definitions from the independently verified preceding analysis.
tree=ast.parse((E/'analysis/sample_identity.py').read_text())
keep=[n for n in tree.body if (isinstance(n,ast.FunctionDef) and n.name=='metrics') or
      (isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='models' for t in n.targets))]
exec(compile(ast.Module(body=keep,type_ignores=[]),'verified_models_and_metrics','exec'))
models={k:v for k,v in models.items() if k in ['Random Forest','Gradient Boosting','Logistic Regression','XGBoost']}
scenarios={'training_prevalence_at_least40':set()}
manifest=dict(scope='Nested fixed-matrix sensitivity adding at least40 carriers in each training fold to exact train-only purity and >=3 training genera. Original426 rows, labels, folds and model defaults frozen; no retuning or inventory reranking.',
              identifiability='A family globally present in fewer than40 genomes cannot have40 training carriers, so the nested mask is identifiable within the released2286-family matrix. Global clustering/family definitions remain fixed and are not reconstructed.',
              interpretation='Absolute40 in a smaller training cohort is a stricter relative prevalence criterion than40 in426; this is altered eligibility sensitivity, not a replacement primary estimate or isolated proof of selection bias.',
              model_parameters={k:v().get_params() for k,v in models.items()},
              python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,xgboost=xgboost.__version__,threads=1,
              input_sha256={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'gene_family_matrix.csv',P/'complete_genomes_with_proteins.csv']})
(D/'design_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
preds=[];foldresults=[];features=[];assignments=[];counts=[];coefficients=[];summary=[]
for scenario,removed in scenarios.items():
    print('START',scenario,flush=True)
    eligible=[a for a in ids if a not in removed]
    mm=M.loc[eligible]; yy=y.loc[eligible]
    for fold in range(1,6):
        test=[a for a in eligible if original.loc[a,'fold']==fold]; train=[a for a in eligible if a not in set(test)]
        assert not set(mm.loc[test,'genus'])&set(mm.loc[train,'genus'])
        present=X.loc[train].eq(1);n=present.sum();k=present.mul(yy.loc[train],axis=0).sum()
        breadth=present.groupby(mm.loc[train,'genus']).any().sum()
        cols=sorted(n.index[(n>=40)&(10*k>n)&(10*k<9*n)&(breadth>=3)])
        features.extend(dict(scenario=scenario,fold=fold,gene_family=c) for c in cols)
        assignments.extend(dict(scenario=scenario,fold=fold,assembly_accession=a,stored_genus=mm.loc[a,'genus'],audit_group=mm.loc[a,'genus'],y_true=int(yy[a])) for a in test)
        counts.append(dict(scenario=scenario,fold=fold,n_train=len(train),n_test=len(test),n_train_pos=int(yy.loc[train].sum()),n_test_pos=int(yy.loc[test].sum()),n_train_genera=mm.loc[train,'genus'].nunique(),n_test_genera=mm.loc[test,'genus'].nunique(),n_features=len(cols)))
        for model_name,factory in models.items():
            started=time.perf_counter()
            if model_name=='assembly_length_only_LR':z=mm[['total_ungapped_length']].astype(float)
            elif model_name=='retained_family_count_only_LR':z=X.loc[eligible].eq(1).sum(axis=1).to_frame('retained_family_count')
            else:z=X.loc[eligible,cols]
            xt=z.loc[train];xv=z.loc[test];scaler=None
            if model_name=='Logistic Regression' or model_name.endswith('_LR'):
                scaler=StandardScaler();xt=scaler.fit_transform(xt);xv=scaler.transform(xv)
            model=factory()
            with threadpool_limits(limits=1):model.fit(xt,yy.loc[train])
            p=pd.DataFrame(dict(scenario=scenario,model=model_name,fold=fold,assembly_accession=test,genus=mm.loc[test,'genus'].to_numpy(),audit_group=mm.loc[test,'genus'].to_numpy(),y_true=yy.loc[test].to_numpy(),y_prob=model.predict_proba(xv)[:,1],y_pred=model.predict(xv)))
            preds.append(p)
            row=metrics(p);row.update(scenario=scenario,model=model_name,fold=fold,n_train=len(train),n_test=len(test),n_train_pos=int(yy.loc[train].sum()),n_test_pos=int(yy.loc[test].sum()),n_features=z.shape[1],seconds=time.perf_counter()-started);foldresults.append(row)
            if scaler is not None:
                coefficients.extend(dict(scenario=scenario,model=model_name,fold=fold,feature=c,train_mean=float(mu),train_scale=float(sd),coefficient=float(co),intercept=float(model.intercept_[0])) for c,mu,sd,co in zip(z.columns,scaler.mean_,scaler.scale_,model.coef_[0]))
    for model_name in models:
        p=pd.concat([p for p in preds if p.scenario.iloc[0]==scenario and p.model.iloc[0]==model_name])
        f=pd.DataFrame(foldresults).query('scenario==@scenario and model==@model_name')
        row=dict(scenario=scenario,model=model_name,n_genomes=len(eligible),n_positive=int(yy.sum()),n_negative=int((yy==0).sum()),n_stored_genera=mm.genus.nunique(),auc_defined_folds=int(f.roc_auc.notna().sum()))
        for key in metrics(p):row[key+'_mean']=f[key].mean();row[key+'_sd']=f[key].std(ddof=1)
        row.update({'pooled_'+key:v for key,v in metrics(p).items()});summary.append(row)
    pd.concat(preds).to_csv(D/'oof_predictions.csv',index=False)
    pd.DataFrame(foldresults).to_csv(D/'fold_metrics.csv',index=False)
    pd.DataFrame(features).to_csv(D/'fold_family_features.csv',index=False)
    pd.DataFrame(assignments).to_csv(D/'fold_assignments.csv',index=False)
    pd.DataFrame(counts).to_csv(D/'fold_cohort_and_mask_counts.csv',index=False)
    pd.DataFrame(coefficients).to_csv(D/'lr_scalers_and_coefficients.csv',index=False)
    pd.DataFrame(summary).to_csv(D/'summary.csv',index=False)
    print(pd.DataFrame(summary).query('scenario==@scenario')[['scenario','model','n_genomes','n_positive','roc_auc_mean','roc_auc_sd','genus_balanced_auc_mean','pooled_genus_balanced_auc']].to_string(index=False),flush=True)

# Separate evaluation-population change from refitting using original predictions on identical surviving rows.
original_oof=pd.read_csv(E/'reference/sample_identity/oof_predictions.csv').query("scenario=='original_control'")
base=pd.read_csv(E/'reference/sample_identity/summary.csv').query("scenario=='original_control'").set_index('model')
references=[];deltas=[]
for s in summary:
    p=original_oof[(original_oof.model==s['model'])&~original_oof.assembly_accession.isin(scenarios[s['scenario']])]
    f=pd.DataFrame([metrics(q) for _,q in p.groupby('fold')])
    ref=dict(scenario=s['scenario'],model=s['model'],reference='original_OOF_scores_restricted_to_same_surviving_rows')
    for key in metrics(p):ref[key+'_mean']=f[key].mean()
    ref.update({'pooled_'+key:v for key,v in metrics(p).items()});references.append(ref)
    d=dict(scenario=s['scenario'],model=s['model'])
    for key in ['roc_auc_mean','genus_balanced_auc_mean','average_precision_mean','pooled_roc_auc','pooled_genus_balanced_auc']:
        d[key+'_minus_original426']=s[key]-base.loc[s['model'],key]
        d[key+'_minus_same_rows_original_OOF']=s[key]-ref[key]
    deltas.append(d)
pd.DataFrame(references).to_csv(D/'original_OOF_same_rows_reference.csv',index=False)
pd.DataFrame(deltas).to_csv(D/'sensitivity_differences.csv',index=False)
manifest['completed_fits']=len(foldresults)
(D/'design_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print('DONE: 20 fits; training prevalence>=40 only; canonical data/resources untouched.',flush=True)
