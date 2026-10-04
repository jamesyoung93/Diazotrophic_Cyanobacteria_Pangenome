"""Bounded sample-identity/grouping sensitivity; fixed matrix, no canonical edits.

Original control replay; remove either flagged row while retaining original folds;
separately connect explicitly named genus labels and repartition with seed42.
"""
from pathlib import Path
import os,ast,json,sys,hashlib,time
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import numpy as np
import pandas as pd
import sklearn,scipy,xgboost
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score,precision_score,recall_score,f1_score,roc_auc_score,average_precision_score
from threadpoolctl import threadpool_limits
OUT=Path(__file__).resolve().parents[1];D=Path(os.environ['MODEL_AUDIT_OUTPUT'])/'sample_identity';D.mkdir(parents=True,exist_ok=True)
REPO=OUT;RUN=OUT/'inputs'
CODE=OUT/'source/04_classify_original.py'
X=pd.read_csv(RUN/'gene_family_matrix.csv',index_col=0)
meta=pd.read_csv(RUN/'complete_genomes_with_proteins.csv').drop_duplicates('assembly_accession').set_index('assembly_accession')
common=sorted(set(X.index)&set(meta.index));X=X.loc[common];meta=meta.loc[common];y=meta.is_diazotroph.astype(int)
oldassign=pd.read_csv(OUT/'inputs/model_fold_assignments.csv').set_index('assembly_accession')
node=next(n for n in ast.parse(CODE.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='genus_cross_validation_split')
ns={'np':np};exec(compile(ast.Module(body=[node],type_ignores=[]),str(CODE),'exec'),ns);split=ns['genus_cross_validation_split']
original_folds=split(meta,5,42)
for i,test in enumerate(original_folds,1):assert test==oldassign.index[oldassign.fold==i].tolist()
PAIR=['GCF_047302775.1','GCF_047302855.1']
assert (X.loc[PAIR[0]]==X.loc[PAIR[1]]).all() and X.loc[PAIR].sum(axis=1).eq(178).all()
assert y.loc[PAIR].eq(0).all()
pairsets={a:set() for a in PAIR}
for chunk in pd.read_csv(RUN/'genome_protein_family_map.tsv',sep='\t',chunksize=100000):
    for a in PAIR:
        sub=chunk[chunk.genome_accession==a];pairsets[a].update(zip(sub.gene_family,sub.protein_accession))
assert pairsets[PAIR[0]]==pairsets[PAIR[1]] and len(pairsets[PAIR[0]])==179
for field in ['total_ungapped_length','gc_percent','number_of_contigs','contig_n50','best_hit_ids','nifH_best_evalue','nifD_best_evalue','nifK_best_evalue']:
    assert meta.loc[PAIR[0],field]==meta.loc[PAIR[1],field]
meta.loc[PAIR,['organism_full','genus','is_diazotroph','total_ungapped_length','gc_percent','number_of_contigs','contig_n50','best_hit_ids']].to_csv(D/'flagged_pair_evidence.csv')

# Connected components only for the specifically requested sensitivity edges.
edges=[('Synechococcus','Prochlorococcus'),('Synechococcus','[Synechococcus]'),('Leptolyngbya','[Leptolyngbya]'),('Phormidium','[Phormidium]')]
labels=sorted(meta.genus.unique());parent={g:g for g in labels}
def find(a):
    while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
    return a
for a,b in edges:assert a in parent and b in parent;parent[find(b)]=find(a)
components={}
for g in labels:components.setdefault(find(g),[]).append(g)
mapping={g:(' + '.join(sorted(v)) if len(v)>1 else g) for v in components.values() for g in v}
assert len(set(mapping.values()))==71
pd.DataFrame([dict(stored_genus=g,audit_group=mapping[g]) for g in labels]).to_csv(D/'stored_genus_to_audit_group.csv',index=False)
models={
    'Random Forest':lambda:RandomForestClassifier(n_estimators=100,max_depth=10,min_samples_leaf=5,random_state=42,n_jobs=1),
    'Gradient Boosting':lambda:GradientBoostingClassifier(n_estimators=100,max_depth=5,learning_rate=0.1,random_state=42),
    'Logistic Regression':lambda:LogisticRegression(max_iter=1000,random_state=42),
    'XGBoost':lambda:xgboost.XGBClassifier(n_estimators=200,max_depth=5,learning_rate=0.1,subsample=0.8,colsample_bytree=0.8,random_state=42,n_jobs=1,use_label_encoder=False,eval_metric='logloss'),
    'assembly_length_only_LR':lambda:LogisticRegression(max_iter=1000,random_state=42),
    'retained_family_count_only_LR':lambda:LogisticRegression(max_iter=1000,random_state=42),
}
scenarios=[('original_control',None,False),('remove_MIT_S9506',PAIR[0],False),('remove_SS51',PAIR[1],False),('connected_group_repartition',None,True)]
manifest={'scenarios':[dict(name=s,removed_accession=drop,grouping_changed=merge) for s,drop,merge in scenarios],
          'scope':'Fixed original2286-family matrix. Single-row exclusions preserve original saved fold IDs. Connected-group sensitivity recomputes five round-robin folds with seed42 and71 groups. Not whole-genome deduplication or curated taxonomy.',
          'grouping_edges':edges,'training_filter':'Exact10%/90% purity and>=3 training audit groups. Audit groups equal stored genus labels except connected_group_repartition.',
          'metric_weights':'genus_balanced uses unchanged stored-genus tokens for comparability; audit_group_balanced uses connected group units. Both pooled and within-fold forms are exported.',
          'model_parameters':{name:factory().get_params() for name,factory in models.items()},
          'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'xgboost':xgboost.__version__,'threads':1,
          'input_sha256':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [RUN/'gene_family_matrix.csv',RUN/'complete_genomes_with_proteins.csv',RUN/'genome_protein_family_map.tsv',CODE]}}
(D/'design_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')

def metrics(p):
    yy=p.y_true.to_numpy();pr=p.y_prob.to_numpy();hh=p.y_pred.to_numpy()
    both=len(np.unique(yy))==2
    result=dict(accuracy=accuracy_score(yy,hh),precision=precision_score(yy,hh,zero_division=0),recall=recall_score(yy,hh,zero_division=0),f1=f1_score(yy,hh,zero_division=0),
                roc_auc=roc_auc_score(yy,pr) if both else np.nan,average_precision=average_precision_score(yy,pr) if both else np.nan)
    for name,column in [('genus','genus'),('audit_group','audit_group')]:
        weights=1/p.groupby(column)[column].transform('size')
        result[name+'_balanced_auc']=roc_auc_score(yy,pr,sample_weight=weights) if both else np.nan
        result[name+'_balanced_average_precision']=average_precision_score(yy,pr,sample_weight=weights) if both else np.nan
    return result

allfolds=[];allpred=[];features=[];assignments=[];maskcounts=[];lrparams=[];summaries=[]
for scenario,remove,merge in scenarios:
    print('START',scenario,flush=True)
    ids=[a for a in X.index if a!=remove];xx=X.loc[ids];mm=meta.loc[ids].copy();yy=y.loc[ids]
    mm['stored_genus']=mm.genus;mm['audit_group']=mm.genus.map(mapping) if merge else mm.genus
    splitmeta=mm.copy();splitmeta['genus']=mm.audit_group
    folds=split(splitmeta,5,42) if merge else [[a for a in test if a!=remove] for test in original_folds]
    assert sum(map(len,folds))==len(ids)==len(set(a for test in folds for a in test))
    for fold,test in enumerate(folds,1):
        train=[a for a in ids if a not in set(test)]
        assert not set(mm.loc[train,'audit_group'])&set(mm.loc[test,'audit_group'])
        present=xx.loc[train].eq(1);n=present.sum();k=present.mul(yy.loc[train],axis=0).sum()
        purity=(n>0)&(10*k>n)&(10*k<9*n)
        audit_breadth=present.groupby(mm.loc[train,'audit_group']).any().sum()
        stored_breadth=present.groupby(mm.loc[train,'stored_genus']).any().sum()
        cols=sorted(n.index[purity&(audit_breadth>=3)])
        stored_cols=sorted(n.index[purity&(stored_breadth>=3)])
        maskcounts.append(dict(scenario=scenario,fold=fold,n_train=len(train),n_test=len(test),n_train_pos=int(yy.loc[train].sum()),n_test_pos=int(yy.loc[test].sum()),
                              n_test_stored_genera=mm.loc[test,'stored_genus'].nunique(),n_test_audit_groups=mm.loc[test,'audit_group'].nunique(),n_features=len(cols),
                              n_features_if_stored_genus_breadth=len(stored_cols),removed_by_conservative_breadth=';'.join(sorted(set(stored_cols)-set(cols)))))
        features.extend(dict(scenario=scenario,fold=fold,gene_family=c) for c in cols)
        assignments.extend(dict(scenario=scenario,fold=fold,assembly_accession=a,stored_genus=mm.loc[a,'stored_genus'],audit_group=mm.loc[a,'audit_group'],y_true=int(yy.loc[a])) for a in test)
        for name,factory in models.items():
            start=time.perf_counter()
            if name=='assembly_length_only_LR':z=mm[['total_ungapped_length']].astype(float)
            elif name=='retained_family_count_only_LR':z=xx.eq(1).sum(axis=1).to_frame('retained_family_count')
            else:z=xx[cols]
            xt,xv=z.loc[train],z.loc[test]
            scaler=None
            if name=='Logistic Regression' or name.endswith('_LR'):
                scaler=StandardScaler();xt=scaler.fit_transform(xt);xv=scaler.transform(xv)
            model=factory()
            with threadpool_limits(limits=1):model.fit(xt,yy.loc[train])
            p=pd.DataFrame(dict(scenario=scenario,model=name,fold=fold,assembly_accession=test,genus=mm.loc[test,'stored_genus'].to_numpy(),audit_group=mm.loc[test,'audit_group'].to_numpy(),y_true=yy.loc[test].to_numpy(),y_prob=model.predict_proba(xv)[:,1],y_pred=model.predict(xv)))
            row=metrics(p);row.update(scenario=scenario,model=name,fold=fold,n_train=len(train),n_test=len(test),n_train_pos=int(yy.loc[train].sum()),n_test_pos=int(yy.loc[test].sum()),n_features=z.shape[1],seconds=time.perf_counter()-start)
            allfolds.append(row);allpred.append(p)
            if scaler is not None:
                lrparams.extend(dict(scenario=scenario,model=name,fold=fold,feature=c,train_mean=float(mu),train_scale=float(sd),coefficient=float(co),intercept=float(model.intercept_[0])) for c,mu,sd,co in zip(z.columns,scaler.mean_,scaler.scale_,model.coef_[0]))
    for name in models:
        p=pd.concat([b for b in allpred if b.scenario.iloc[0]==scenario and b.model.iloc[0]==name])
        ff=pd.DataFrame(allfolds).query('scenario==@scenario and model==@name')
        entry=dict(scenario=scenario,model=name,n_genomes=len(ids),n_stored_genera=mm.stored_genus.nunique(),n_audit_groups=mm.audit_group.nunique(),auc_defined_folds=int(ff.roc_auc.notna().sum()))
        for key in metrics(p):entry[key+'_mean']=ff[key].mean();entry[key+'_sd']=ff[key].std(ddof=1)
        entry.update({'pooled_'+key:value for key,value in metrics(p).items()});summaries.append(entry)
    pd.DataFrame(allfolds).to_csv(D/'fold_metrics.csv',index=False);pd.concat(allpred).to_csv(D/'oof_predictions.csv',index=False)
    pd.DataFrame(features).to_csv(D/'fold_family_features.csv',index=False);pd.DataFrame(assignments).to_csv(D/'fold_assignments.csv',index=False)
    pd.DataFrame(maskcounts).to_csv(D/'fold_cohort_and_mask_counts.csv',index=False);pd.DataFrame(lrparams).to_csv(D/'lr_scalers_and_coefficients.csv',index=False)
    pd.DataFrame(summaries).to_csv(D/'summary.csv',index=False)
    print(pd.DataFrame(summaries).query('scenario==@scenario')[['scenario','model','roc_auc_mean','roc_auc_sd','genus_balanced_auc_mean','pooled_genus_balanced_auc','pooled_audit_group_balanced_auc']].to_string(index=False),flush=True)

# Check original-control predictions exactly link to the already verified analyses.
allpred=pd.concat(allpred);reference=pd.read_csv(OUT/'reference/seed42/oof_predictions.csv')
size_ref=pd.read_csv(OUT/'reference/repeated/oof_predictions.csv')
size_ref=size_ref[(size_ref.seed==42)&size_ref.model.isin(['assembly_length_only_LR','retained_family_count_only_LR'])].copy()
reference=pd.concat([reference,size_ref],ignore_index=True)
control=allpred[allpred.scenario=='original_control']
compare=control.merge(reference,on=['model','fold','assembly_accession'],suffixes=('_new','_ref'),validate='one_to_one')
assert len(compare)==2556
np.testing.assert_allclose(compare.y_prob_new,compare.y_prob_ref,rtol=1e-12,atol=1e-12)
manifest['original_control_max_prediction_difference']=float((compare.y_prob_new-compare.y_prob_ref).abs().max())
manifest['completed_fits']=len(allfolds)
(D/'design_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')

# Isolate test-population removal from changes due to refitting/filtering.
restricted=[]
for scenario,remove,merge in scenarios:
    if remove is None:continue
    for name in models:
        p=control[(control.model==name)&(control.assembly_accession!=remove)]
        foldresults=pd.DataFrame([metrics(b) for _,b in p.groupby('fold')])
        row=dict(scenario=scenario,model=name,reference='original_OOF_probabilities_restricted_to_surviving_rows')
        for key in metrics(p):row[key+'_mean']=foldresults[key].mean()
        row.update({'pooled_'+key:value for key,value in metrics(p).items()});restricted.append(row)
pd.DataFrame(restricted).to_csv(D/'original_OOF_same_rows_reference.csv',index=False)
base=pd.DataFrame(summaries).query("scenario=='original_control'").set_index('model')
delta=[]
for _,row in pd.DataFrame(summaries).iterrows():
    if row.scenario=='original_control':continue
    item=dict(scenario=row.scenario,model=row.model)
    for key in ['roc_auc_mean','genus_balanced_auc_mean','average_precision_mean','pooled_roc_auc','pooled_genus_balanced_auc','pooled_audit_group_balanced_auc']:
        item[key+'_minus_original']=row[key]-base.loc[row.model,key]
        if row.scenario.startswith('remove_'):
            same=next(r for r in restricted if r['scenario']==row.scenario and r['model']==row.model)
            item[key+'_minus_same_rows_reference']=row[key]-same[key]
    delta.append(item)
pd.DataFrame(delta).to_csv(D/'sensitivity_differences.csv',index=False)
print('DONE:120 fits, unchanged original-control predictions, no canonical data/resource changes.',flush=True)
