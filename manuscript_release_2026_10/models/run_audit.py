#!/usr/bin/env python
"""Portable fixed-matrix audit. All paths are relative to this bundle or --out."""
from pathlib import Path
import argparse,ast,contextlib,hashlib,io,json,os,sys,subprocess,importlib.metadata
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import numpy as np
import pandas as pd
import scipy,sklearn,xgboost,threadpoolctl
from fractions import Fraction
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score,precision_score,recall_score,f1_score,roc_auc_score,average_precision_score,precision_recall_curve,auc
from threadpoolctl import threadpool_limits
BASE=Path(__file__).resolve().parent

def source_functions(path,names):
    tree=ast.parse(path.read_text(encoding='utf8'))
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
    assert len(nodes)==len(names)
    ns=dict(np=np,pd=pd,Fraction=Fraction,StandardScaler=StandardScaler,
            accuracy_score=accuracy_score,precision_score=precision_score,recall_score=recall_score,f1_score=f1_score,roc_auc_score=roc_auc_score)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)
    return ns

def load():
    integrity=json.loads((BASE/'input_manifest.json').read_text())
    for name,expected in integrity['sha256'].items():
        assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==expected, name
    X=pd.read_csv(BASE/'inputs/gene_family_matrix.csv',index_col=0)
    m=pd.read_csv(BASE/'inputs/complete_genomes_with_proteins.csv').drop_duplicates('assembly_accession').set_index('assembly_accession')
    common=sorted(set(X.index)&set(m.index));X=X.loc[common];m=m.loc[common]
    assert X.shape==(426,2286) and m.genus.nunique()==75
    quality=pd.read_csv(BASE/'inputs/assembly_quality.tsv',sep='\t').drop_duplicates('accession').set_index('accession')
    np.testing.assert_array_equal(m.total_ungapped_length,quality.loc[m.index,'total_ungapped_length'])
    y=m.is_diazotroph.astype(int);assert y.sum()==112
    names=['_safe_select','family_purity_statistics','filter_min_genera','select_training_features','genus_cross_validation_split','evaluate_model','run_genus_cv']
    new=source_functions(BASE/'source/04_classify_proposed.py',names)
    old=source_functions(BASE/'source/04_classify_original.py',['genus_cross_validation_split'])
    for seed in range(42,52):assert old['genus_cross_validation_split'](m,5,seed)==new['genus_cross_validation_split'](m,5,seed)
    return X,m,y,new

def metrics(p):
    yy=p.y_true.to_numpy();pr=p.y_prob.to_numpy();hh=p.y_pred.to_numpy()
    w=1/p.groupby('genus').genus.transform('size')
    prec,rec,_=precision_recall_curve(yy,pr)
    return dict(roc_auc=roc_auc_score(yy,pr),average_precision=average_precision_score(yy,pr),
                genus_balanced_auc=roc_auc_score(yy,pr,sample_weight=w),genus_balanced_average_precision=average_precision_score(yy,pr,sample_weight=w),
                pr_auc_trapezoid=auc(rec,prec),accuracy=accuracy_score(yy,hh),precision=precision_score(yy,hh,zero_division=0),recall=recall_score(yy,hh,zero_division=0),f1=f1_score(yy,hh,zero_division=0))

def save(out,folds,preds,features):
    out.mkdir(parents=True,exist_ok=True)
    ff=pd.DataFrame(folds);pp=pd.concat(preds,ignore_index=True)
    ff.to_csv(out/'fold_metrics.csv',index=False);pp.to_csv(out/'oof_predictions.csv',index=False)
    pd.DataFrame(features).to_csv(out/'fold_features.csv',index=False)
    summaries=[]
    for (seed,model),p in pp.groupby(['seed','model'],sort=False):
        f=ff[(ff.seed==seed)&(ff.model==model)];row=dict(seed=seed,model=model)
        for key in metrics(p):row[key+'_mean']=f[key].mean();row[key+'_sd']=f[key].std(ddof=1)
        row.update({'pooled_'+key:value for key,value in metrics(p).items()})
        row.update(n_features_min=f.n_features.min(),n_features_max=f.n_features.max());summaries.append(row)
    pd.DataFrame(summaries).to_csv(out/'summary.csv',index=False)

def seed42(out,X,m,y,new):
    stats=new['family_purity_statistics'](X,y).set_index('gene_family')
    eligible=new['filter_min_genera'](X[stats.index[~stats.is_pure]],m).columns.tolist()
    assert len(eligible)==624
    models={'Random Forest':RandomForestClassifier(n_estimators=100,max_depth=10,min_samples_leaf=5,random_state=42,n_jobs=1),
            'Gradient Boosting':GradientBoostingClassifier(n_estimators=100,max_depth=5,learning_rate=0.1,random_state=42),
            'Logistic Regression':LogisticRegression(max_iter=1000,random_state=42),
            'XGBoost':xgboost.XGBClassifier(n_estimators=200,max_depth=5,learning_rate=0.1,subsample=0.8,colsample_bytree=0.8,random_state=42,n_jobs=1,use_label_encoder=False,eval_metric='logloss')}
    folds=[];preds=[];features=[];out.mkdir(parents=True,exist_ok=True)
    for name,model in models.items():
        log=io.StringIO()
        with threadpool_limits(limits=1),contextlib.redirect_stdout(log):
            result=new['run_genus_cv'](X,y,m,model,name,5,42,inventory_features=eligible)
        (out/(name.replace(' ','_')+'.log')).write_text(log.getvalue(),encoding='utf8')
        p=result['oof_predictions'].copy();p['seed']=42;p['model']=name;preds.append(p)
        for _,source in result['fold_metrics'].iterrows():
            b=p[p.fold==source.fold];row=metrics(b)
            row.update(seed=42,model=name,fold=int(source.fold),n_train=int(source.n_train),n_test=int(source.n_test),
                       n_train_pos=int(source.n_train_positive),n_test_pos=int(source.n_test_positive),n_train_genera=60,n_test_genera=15,genus_overlap=0,n_features=int(source.n_features));folds.append(row)
        features.extend(dict(seed=42,model=name,fold=int(r.fold),gene_family=r.gene_family) for _,r in result['fold_features'].iterrows())
        result['feature_importance'].to_csv(out/f"AUDIT_ONLY_cv_importance_{name.replace(' ','_')}.csv",index=False)
        save(out,folds,preds,features);print('seed42 completed:',name,flush=True)
    stats.to_csv(out/'exact_global_purity_descriptive.csv')

def repeated(out,X,m,y,new):
    arms=['corrected_family_LR','assembly_length_only_LR','retained_family_count_only_LR','joint_family_plus_length_LR']
    count=X.eq(1).sum(axis=1);length=m.total_ungapped_length.astype(float)
    folds=[];preds=[];features=[]
    for seed in range(42,52):
        for fold,test in enumerate(new['genus_cross_validation_split'](m,5,seed),1):
            train=[a for a in X.index if a not in set(test)]
            assert not set(m.loc[train,'genus'])&set(m.loc[test,'genus'])
            present=X.loc[train].eq(1);n=present.sum();k=present.mul(y.loc[train],axis=0).sum();breadth=present.groupby(m.loc[train,'genus']).any().sum()
            cols=sorted(n.index[(n>0)&(10*k>n)&(10*k<9*n)&(breadth>=3)])
            if seed==42:assert cols==new['select_training_features'](X.loc[train],y.loc[train],m.loc[train])
            features.extend(dict(seed=seed,fold=fold,gene_family=c) for c in cols)
            for arm in arms:
                if arm=='corrected_family_LR':z=X[cols]
                elif arm=='assembly_length_only_LR':z=length.to_frame('assembly_length_bp')
                elif arm=='retained_family_count_only_LR':z=count.to_frame('retained_family_count')
                else:z=X[cols].assign(__assembly_length_bp=length)
                scaler=StandardScaler();xt=scaler.fit_transform(z.loc[train]);xv=scaler.transform(z.loc[test])
                model=LogisticRegression(max_iter=1000,random_state=42)
                with threadpool_limits(limits=1):model.fit(xt,y.loc[train])
                p=pd.DataFrame(dict(seed=seed,fold=fold,model=arm,assembly_accession=test,genus=m.loc[test,'genus'].to_numpy(),y_true=y.loc[test].to_numpy(),y_prob=model.predict_proba(xv)[:,1],y_pred=model.predict(xv)))
                preds.append(p);row=metrics(p);row.update(seed=seed,fold=fold,model=arm,n_train=len(train),n_test=len(test),n_train_pos=int(y.loc[train].sum()),n_test_pos=int(y.loc[test].sum()),n_train_genera=60,n_test_genera=15,genus_overlap=0,n_features=len(z.columns));folds.append(row)
        save(out,folds,preds,features);print('repeated completed seed:',seed,flush=True)
    ss=pd.read_csv(out/'summary.csv');pairs=[]
    for seed in range(42,52):
        for a,b,label in [('corrected_family_LR','assembly_length_only_LR','family_minus_length'),('corrected_family_LR','retained_family_count_only_LR','family_minus_count'),('joint_family_plus_length_LR','assembly_length_only_LR','joint_minus_length'),('joint_family_plus_length_LR','corrected_family_LR','joint_minus_family')]:
            aa=ss[(ss.seed==seed)&(ss.model==a)].iloc[0];bb=ss[(ss.seed==seed)&(ss.model==b)].iloc[0]
            row=dict(seed=seed,comparison=label)
            for key in ['roc_auc_mean','genus_balanced_auc_mean','pooled_genus_balanced_auc','average_precision_mean','pooled_average_precision']:row[key+'_delta']=aa[key]-bb[key]
            pairs.append(row)
    pd.DataFrame(pairs).to_csv(out/'paired_seed_differences.csv',index=False)

def independent_auc(y,p,w):
    y=np.asarray(y);p=np.asarray(p);w=np.asarray(w);a=p[y==1];b=p[y==0];wa=w[y==1];wb=w[y==0]
    comparison=(a[:,None]>b[None,:])+0.5*(a[:,None]==b[None,:])
    return np.sum(comparison*wa[:,None]*wb[None,:])/(wa.sum()*wb.sum())

def independent_ap(y,p,w):
    y=np.asarray(y);p=np.asarray(p);w=np.asarray(w);total=w[y==1].sum();prev=0.;answer=0.
    for t in sorted(set(p),reverse=True):
        chosen=p>=t;tp=w[chosen&(y==1)].sum();rec=tp/total
        answer+=(tp/w[chosen].sum())*(rec-prev);prev=rec
    return answer

def verify(out,reference,X,m,y,new):
    p=pd.read_csv(out/'oof_predictions.csv');f=pd.read_csv(out/'fold_metrics.csv');s=pd.read_csv(out/'summary.csv');checks=[]
    for (seed,model,fold),b in p.groupby(['seed','model','fold']):
        expected=new['genus_cross_validation_split'](m,5,int(seed))[int(fold)-1]
        assert b.assembly_accession.tolist()==expected and not b.assembly_accession.duplicated().any()
        np.testing.assert_array_equal(b.y_true,y.loc[expected]);np.testing.assert_array_equal(b.genus,m.loc[expected,'genus'])
        saved=f[(f.seed==seed)&(f.model==model)&(f.fold==fold)].iloc[0]
        w=1/b.groupby('genus').genus.transform('size')
        for metric,value in {'roc_auc':independent_auc(b.y_true,b.y_prob,np.ones(len(b))), 'genus_balanced_auc':independent_auc(b.y_true,b.y_prob,w),
                             'average_precision':independent_ap(b.y_true,b.y_prob,np.ones(len(b))), 'genus_balanced_average_precision':independent_ap(b.y_true,b.y_prob,w)}.items():
            assert abs(value-saved[metric])<1e-12;checks.append(dict(seed=seed,model=model,fold=fold,metric=metric,delta=value-saved[metric]))
    for (seed,model),b in p.groupby(['seed','model']):
        assert len(b)==b.assembly_accession.nunique()==426
        saved=s[(s.seed==seed)&(s.model==model)].iloc[0];w=1/b.groupby('genus').genus.transform('size')
        for metric,value in {'roc_auc':independent_auc(b.y_true,b.y_prob,np.ones(len(b))), 'genus_balanced_auc':independent_auc(b.y_true,b.y_prob,w),
                             'average_precision':independent_ap(b.y_true,b.y_prob,np.ones(len(b))), 'genus_balanced_average_precision':independent_ap(b.y_true,b.y_prob,w)}.items():
            assert abs(value-saved['pooled_'+metric])<1e-12;checks.append(dict(seed=seed,model=model,fold='pooled',metric=metric,delta=value-saved['pooled_'+metric]))
        ff=f[(f.seed==seed)&(f.model==model)]
        for key in ['roc_auc','genus_balanced_auc','average_precision','genus_balanced_average_precision']:
            assert abs(ff[key].mean()-saved[key+'_mean'])<1e-12
            assert abs(ff[key].std(ddof=1)-saved[key+'_sd'])<1e-12
    ref=pd.read_csv(reference/'oof_predictions.csv')
    if 'seed' not in ref:ref['seed']=42
    keys=['seed','model','fold','assembly_accession'];merged=p.merge(ref,on=keys,suffixes=('_new','_ref'),validate='one_to_one')
    assert len(merged)==len(p)==len(ref)
    np.testing.assert_allclose(merged.y_prob_new,merged.y_prob_ref,rtol=1e-10,atol=1e-12)
    np.testing.assert_array_equal(merged.y_true_new,merged.y_true_ref)
    np.testing.assert_array_equal(merged.y_pred_new,merged.y_pred_ref)
    pd.DataFrame(checks).to_csv(out/'independent_metric_checks.csv',index=False)
    result=dict(status='PASS',independent_metric_checks=len(checks),reference_prediction_rows=len(merged),
                maximum_reference_probability_delta=float((merged.y_prob_new-merged.y_prob_ref).abs().max()),maximum_independent_metric_delta=float(pd.DataFrame(checks).delta.abs().max()))
    (out/'verification.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(out.name,result,flush=True)

def contracts(X,m,y,new):
    toy_y=pd.Series([1]*10+[0]*90);toy=pd.DataFrame({'exact10':np.ones(100,dtype=int)})
    assert bool(new['family_purity_statistics'](toy,toy_y).is_pure.iloc[0])
    assert bool(new['family_purity_statistics'](toy,1-toy_y).is_pure.iloc[0])
    test=new['genus_cross_validation_split'](m,5,42)[0];train=[a for a in X.index if a not in set(test)]
    flipped=y.copy();flipped.loc[test]=1-flipped.loc[test]
    assert new['select_training_features'](X.loc[train],y.loc[train],m.loc[train])==new['select_training_features'](X.loc[train],flipped.loc[train],m.loc[train])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['primary','seed42','repeated','sample_identity','current_status','official_taxonomy','sensitivities','all','verify'],default='all')
    parser.add_argument('--out',type=Path,default=BASE/'rerun_outputs')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if args.mode=='primary':args.mode='seed42'
    X,m,y,new=load();contracts(X,m,y,new)
    runtime=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,xgboost=xgboost.__version__,threadpoolctl=threadpoolctl.__version__,threads=1,mode=args.mode)
    runtime['distribution_metadata']={name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','scikit-learn','xgboost','threadpoolctl']}
    (args.out/'runtime.json').write_text(json.dumps(runtime,indent=2),encoding='utf8')
    if args.mode in ['seed42','all']:seed42(args.out/'seed42',X,m,y,new)
    if args.mode in ['repeated','all']:repeated(args.out/'repeated',X,m,y,new)
    if args.mode=='verify' and not any((args.out/mode/'oof_predictions.csv').exists() for mode in ['seed42','repeated','sample_identity','current_status','official_taxonomy']):
        parser.error('No saved predictions under --out; run seed42/repeated or use --out reference.')
    for mode in ['seed42','repeated']:
        if args.mode in [mode,'all','verify'] and (args.out/mode/'oof_predictions.csv').exists():verify(args.out/mode,BASE/'reference'/mode,X,m,y,new)
    if args.mode in ['sample_identity','current_status','official_taxonomy','sensitivities','all','verify']:
        subprocess.run([sys.executable,str(BASE/'run_sensitivities.py'),'--mode',('all' if args.mode in ['sensitivities','all'] else args.mode),'--out',str(args.out.resolve())],check=True)
    print('Complete. Released inventory/ranking files were not read as model inputs or modified.',flush=True)
if __name__=='__main__':main()
