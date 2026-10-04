"""Read-only comparison of saved model outputs across numerical backends.

The original strict probability tolerance is retained and its failure is reported.
This diagnostic does not relabel a failed reproduction as a successful strict one.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference',type=Path,required=True)
parser.add_argument('--clean',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
parser.add_argument('--assert-strict',action='store_true',help='After writing diagnostics, enforce the unchanged original strict probability assertion (expected to fail for the saved clean run).')
args=parser.parse_args();D=args.out;D.mkdir(parents=True,exist_ok=True)
ref=pd.read_csv(args.reference/'oof_predictions.csv');clean=pd.read_csv(args.clean/'oof_predictions.csv')
keys=['seed','model','fold','assembly_accession'];merged=clean.merge(ref,on=keys,suffixes=('_clean','_reference'),validate='one_to_one')
assert len(merged)==len(ref)==len(clean)==1704
for column in ['genus','y_true']:np.testing.assert_array_equal(merged[column+'_clean'],merged[column+'_reference'])
assert np.isfinite(merged.y_prob_clean).all() and merged.y_prob_clean.between(0,1).all()
merged['absolute_probability_difference']=(merged.y_prob_clean-merged.y_prob_reference).abs()
merged['original_strict_allowance']=1e-12+1e-10*merged.y_prob_reference.abs()
merged['fails_original_strict_tolerance']=merged.absolute_probability_difference>merged.original_strict_allowance
merged['decision_changed']=merged.y_pred_clean!=merged.y_pred_reference
merged['clean_margin_from_probability_half']=(merged.y_prob_clean-0.5).abs()
merged.to_csv(D/'per_prediction_comparison.csv',index=False)
rows=[]
for model,p in merged.groupby('model'):
    rows.append(dict(model=model,n=len(p),n_exactly_equal=int(p.absolute_probability_difference.eq(0).sum()),strict_failures=int(p.fails_original_strict_tolerance.sum()),max_abs_probability_difference=p.absolute_probability_difference.max(),max_difference_over_original_allowance=(p.absolute_probability_difference/p.original_strict_allowance).max(),decision_changes=int(p.decision_changed.sum()),minimum_clean_margin_from_half=p.clean_margin_from_probability_half.min()))
by_model=pd.DataFrame(rows);by_model.to_csv(D/'prediction_comparison_by_model.csv',index=False)
feature_keys=['seed','model','fold','gene_family']
fa=pd.read_csv(args.reference/'fold_features.csv')[feature_keys].sort_values(feature_keys).reset_index(drop=True)
fb=pd.read_csv(args.clean/'fold_features.csv')[feature_keys].sort_values(feature_keys).reset_index(drop=True)
assert not fa.duplicated().any() and not fb.duplicated().any()
pd.testing.assert_frame_equal(fa,fb,check_exact=True)
pd.testing.assert_frame_equal(pd.read_csv(args.reference/'exact_global_purity_descriptive.csv'),pd.read_csv(args.clean/'exact_global_purity_descriptive.csv'))

# Compare every pairwise score ordering, including exact ties, within folds and pooled.
orders=[];changed_pairs=[]
for model,p in merged.groupby('model'):
    for fold,q in [('pooled',p)]+list(p.groupby('fold')):
        a=q.y_prob_reference.to_numpy();b=q.y_prob_clean.to_numpy();ii=np.triu_indices(len(q),1)
        old=np.sign(a[:,None]-a[None,:])[ii];new=np.sign(b[:,None]-b[None,:])[ii]
        for i,j in zip(ii[0][old!=new],ii[1][old!=new]):
            aa=q.iloc[i];bb=q.iloc[j]
            changed_pairs.append(dict(model=model,scope_fold=fold,accession_a=aa.assembly_accession,accession_b=bb.assembly_accession,y_a=int(aa.y_true_reference),y_b=int(bb.y_true_reference),reference_probability_a=aa.y_prob_reference,reference_probability_b=bb.y_prob_reference,clean_probability_a=aa.y_prob_clean,clean_probability_b=bb.y_prob_clean))
        orders.append(dict(model=model,fold=fold,n=len(q),n_unordered_pairs=len(old),changed_pair_order_or_tie=int((old!=new).sum()),reference_tied_pairs=int((old==0).sum()),clean_tied_pairs=int((new==0).sum())))
orders=pd.DataFrame(orders);orders.to_csv(D/'score_ordering_comparison.csv',index=False)
pd.DataFrame(changed_pairs).to_csv(D/'changed_score_order_pairs.csv',index=False)

def auc_pairwise(y,p,w):
    y=np.asarray(y);p=np.asarray(p);w=np.asarray(w);pos=p[y==1];neg=p[y==0];wp=w[y==1];wn=w[y==0]
    return (((pos[:,None]>neg[None,:])+0.5*(pos[:,None]==neg[None,:]))*wp[:,None]*wn[None,:]).sum()/(wp.sum()*wn.sum())
def ap_threshold(y,p,w):
    y=np.asarray(y);p=np.asarray(p);w=np.asarray(w);total=w[y==1].sum();prev=0.;area=0.
    for t in sorted(set(p),reverse=True):
        chosen=p>=t;tp=w[chosen&(y==1)].sum();rec=tp/total
        area+=(tp/w[chosen].sum())*(rec-prev);prev=rec
    return area
def independent_metrics(p):
    y=p.y_true.to_numpy();score=p.y_prob.to_numpy();h=p.y_pred.to_numpy();w=1/p.groupby('genus').genus.transform('size')
    tp=int(((y==1)&(h==1)).sum());fp=int(((y==0)&(h==1)).sum());fn=int(((y==1)&(h==0)).sum())
    precision=tp/(tp+fp) if tp+fp else 0.;recall=tp/(tp+fn) if tp+fn else 0.
    pr=[];rec=[]
    for t in sorted(set(score)):
        chosen=score>=t;pr.append(y[chosen].sum()/chosen.sum());rec.append(y[chosen].sum()/y.sum())
    pr.append(1.);rec.append(0.)
    return dict(roc_auc=auc_pairwise(y,score,np.ones(len(y))),average_precision=ap_threshold(y,score,np.ones(len(y))),genus_balanced_auc=auc_pairwise(y,score,w),genus_balanced_average_precision=ap_threshold(y,score,w),pr_auc_trapezoid=-float(np.trapezoid(pr,rec)),accuracy=float((y==h).mean()),precision=precision,recall=recall,f1=2*precision*recall/(precision+recall) if precision+recall else 0.)

# Normalize documented export-schema differences; never alter supplied raw tables.
# Older references use n_train_positive/test_positive and omit pooled classification.
normalization=[]
def normalized_table(name,path,pred):
    table=pd.read_csv(path/name).rename(columns={'n_train_positive':'n_train_pos','n_test_positive':'n_test_pos'})
    if name=='fold_metrics.csv' and 'genus_overlap' not in table:
        assert pred.groupby(['model','genus']).fold.nunique().eq(1).all()
        table['genus_overlap']=0;normalization.append(dict(source=('reference' if path==args.reference else 'clean'),table=name,derived_field='genus_overlap',definition='checked literal-genus fold disjointness'))
    if name=='summary.csv':
        if 'pooled_minus_mean_auc' not in table:
            table['pooled_minus_mean_auc']=table.pooled_roc_auc-table.roc_auc_mean
            normalization.append(dict(source=('reference' if path==args.reference else 'clean'),table=name,derived_field='pooled_minus_mean_auc',definition='pooled_roc_auc minus roc_auc_mean'))
        for key in ['accuracy','precision','recall','f1']:
            if 'pooled_'+key not in table:
                table['pooled_'+key]=[independent_metrics(pred[pred.model==model])[key] for model in table.model]
                normalization.append(dict(source=('reference' if path==args.reference else 'clean'),table=name,derived_field='pooled_'+key,definition='independently calculated from saved OOF labels/predicted classes'))
    return table

# Every saved numerical metric or count, including fold mean/SD and PR trapezoid AUC.
metric_rows=[]
for name,index in [('fold_metrics.csv',['seed','model','fold']),('summary.csv',['seed','model'])]:
    a=normalized_table(name,args.reference,ref);b=normalized_table(name,args.clean,clean)
    assert set(a.columns)==set(b.columns)
    joined=b.merge(a,on=index,suffixes=('_clean','_reference'),validate='one_to_one');assert len(joined)==len(a)==len(b)
    for column in [c for c in a.columns if c not in index and pd.api.types.is_numeric_dtype(a[c])]:
        for _,row in joined.iterrows():
            av=float(row[column+'_reference']);bv=float(row[column+'_clean'])
            metric_rows.append(dict(table=name,**{k:row[k] for k in index},metric=column,reference=av,clean=bv,delta=bv-av,same_rounded_4_decimals=round(av,4)==round(bv,4)))
metric_rows=pd.DataFrame(metric_rows);metric_rows.to_csv(D/'all_saved_metric_comparisons.csv',index=False)
independent=[]
for label,path,pred in [('reference',args.reference,ref),('clean',args.clean,clean)]:
    foldmetrics=normalized_table('fold_metrics.csv',path,pred);summary=normalized_table('summary.csv',path,pred)
    for model,p in pred.groupby('model'):
        for fold,q in [('pooled',p)]+list(p.groupby('fold')):
            saved=summary[summary.model==model].iloc[0] if fold=='pooled' else foldmetrics[(foldmetrics.model==model)&(foldmetrics.fold==fold)].iloc[0]
            for key,value in independent_metrics(q).items():
                expected=saved['pooled_'+key] if fold=='pooled' else saved[key]
                assert abs(value-expected)<1e-12
                independent.append(dict(environment=label,model=model,fold=fold,metric=key,independent_value=value,saved_value=expected,delta=value-expected))
independent=pd.DataFrame(independent);independent.to_csv(D/'independent_metric_reconstruction.csv',index=False)

# Compare audit-only family importance ranks, without writing a candidate inventory.
importance=[];old={};new={}
for model in sorted(ref.model.unique()):
    name='AUDIT_ONLY_cv_importance_'+model.replace(' ','_')+'.csv'
    a=pd.read_csv(args.reference/name).set_index('gene_family').importance.abs().sort_index()
    b=pd.read_csv(args.clean/name).set_index('gene_family').importance.abs().sort_index()
    assert a.index.equals(b.index);old[model]=a;new[model]=b
def top(s,k):return set(pd.DataFrame({'family':s.index,'importance':s.values}).sort_values(['importance','family'],ascending=[False,True]).head(k).family)
old['four_model_percentile_consensus']=pd.concat({m:s.rank(pct=True) for m,s in old.items()},axis=1).mean(axis=1)
new['four_model_percentile_consensus']=pd.concat({m:s.rank(pct=True) for m,s in new.items()},axis=1).mean(axis=1)
importance_details=[]
for model in old:
    a=old[model];b=new[model];ar=a.rank(ascending=False);br=b.rank(ascending=False)
    importance.append(dict(model=model,n_families=len(a),max_absolute_importance_difference=(a-b).abs().max(),changed_average_ranks=int((ar!=br).sum()),spearman=ar.corr(br),top50_overlap=len(top(a,50)&top(b,50)),top100_overlap=len(top(a,100)&top(b,100))))
    importance_details.extend(dict(model=model,gene_family=f,reference_importance=a[f],clean_importance=b[f],reference_rank=ar[f],clean_rank=br[f]) for f in a.index)
importance=pd.DataFrame(importance);importance.to_csv(D/'audit_only_importance_rank_comparison.csv',index=False)
pd.DataFrame(importance_details).to_csv(D/'audit_only_importance_family_detail.csv',index=False)
pd.DataFrame(normalization).drop_duplicates().to_csv(D/'export_schema_normalization.csv',index=False)
result=dict(strict_probability_comparison='FAIL' if merged.fails_original_strict_tolerance.any() else 'PASS',strict_rtol=1e-10,strict_atol=1e-12,strict_failure_preserved=True,
            prediction_rows=len(merged),strict_probability_mismatches=int(merged.fails_original_strict_tolerance.sum()),maximum_absolute_probability_difference=float(merged.absolute_probability_difference.max()),
            decision_label_changes=int(merged.decision_changed.sum()),exact_feature_masks_equal=True,exact_global_purity_descriptive_equal=True,
            changed_score_order_or_tie_pairs=int(orders.changed_pair_order_or_tie.sum()),saved_numerical_values_compared=len(metric_rows),maximum_saved_metric_difference=float(metric_rows.delta.abs().max()),all_saved_metrics_same_to_4_decimals=bool(metric_rows.same_rounded_4_decimals.all()),
            unique_changed_score_pairs=len({(x['model'],x['accession_a'],x['accession_b']) for x in changed_pairs}),changed_score_pairs_all_same_class=all(x['y_a']==x['y_b'] for x in changed_pairs),
            independent_metric_checks=len(independent),maximum_independent_metric_difference=float(independent.delta.abs().max()),audit_importance_ranks_unchanged=bool(importance.changed_average_ranks.eq(0).all()),
            no_retraining=True,no_tolerance_relaxation=True,scientific_parameters_changed=False,
            input_sha256={label+'/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for label,path in [('reference',args.reference),('clean',args.clean)] for p in path.glob('*.csv')})
(D/'comparison_summary.json').write_text(json.dumps(result,indent=2),encoding='utf8')
print(by_model.to_string(index=False));print(importance.to_string(index=False));print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))
if args.assert_strict:
    np.testing.assert_allclose(merged.y_prob_clean,merged.y_prob_reference,rtol=1e-10,atol=1e-12)
