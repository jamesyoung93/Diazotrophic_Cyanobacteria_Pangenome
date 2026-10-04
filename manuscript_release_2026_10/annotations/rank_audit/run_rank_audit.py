"""Regenerate rank comparisons from frozen saved model outputs, without fitting."""
from pathlib import Path
import argparse, json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,default=Path(__file__).resolve().parent/'inputs')
parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results')
parser.add_argument('--verify',action='store_true',help='Check against included reference tables after recomputing.')
args=parser.parse_args();I=args.inputs;D=args.out;D.mkdir(parents=True,exist_ok=True)
manifest=json.loads((Path(__file__).resolve().parent/'input_manifest.json').read_text())
for name,data in manifest.items():assert hashlib.sha256((I/name).read_bytes()).hexdigest()==data['staged_sha256']
inventory=pd.read_csv(I/'frozen_model_supported476.tsv',sep='\t').set_index('gene_family');assert len(inventory)==476
original={};corrected={};replay={};runtime=pd.read_csv(I/'current_global_replay_fold_importances.csv')
models=['Random Forest','Gradient Boosting','Logistic Regression','XGBoost']
for model in models:
    original[model]=pd.read_csv(I/f"released_{model.replace(' ','_')}.csv").set_index('gene_family').importance.abs()
    corrected[model]=pd.read_csv(I/f"corrected_{model.replace(' ','_')}.csv").set_index('gene_family').importance.abs()
    replay[model]=runtime[runtime.model==model].groupby('gene_family').importance.mean().abs()
original_consensus=pd.concat({m:s.rank(pct=True,ascending=True) for m,s in original.items()},axis=1).mean(axis=1)
corrected_consensus=pd.concat({m:s.rank(pct=True,ascending=True) for m,s in corrected.items()},axis=1).mean(axis=1)
np.testing.assert_allclose(original_consensus.loc[inventory.index],inventory.consensus_rank_pct_mean,rtol=1e-12,atol=1e-12)
comparison=[];details=[]
def top(series,k):
    return set(pd.DataFrame({'gene_family':series.index,'value':series.values}).sort_values(['value','gene_family'],ascending=[False,True]).head(k).gene_family)
def compare(old,new,source,model,scope):
    ids=sorted(set(old.index)&set(new.index))
    if scope=='positive_inventory':ids=sorted(set(ids)&set(inventory.index))
    a=old.loc[ids];b=new.loc[ids]
    row=dict(comparison=source,model=model,scope=scope,n_shared=len(ids),spearman=float(spearmanr(a,b).statistic))
    for k in [50,100]:
        aa=top(a,k);bb=top(b,k);row[f'top{k}_overlap']=len(aa&bb);row[f'top{k}_jaccard']=len(aa&bb)/len(aa|bb)
    comparison.append(row)
    ar=a.rank(ascending=False);br=b.rank(ascending=False)
    for family in ids:details.append(dict(comparison=source,model=model,scope=scope,gene_family=family,old_importance=old[family],new_importance=new[family],old_rank=float(ar[family]),new_rank=float(br[family])))
for model in models:
    for scope in ['all_shared_residual','positive_inventory']:
        compare(original[model],corrected[model],'released_vs_corrected',model,scope)
        compare(replay[model],corrected[model],'current_global_replay_vs_corrected',model,scope)
for scope in ['all_shared_residual','positive_inventory']:compare(original_consensus,corrected_consensus,'released_vs_corrected','four_model_consensus',scope)
pd.DataFrame(comparison).to_csv(D/'rank_comparison_summary.csv',index=False)
pd.DataFrame(details).to_csv(D/'rank_comparison_family_detail.csv',index=False)
coef=pd.read_csv(I/'exact_corrected_family_LR_coefficients_50fits.csv');coef['abs_coef']=coef.coefficient.abs()
frozen_ids=sorted(inventory.index)
index=pd.MultiIndex.from_product([range(42,52),range(1,6),frozen_ids],names=['seed','fold','feature'])
full=coef.set_index(['seed','fold','feature']).reindex(index)
full['selected']=full.coefficient.notna();full['coefficient_filled']=full.coefficient.fillna(0);full['abs_coef_filled']=full.abs_coef.fillna(0)
means=full.groupby(level=['seed','feature']).agg(mean_abs_coef=('abs_coef_filled','mean'),mean_signed_coef=('coefficient_filled','mean'),selected_folds=('selected','sum')).reset_index()
means['rank_within_frozen476']=means.groupby('seed').mean_abs_coef.rank(ascending=False,method='average')
seq=[]
consensus_ranks=inventory.consensus_rank_pct_mean.rank(ascending=False)
lr_ranks=original['Logistic Regression'].loc[frozen_ids].rank(ascending=False)
for family in inventory.index:
    folds=full.xs(family,level='feature');chosen=folds[folds.selected];ranks=means[means.feature==family]
    seq.append(dict(gene_family=family,product=inventory.loc[family,'product'],released_consensus_rank_within476=float(consensus_ranks[family]),released_LR_rank_within476=float(lr_ranks[family]),selected_of50=int(folds.selected.sum()),positive_coefficient_fits=int((chosen.coefficient>0).sum()),negative_coefficient_fits=int((chosen.coefficient<0).sum()),zero_coefficient_fits=int((chosen.coefficient==0).sum()),positive_fraction_when_selected=float((chosen.coefficient>0).mean()),mean_signed_coefficient_including_unselected_zero=float(folds.coefficient_filled.mean()),min_selected_coefficient=float(chosen.coefficient.min()),max_selected_coefficient=float(chosen.coefficient.max()),median_rank_within476=float(ranks.rank_within_frozen476.median()),min_rank_within476=float(ranks.rank_within_frozen476.min()),max_rank_within476=float(ranks.rank_within_frozen476.max()),top50_partitions=sum(family in top(g.set_index('feature').mean_abs_coef,50) for _,g in means.groupby('seed')),top100_partitions=sum(family in top(g.set_index('feature').mean_abs_coef,100) for _,g in means.groupby('seed'))))
pd.DataFrame(seq).to_csv(D/'all_frozen476_repeated_LR_stability.csv',index=False)
verification={'input_hashes_pass':True,'no_model_fitting':True,'outputs':3}
if args.verify:
    for name in ['rank_comparison_summary.csv','rank_comparison_family_detail.csv','all_frozen476_repeated_LR_stability.csv']:
        actual=pd.read_csv(D/name);expected=pd.read_csv(Path(__file__).resolve().parent/'reference'/name)
        pd.testing.assert_frame_equal(actual,expected,check_exact=False,rtol=1e-11,atol=1e-12)
    verification['reference_tables_match']=True
(D/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf8')
print(json.dumps(verification,indent=2))
