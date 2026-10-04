"""Verify saved identity/grouping sensitivity outputs without fitting models."""
from pathlib import Path
import ast,json,os
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parents[1];D=Path(os.environ['MODEL_AUDIT_OUTPUT'])/os.environ['MODEL_AUDIT_SCOPE']
helpers=[n for n in ast.parse((OUT/'analysis/independent_helpers.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name in ['pairwise_auc','threshold_ap']]
ns={'np':np};exec(compile(ast.Module(body=helpers,type_ignores=[]),'independent_scores','exec'),ns)
X=pd.read_csv(OUT/'inputs/gene_family_matrix.csv',index_col=0)
meta=pd.read_csv(OUT/'inputs/complete_genomes_with_proteins.csv').set_index('assembly_accession')
common=sorted(set(X.index)&set(meta.index));X=X.loc[common];meta=meta.loc[common];y=meta.is_diazotroph.astype(int)
pred=pd.read_csv(D/'oof_predictions.csv');folds=pd.read_csv(D/'fold_metrics.csv');summary=pd.read_csv(D/'summary.csv')
assign=pd.read_csv(D/'fold_assignments.csv');features=pd.read_csv(D/'fold_family_features.csv');params=pd.read_csv(D/'lr_scalers_and_coefficients.csv')
old=pd.read_csv(OUT/'inputs/model_fold_assignments.csv').set_index('assembly_accession')
counts=pd.read_csv(D/'fold_cohort_and_mask_counts.csv');rows=[];reconstructed=0
for scenario,part in assign.groupby('scenario'):
    assert not part.assembly_accession.duplicated().any()
    lookup=part.set_index('assembly_accession');ids=sorted(lookup.index)
    if scenario!='connected_group_repartition':
        np.testing.assert_array_equal(lookup.fold,old.loc[lookup.index,'fold'])
    else:
        assert lookup.audit_group.nunique()==71
        assert lookup.loc['GCF_047302775.1','fold']==lookup.loc['GCF_047302855.1','fold']
        for a,b in [('Leptolyngbya','[Leptolyngbya]'),('Phormidium','[Phormidium]'),('Synechococcus','[Synechococcus]')]:
            assert part[part.stored_genus.isin([a,b])].fold.nunique()==1
    for fold,bucket in part.groupby('fold'):
        test=bucket.assembly_accession.tolist();train=[a for a in ids if a not in set(test)]
        assert not set(lookup.loc[test,'audit_group'])&set(lookup.loc[train,'audit_group'])
        pp=X.loc[train].eq(1);n=pp.sum();k=pp.mul(y.loc[train],axis=0).sum()
        possible=sorted(n.index[(n>0)&(k/n>0.1)&(k/n<0.9)])
        expected=[c for c in possible if lookup.loc[pp.index[pp[c]],'audit_group'].nunique()>=3]
        saved=features[(features.scenario==scenario)&(features.fold==fold)].gene_family.tolist()
        assert expected==saved
        for model in summary.model.unique():
            p=pred[(pred.scenario==scenario)&(pred.fold==fold)&(pred.model==model)]
            assert p.assembly_accession.tolist()==test
            np.testing.assert_array_equal(p.y_true,y.loc[test]);np.testing.assert_array_equal(p.genus,meta.loc[test,'genus'])
            np.testing.assert_array_equal(p.audit_group,lookup.loc[test,'audit_group'])
            f=folds[(folds.scenario==scenario)&(folds.fold==fold)&(folds.model==model)].iloc[0]
            if model=='Logistic Regression' or model.endswith('_LR'):
                if model=='assembly_length_only_LR':z=meta.loc[ids,['total_ungapped_length']].astype(float)
                elif model=='retained_family_count_only_LR':z=X.loc[ids].eq(1).sum(axis=1).to_frame('retained_family_count')
                else:z=X.loc[ids,saved]
                coeff=params[(params.scenario==scenario)&(params.fold==fold)&(params.model==model)]
                assert coeff.feature.tolist()==z.columns.tolist()
                means=z.loc[train].mean().to_numpy();scales=z.loc[train].std(ddof=0).to_numpy();scales[z.loc[train].nunique().to_numpy()==1]=1
                np.testing.assert_allclose(means,coeff.train_mean,rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(scales,coeff.train_scale,rtol=1e-12,atol=1e-12)
                logits=((z.loc[test].to_numpy()-coeff.train_mean.to_numpy())/coeff.train_scale.to_numpy())@coeff.coefficient.to_numpy()+coeff.intercept.iloc[0]
                probability=1/(1+np.exp(-logits))
                np.testing.assert_allclose(probability,p.y_prob,rtol=1e-11,atol=1e-12);reconstructed+=len(p)
            if p.y_true.nunique()==2:
                for weight_label,w in [('genome',np.ones(len(p))),('genus',1/p.groupby('genus').genus.transform('size')),('audit_group',1/p.groupby('audit_group').audit_group.transform('size'))]:
                    for suffix,helper in [('auc','pairwise_auc'),('average_precision','threshold_ap')]:
                        metric=('roc_auc' if suffix=='auc' else suffix) if weight_label=='genome' else weight_label+'_balanced_'+suffix
                        value=ns[helper](p.y_true,p.y_prob,w);assert abs(value-f[metric])<1e-12
                        rows.append(dict(scenario=scenario,model=model,fold=fold,metric=metric,delta=value-f[metric]))
for (scenario,model),p in pred.groupby(['scenario','model']):
    population=assign[assign.scenario==scenario]
    assert len(p)==p.assembly_accession.nunique()==len(population)
    s=summary[(summary.scenario==scenario)&(summary.model==model)].iloc[0]
    f=folds[(folds.scenario==scenario)&(folds.model==model)]
    assert int(s.auc_defined_folds)==int(f.roc_auc.notna().sum())
    for weight_label,w in [('genome',np.ones(len(p))),('genus',1/p.groupby('genus').genus.transform('size')),('audit_group',1/p.groupby('audit_group').audit_group.transform('size'))]:
        for suffix,helper in [('auc','pairwise_auc'),('average_precision','threshold_ap')]:
            metric=('roc_auc' if suffix=='auc' else suffix) if weight_label=='genome' else weight_label+'_balanced_'+suffix
            value=ns[helper](p.y_true,p.y_prob,w);assert abs(value-s['pooled_'+metric])<1e-12
            assert abs(f[metric].mean()-s[metric+'_mean'])<1e-12
            assert abs(f[metric].std(ddof=1)-s[metric+'_sd'])<1e-12
            rows.append(dict(scenario=scenario,model=model,fold='pooled',metric=metric,delta=value-s['pooled_'+metric]))
pd.DataFrame(rows).to_csv(D/'independent_metric_verification.csv',index=False)
result=dict(status='PASS',independent_score_checks=len(rows),maximum_score_discrepancy=float(pd.DataFrame(rows).delta.abs().max()),
            original_fold_ids_retained_for_all_local_exclusion_arms=True,group_and_sample_overlap_checks_passed=True,independently_verified_masks=assign.scenario.nunique()*5,
            LR_predictions_reconstructed=reconstructed,total_prediction_rows=len(pred),single_class_test_folds=int((folds.groupby(['scenario','fold']).first().n_test_pos==0).sum()))
(D/'verification_summary.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result,indent=2))
