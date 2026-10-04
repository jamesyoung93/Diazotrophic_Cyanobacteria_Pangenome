"""Subset independently computed ranking diagnostics to the current score shortlist.

This does not refit models or alter any score, tier, or released ranking.
"""
from pathlib import Path
import json, argparse
import pandas as pd

P=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--all-composite-csv',type=Path,default=P/'results/current167/corrected_AllComposite.csv')
ap.add_argument('--rank-dir',type=Path,default=P/'results/model_rank')
ap.add_argument('--output-dir',type=Path,default=P/'results/rank_sensitivity')
args=ap.parse_args();OUT=args.output_dir;OUT.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(args.all_composite_csv)
lr=pd.read_csv(args.rank_dir/'all_frozen476_repeated_LR_stability.csv')
rank=pd.read_csv(args.rank_dir/'rank_comparison_family_detail.csv')

cons = rank.loc[rank.comparison.eq('released_vs_corrected') & rank.model.eq('four_model_consensus') & rank.scope.eq('positive_inventory')]
merged = df[['gene_family','product','priority_tier','candidate_set','model_score','consensus_rank_pct_mean','scope_adjusted_story_score']].merge(lr.drop(columns='product'),on='gene_family',validate='one_to_one').merge(cons[['gene_family','old_rank','new_rank']],on='gene_family',validate='one_to_one')
merged['coefficient_has_both_signs'] = merged.positive_coefficient_fits.gt(0) & merged.negative_coefficient_fits.gt(0)
merged['rank_shift'] = merged.new_rank - merged.old_rank
summaries = {}
for name, mask in [('current_TierA22',merged.priority_tier.str.startswith('Tier A:')),('current_TierB140',merged.priority_tier.str.startswith('Tier B:')),('current_MS_shortlist162',merged.priority_tier.str.startswith(('Tier A:','Tier B:'))),('frozen476',pd.Series(True,index=merged.index))]:
    x = merged.loc[mask].copy()
    summaries[name] = {
        'n':len(x), 'selected_all50':int(x.selected_of50.eq(50).sum()),
        'both_coefficient_signs':int(x.coefficient_has_both_signs.sum()),
        'positive_only_when_selected':int((x.positive_coefficient_fits.gt(0)&x.negative_coefficient_fits.eq(0)).sum()),
        'negative_only_when_selected':int((x.negative_coefficient_fits.gt(0)&x.positive_coefficient_fits.eq(0)).sum()),
        'negative_all50_ids':x.loc[x.negative_coefficient_fits.eq(50),'gene_family'].tolist(),
        'positive_all50_ids':x.loc[x.positive_coefficient_fits.eq(50),'gene_family'].tolist(),
        'old_consensus_top50_n':int(x.old_rank.le(50).sum()),
        'new_consensus_top50_n':int(x.new_rank.le(50).sum()),
        'old_consensus_median_rank':float(x.old_rank.median()),
        'new_consensus_median_rank':float(x.new_rank.median()),
        'max_absolute_rank_shift':float(x.rank_shift.abs().max()),
    }
    x.to_csv(OUT / (name+'_rank_stability.csv'),index=False)
result = {'status':'Read-only sensitivity analysis; score unchanged','summaries':summaries,
          'model_score_rule':'2.0 for frozen Model-Supported inventory membership; 1.0 for Highly Pure membership. No importance or rank enters the model_score component.',
          'displayed_released_rank_rule':'Per model: absolute raw importance ranked ascending with pandas average tied ranks and pct=True (highest importance gives percentile 1); consensus is the mean available within-model percentile. The released per-model universe is 627 families, not the 476-family positive inventory.',
          'interpretation':'Ranks and signs are sensitivity evidence. Logistic-regression coefficients are conditional on other selected, correlated binary family features and standardized by the model pipeline; mixed signs do not alone disprove a marginal association or show a protein lacks biological function. Heuristic score bands are not independent validation.',
          'recommendation':'Retain the explicit historical-inventory model weight and label Figure3 importance as released/historical. Supply the corrected model rank comparison and current-shortlist coefficient stability table. Changing the flat inventory component into corrected rank weighting would define a new, unvalidated prioritization method and should be an explicit author decision, with prespecified weights/thresholds and sensitivity analysis rather than an unannounced bug correction.'}
(OUT/'rank_stability_summary.json').write_text(json.dumps(result,indent=2),encoding='utf8')
(OUT/'rank_stability_interpretation.txt').write_text('\n'.join([result['status'],result['model_score_rule'],result['displayed_released_rank_rule'],result['interpretation'],result['recommendation'],'','Exact subset results:',json.dumps(summaries,indent=2)]),encoding='utf8')
print(json.dumps(result,indent=2))
