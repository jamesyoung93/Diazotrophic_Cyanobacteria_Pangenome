"""Exclude ubiquitous controls after the verified cross-atlas correction.

Run with --workbook SOURCE.xlsx --matrix MATRIX.csv --output-dir OUTPUT_DIR.
The companion score_rules.py supplies the already audited first pass.
All source files, discovery inventory membership and numeric scores are preserved.
The additional rule changes tier eligibility, not score weights or thresholds.
"""
from pathlib import Path
import argparse
import json
import importlib.util
import numpy as np
import pandas as pd
# Load the exact packaged companion explicitly. Windows FileFinder can reject
# normal sibling imports near MAX_PATH even when Python can open the file.
_spec=importlib.util.spec_from_file_location('score_rules',Path(__file__).with_name('score_rules.py'))
_rules=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_rules)
DEFAULT,ABUNDANT,TIERS,flag,sha,propose,shortlist,summarize=(_rules.DEFAULT,_rules.ABUNDANT,_rules.TIERS,_rules.flag,_rules.sha,_rules.propose,_rules.shortlist,_rules.summarize)

HERE = Path(__file__).resolve().parent
RUN = HERE.parent/'inputs/raw'
CONTROL_LABEL = 'Ubiquitous control: no presence/absence contrast'


def tier_independent(row):
    # This new exclusion precedes both pre-existing control exclusions and score bands.
    if row['ubiquitous_control_flag']:
        return CONTROL_LABEL
    if row['nif_or_nitrogenase_audit']:
        return 'Core/diagnostic control'
    if row['housekeeping_lineage_flag']:
        return 'De-emphasize: housekeeping/lineage'
    score = row['scope_adjusted_story_score']
    if row['candidate_set'] == 'Model-Supported':
        if score >= 6:
            return TIERS[0]
        if score >= 4.5:
            return TIERS[1]
    if row['candidate_set'] == 'Highly Pure' and score >= 4.5:
        return TIERS[2]
    if score >= 4.5:
        return TIERS[3]
    if row['n_literature_active_phase_up_studies'] >= 2:
        return 'Tier D: active-phase-up proteomics background'
    return 'Tier E: lower current evidence'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workbook', type=Path, default=DEFAULT)
    ap.add_argument('--matrix', type=Path, default=RUN/'gene_family_matrix.csv')
    ap.add_argument('--prior-corrected-csv', type=Path, default=None)
    ap.add_argument('--independent-reconstruction', type=Path, default=HERE.parent/'results/reconstruction/score_family_reconstruction.json')
    ap.add_argument('--output-dir', type=Path, default=HERE.parent/'results/current167')
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    original_hash = sha(args.workbook)
    released = pd.read_excel(args.workbook, sheet_name='AllComposite')
    released_top = pd.read_excel(args.workbook, sheet_name='TopFamilies')
    prior = propose(released, True)
    if args.prior_corrected_csv is not None:
        prior_snapshot = pd.read_csv(args.prior_corrected_csv).set_index('gene_family')
        aligned = prior.set_index('gene_family').loc[prior_snapshot.index]
        assert np.allclose(aligned.scope_adjusted_story_score, prior_snapshot.scope_adjusted_story_score)
        assert aligned.priority_tier.eq(prior_snapshot.priority_tier).all()

    matrix = pd.read_csv(args.matrix, usecols=lambda c: c.startswith('GF_'))
    assert len(matrix) == 426 and matrix.shape[1] == 2286
    assert matrix.isin([0, 1]).all().all()
    counts = matrix.sum(axis=0)
    revised = prior.copy()
    revised['prior169_priority_tier'] = revised.priority_tier
    revised['verified_carrier_count'] = revised.gene_family.map(counts).astype(int)
    revised['noncarrier_genome_count'] = len(matrix)-revised.verified_carrier_count
    revised['ubiquitous_control_flag'] = revised.noncarrier_genome_count.eq(0)
    revised['directionality_contrast_estimable'] = revised.verified_carrier_count.gt(0) & revised.noncarrier_genome_count.gt(0)
    revised['prioritization_exclusion_reason'] = np.where(revised.ubiquitous_control_flag,
        'All 426 genomes carry this family; the noncarrier diazotrophy rate and presence/absence contrast are undefined.', '')
    controls = revised[revised.ubiquitous_control_flag].copy()
    ids = sorted(controls.gene_family)
    assert ids == ['GF_01297', 'GF_01899', 'GF_01945']
    assert controls.candidate_set.eq('Model-Supported').all()
    assert revised.verified_carrier_count.eq(revised.n_member_genomes).all()
    revised['priority_tier'] = revised.apply(tier_independent, axis=1)
    revised['score_revision_rule'] = ('Remove raw-ID atlas bonus and imported HGT penalty; retain HP offset 1.0; '
        'exclude ubiquitous families before tier assignment without changing their numeric scores or frozen inventory membership')
    revised_top = shortlist(revised)
    assert len(revised_top) == 167
    assert set(shortlist(prior).gene_family)-set(revised_top.gene_family) == {'GF_01297', 'GF_01945'}
    assert set(revised_top.gene_family)-set(shortlist(prior).gene_family) == set()
    assert revised.groupby('candidate_set').size().to_dict() == {'Highly Pure':981, 'Model-Supported':476}
    changed = revised[['gene_family','product','candidate_set','released_priority_tier','prior169_priority_tier','priority_tier',
        'released_scope_adjusted_story_score','scope_adjusted_story_score','ubiquitous_control_flag',
        'directionality_contrast_estimable','verified_carrier_count','noncarrier_genome_count','prioritization_exclusion_reason']].copy()
    changed['score_delta'] = changed.scope_adjusted_story_score-changed.released_scope_adjusted_story_score
    changed['tier_changed'] = changed.priority_tier.ne(changed.released_priority_tier)
    changed['tier_changed_vs169'] = changed.priority_tier.ne(changed.prior169_priority_tier)
    changed['released_shortlist'] = changed.released_priority_tier.isin(TIERS)
    changed['corrected_shortlist'] = changed.priority_tier.isin(TIERS)
    changed['shortlist_changed'] = changed.released_shortlist.ne(changed.corrected_shortlist)
    migrations = changed[changed.tier_changed | changed.score_delta.ne(0)].copy()
    controls = revised[revised.ubiquitous_control_flag].copy()
    ms = revised[revised.candidate_set.eq('Model-Supported')]
    ms473 = ms[ms.directionality_contrast_estimable]
    summary = summarize(revised)
    summary['model_importance_vs_score_spearman_scope'] = 'All 476 frozen Model-Supported families, including the three excluded ubiquitous controls; numeric scores unchanged.'
    summary['model_importance_vs_score_spearman_estimable473'] = float(ms473.consensus_rank_pct_mean.rank().corr(ms473.scope_adjusted_story_score.rank()))
    summary['model_importance_vs_score_spearman_estimable_n'] = len(ms473)
    summary['frozen_MS_count'] = len(ms)
    summary['ubiquitous_control_count'] = len(controls)
    summary['activeup_control_count'] = int(controls.n_literature_active_phase_up_studies.gt(0).sum())
    summary['Figure3_TierA_unchanged'] = set(revised_top[revised_top.priority_tier.eq(TIERS[0])].gene_family) == set(prior[prior.priority_tier.eq(TIERS[0])].gene_family)
    summary['Figure5_scope'] = 'Keep the frozen 476 Model-Supported families in diagnostic strata; 16/22 Tier A composition is unchanged. These full-inventory diagnostics include three ubiquitous controls excluded from shortlist eligibility.'

    independent = pd.read_json(args.independent_reconstruction).set_index('gene_family')
    check = revised.set_index('gene_family').loc[independent.index]
    # Raw reconstruction used source study measurements and membership mappings, not workbook scores.
    expected = independent.adjusted_score_if_omit_atlas_component_only + independent.hgt_passenger_flag.astype(int)*.75
    assert np.allclose(check.scope_adjusted_story_score, expected, atol=1e-9)
    expected_tiers = independent.tier_if_omit_atlas_component_and_imported_hgt_flag.copy()
    expected_tiers.loc[ids] = CONTROL_LABEL
    assert check.priority_tier.eq(expected_tiers).all()
    importance_records = []
    for name in ['Random_Forest','Gradient_Boosting','Logistic_Regression','XGBoost']:
        fp = args.matrix.parent/f'feature_importance_{name}.csv'
        imp = pd.read_csv(fp).set_index('gene_family').loc[ids,'importance']
        assert imp.eq(0).all()
        importance_records.extend({'gene_family':gf, 'model':name, 'importance':float(value)} for gf,value in imp.items())
    result = {'status':'Proposed revision for author review — ubiquitous-control exclusion',
        'source_workbook_name':args.workbook.name,'source_workbook_sha256':original_hash,
        'original_workbook_unchanged':sha(args.workbook)==original_hash,
        'matrix_sha256':sha(args.matrix),'prior169_csv_sha256':sha(args.prior_corrected_csv) if args.prior_corrected_csv else None,
        'independent_reconstruction_sha256':sha(args.independent_reconstruction),
        'rule':revised.score_revision_rule.iloc[0], 'composition_regex':ABUNDANT,
        'released':summarize(released),'prior169':summarize(prior),'corrected':summary,
        'score_changed_count':int(changed.score_delta.ne(0).sum()),
        'tier_changed_count':int(changed.tier_changed.sum()),
        'tier_changed_vs169_count':int(changed.tier_changed_vs169.sum()),
        'numeric_score_changed_vs169_count':0,
        'removed_ids':changed.loc[changed.released_shortlist & ~changed.corrected_shortlist,'gene_family'].tolist(),
        'removed_vs169_ids':['GF_01297','GF_01945'],'added_ids':[],
        'control_records':json.loads(controls[['gene_family','product','candidate_set','prior169_priority_tier','priority_tier',
            'scope_adjusted_story_score','n_literature_active_phase_up_studies','n_literature_mapped_studies',
            'consensus_rank_pct_mean','verified_carrier_count','noncarrier_genome_count']].to_json(orient='records')),
        'independent_raw_input_reconstruction_agreement':'All 1457 numeric scores match the existing independent raw-input reconstruction; tiers match after the independently matrix-verified ubiquitous-control exclusion.',
        'importance_checks':importance_records}
    datasets = {'released_AllComposite':released, 'released_TopFamilies':released_top,
        'prior169_AllComposite':prior,'prior169_TopFamilies':shortlist(prior),
        'corrected_AllComposite':revised,'corrected_TopFamilies':revised_top,
        'score_migrations':migrations,'tier_migrations':changed[changed.tier_changed],
        'ubiquitous_controls':controls,'control_migrations_vs169':changed[changed.tier_changed_vs169]}
    for name,df in datasets.items():
        df.to_csv(out/(name+'.csv'), index=False, encoding='utf-8-sig')
    payload = {k:{'columns':list(v.columns),'rows':json.loads(v.to_json(orient='values',double_precision=15))} for k,v in datasets.items()}
    payload['summary'] = result
    (out/'workbook_payload.json').write_text(json.dumps(payload,allow_nan=False),encoding='utf-8')
    (out/'score_revision_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    pd.DataFrame(importance_records).to_csv(out/'ubiquitous_control_importance_check.csv',index=False)
    (out/'README.txt').write_text('Proposed revision for author review. The original workbook, 476/981 discovery inventories and previous 169-family correction snapshot are unchanged.\n'
        'This pass excludes all three matrix-verified ubiquitous families before tier assignment. Numeric scores are retained for traceability.\n'
        'Corrected shortlist:167 =22 A +140 B +5 HP;65 have an active-up annotation. The prior169 pass and raw release remain available alongside corrected outputs.\n'
        'Use score_current.py with its companion score_rules.py; pass --workbook, --matrix, --prior-corrected-csv, --independent-reconstruction and --output-dir as needed.\n'
        'Full-inventory diagnostics use frozen476 MS; score decomposition uses162 shortlisted MS A+B.\n',encoding='utf-8')
    print(json.dumps({'corrected':summary,'controls':result['control_records'],
        'tier_changed_vs169_count':result['tier_changed_vs169_count'],
        'independent_check':result['independent_raw_input_reconstruction_agreement']},indent=2))

if __name__ == '__main__':
    main()
