"""Proposed revision for author review: remove unverifiable cross-atlas ID inputs.

Usage: python score_proposed_revision.py [--workbook SOURCE.xlsx] [--output-dir DIR]
Requires pandas, numpy and openpyxl for reading; never modifies the source workbook.
All released columns/rows remain in released CSVs. Changed columns are also retained
as released_* fields in revised CSVs. Fixed thresholds and other weights are retained.
The shortlist is regenerated from tiers, with gene_family as a final display-only tie
breaker. The original sorter left ties unspecified; this does not alter membership.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT = HERE.parent / 'inputs/released/proteomics_composite_family_evidence.xlsx'
ABUNDANT = r'clp|chaperon|dnak|groe|grol|dnaj|grpe|release factor|translation|elongation factor|typa|trna|nyn|ribosom|isopropylmalate|threonine synthase|ketol-acid|dihydroxy-acid|aminomutase|amino acid'
TIERS = ['Tier A: story-leading accessory', 'Tier B: strong cross-evidence accessory',
         'Tier HP: evidence-supported high-purity marker', 'Tier C: cross-evidence diagnostic/support']

def flag(s):
    return s.astype(str).str.lower().isin(['true', '1', '1.0'])

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def tier(row):
    if row['nif_or_nitrogenase_audit']:
        return 'Core/diagnostic control'
    if row['housekeeping_lineage_flag']:
        return 'De-emphasize: housekeeping/lineage'
    units = round(20 * row['scope_adjusted_story_score'])
    if row['candidate_set'] == 'Model-Supported' and units >= 120:
        return TIERS[0]
    if row['candidate_set'] == 'Model-Supported' and units >= 90:
        return TIERS[1]
    if row['candidate_set'] == 'Highly Pure' and units >= 90:
        return TIERS[2]
    if units >= 90:
        return TIERS[3]
    if row['n_literature_active_phase_up_studies'] >= 2:
        return 'Tier D: active-phase-up proteomics background'
    return 'Tier E: lower current evidence'

def alignment(row):
    pieces = [row['candidate_set']]
    up = int(row['n_literature_active_phase_up_studies'])
    if up:
        pieces.append(f'{up} active-phase-up proteomics studies')
    elif row['n_literature_nfix_response_studies'] > 0:
        pieces.append(f"{int(row['n_literature_nfix_response_studies'])} broad N-fix-responsive studies (not scored)")
    if row['condensate_top20']:
        pieces.append('condensate-ranked')
    if row['primary_bridge_broad_ge1u_ge1f']:
        pieces.append('genus-proxy breadth')
    return '; '.join(pieces)

def shortlist(df):
    out = df[df.priority_tier.isin(TIERS)].copy()
    out['_tier_order'] = out.priority_tier.map(dict(zip(TIERS, range(4))))
    return out.sort_values(['_tier_order','scope_adjusted_story_score','accessory_story_score',
                            'composite_evidence_score','gene_family'], ascending=[True,False,False,False,True]).drop(columns='_tier_order')

def crosstab(df):
    return pd.crosstab(df.priority_tier, df.candidate_set).to_dict()

def summarize(df):
    top = shortlist(df)
    ms = df[df.candidate_set.eq('Model-Supported')]
    ab = top[top.candidate_set.eq('Model-Supported')]
    a = top[top.priority_tier.eq(TIERS[0])]
    b = top[top.priority_tier.eq(TIERS[1])]
    components = ['model_score','literature_proteomics_score','related_atlas_score',
                  'condensate_score','morphotype_breadth_score']
    avg = float(ab.scope_adjusted_story_score.mean())
    matched = a['product'].str.lower().fillna('').str.contains(ABUNDANT)
    rho = ms.consensus_rank_pct_mean.rank().corr(ms.scope_adjusted_story_score.rank())
    return {
        'all_families':len(df),'all_crosstab':crosstab(df),'shortlist_count':len(top),
        'shortlist_crosstab':crosstab(top),'activeup_shortlist':int(top.n_literature_active_phase_up_studies.gt(0).sum()),
        'activeup_shortlist_crosstab':crosstab(top[top.n_literature_active_phase_up_studies.gt(0)]),
        'MS_AplusB_count':len(ab),'MS_AplusB_mean_score':avg,
        'MS_AplusB_score_decomposition':{c:{'mean':float(ab[c].mean()),'fraction_of_mean_score':float(ab[c].mean()/avg)} for c in components},
        'TierB_zero_scored_external_components':int((b[['literature_proteomics_score','related_atlas_score','condensate_score']].sum(axis=1).eq(0)).sum()),
        'TierA_both_morphotype_flags':int((flag(a.primary_bridge_broad_ge1u_ge1f)&flag(a.strict_unicellular_breadth_ge3u_ge1f)).sum()),
        'TierA_regex_composition_count':int(matched.sum()),'TierA_regex_composition_denominator':len(a),
        'TierA_regex_composition_fraction':float(matched.mean()),
        'TierA_regex_composition_ids':a.loc[matched,'gene_family'].tolist(),
        'TierA_story_role_counts':a.story_role.value_counts().to_dict(),
        'TierA_module_bin_counts':a.module_bin.value_counts().to_dict(),
        'model_importance_vs_score_spearman':float(rho),
        'median_model_importance_percentile_TierA':float(a.consensus_rank_pct_mean.median()),
        'median_model_importance_percentile_other_than_A_B':float(ms.loc[~ms.priority_tier.isin(TIERS[:2]),'consensus_rank_pct_mean'].median()),
        'shortlist_top10_condensate':int(flag(top.condensate_top10).sum()),
        'shortlist_broad_flag':int(flag(top.primary_bridge_broad_ge1u_ge1f).sum()),
        'shortlist_strict_flag':int(flag(top.strict_unicellular_breadth_ge3u_ge1f).sum()),
    }

def propose(released, drop_hgt=True):
    d = released.copy()
    changed = ['priority_tier','composite_evidence_score','accessory_story_score','scope_adjusted_story_score',
               'related_atlas_score','related_atlas_scope_adjustment_score','related_atlas_applicability',
               'story_penalty','related_atlas_match','exact_product_match','evidence_alignment_summary',
               'aryal_product','aryal_consensus_rank_pct_mean','aryal_gene','aryal_module_bin',
               'frac_assemblies_proximal','fisher_p_vs_background']
    if drop_hgt:
        changed.append('hgt_passenger_flag')
    for c in changed:
        if c in d:
            d['released_'+c] = d[c]
    for c in ['nif_or_nitrogenase_audit','housekeeping_lineage_flag','hgt_passenger_flag','condensate_top20','primary_bridge_broad_ge1u_ge1f']:
        d[c] = flag(d[c])
    d['related_atlas_score'] = 0.0
    d['related_atlas_scope_adjustment_score'] = 0.0
    d['hp_inventory_offset'] = d.candidate_set.eq('Highly Pure').astype(float)
    d['related_atlas_applicability'] = 'disabled: no validated cross-atlas family mapping'
    d['related_atlas_match'] = False
    d['exact_product_match'] = False
    # Imported metadata must not masquerade as verified evidence in corrected columns.
    for c in ['aryal_product','aryal_consensus_rank_pct_mean','aryal_gene','aryal_module_bin',
              'frac_assemblies_proximal','fisher_p_vs_background']:
        d[c] = np.nan
    if drop_hgt:
        d['hgt_passenger_flag'] = False
    components = ['model_score','literature_proteomics_score','condensate_score','morphotype_breadth_score']
    units = (d[components] * 20).round().sum(axis=1).astype(int)
    penalties = (d.nif_or_nitrogenase_audit.astype(int)*20 + d.housekeeping_lineage_flag.astype(int)*20
                 + d.hgt_passenger_flag.astype(int)*15)
    d['composite_evidence_score'] = units/20
    d['story_penalty'] = penalties/20
    d['accessory_story_score'] = (units-penalties)/20
    d['scope_adjusted_story_score'] = (units-penalties+d.candidate_set.eq('Highly Pure').astype(int)*20)/20
    d['priority_tier'] = d.apply(tier,axis=1)
    d['evidence_alignment_summary'] = d.apply(alignment,axis=1)
    d['score_revision_status'] = 'Proposed revision for author review'
    d['score_revision_rule'] = 'Remove raw-ID atlas bonus and imported HGT penalty; retain HP inventory offset 1.0' if drop_hgt else 'Atlas-only sensitivity: remove raw-ID atlas bonus; retain imported HGT penalty and HP offset 1.0'
    return d.sort_values(['scope_adjusted_story_score','accessory_story_score','composite_evidence_score',
                         'n_literature_active_phase_up_studies','gene_family'],ascending=[False,False,False,False,True])
