"""Independently reconstruct released scores from pre-composite source inputs.

No production scoring function is imported and no released score, tier, active-up
flag, morphotype flag, or condensate top-percentile flag is used as an input.
The workbook is opened only after reconstruction, for external comparison.

Usage: python scripts/reconstruct_annotations.py [--inputs inputs] [--output-dir results/reconstruction]
Requires pandas, numpy, openpyxl. All inputs are read-only.
"""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

AP = argparse.ArgumentParser(description=__doc__)
PACKAGE = Path(__file__).resolve().parents[1]
AP.add_argument('--inputs', type=Path, default=PACKAGE/'inputs')
AP.add_argument('--output-dir', type=Path, default=PACKAGE/'results/reconstruction')
args = AP.parse_args()
ROOT = PACKAGE
INPUTS = args.inputs.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)

def boolval(x):
    return str(x).strip().lower() in {'true', '1', 'yes', 'y'}

def metric_dict(text):
    # Parse the already-normalized source-study measurement text, not outcome flags.
    result = {}
    for token in str(text).split(';'):
        if '=' in token:
            key, value = token.strip().split('=', 1)
            result[key.strip()] = value.strip()
    return result

def val(metrics, key):
    try:
        return float(metrics.get(key, 'nan'))
    except ValueError:
        return float('nan')

def active_up(row):
    m = metric_dict(row['evidence_metric'])
    study = row['paper_short']
    if study == 'Panda et al. 2025':
        return (min(val(m, 'Two-way ANOVA p value nitrate'),
                    val(m, 'Two-way ANOVA p value Interaction')) <= .05
                or val(m, 'Two-way ANOVA p value Interaction') <= .05) and val(m, 'D-_D+(logFC)') >= .25
    if study == 'Sandh et al. 2014':
        ratio = val(m, 'Log2 Ratio (Heterocyst/Filaments)')
        if math.isnan(ratio):
            ratio = val(m, '24h Log2 Ratio (Het/Fil) (Present study)')
        return ratio >= 1
    if study == 'Welkie et al. 2014':
        ordered = [(k, val(m, k)) for k in ['CT_D0', 'CT_D3', 'CT_L0', 'CT_L3']]
        ordered = [(k, v) for k, v in ordered if math.isfinite(v)]
        if not ordered:
            return False
        high = max(ordered, key=lambda kv: kv[1])
        low = min(v for k, v in ordered)
        return high[1] - low >= 1 and high[0] in {'CT_D0', 'CT_D3'}
    if study == 'Held et al. 2022':
        abundance, span = val(m, 'mean_relative_abundance'), val(m, 'dynamic_range')
        match = re.search(r'@([0-9.]+)h', m.get('max', ''))
        hour = float(match[1]) if match else float('nan')
        return abundance != 0 and span / abundance >= .5 and 6 <= hour <= 15
    return False

def normalize(s):
    if pd.isna(s):
        return ''
    return ' '.join(str(s).strip().lower().split()).replace('"', '')

paths = {
 'inventory':INPUTS/'raw/morphotype_bridge_family_inventory.tsv',
 'normalized_proteomics':INPUTS/'raw/normalized_proteins.csv',
 'related_atlas_cache':INPUTS/'raw/aryal_task7_hgt_proximity_from_github.csv',
 'raw_condensate_ranking':INPUTS/'raw/condensate_protein_ranking.csv',
 'protein_family_map':INPUTS/'raw/genome_protein_family_map.tsv',
 'family_matrix':INPUTS/'raw/gene_family_matrix.csv',
 'genome_metadata':INPUTS/'raw/genome_metadata.csv',
 'workbook':INPUTS/'released/proteomics_composite_family_evidence.xlsx',
 'released_scorer':PACKAGE/'provenance/released_composite_scorer.py',
 'atlas_join_generator':PACKAGE/'provenance/released_atlas_join_generator.py',
}

hashes = {key: {'path': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
          for key, path in paths.items()}

inventory = pd.read_csv(paths['inventory'], sep='\t')
inventory = inventory[inventory.candidate_set.isin(['Model-Supported', 'Highly Pure'])]
assert inventory.gene_family.is_unique
inventory = inventory.set_index('gene_family')

# Recount morphology from genome-level rows, excluding workbook/inventory flags.
matrix = pd.read_csv(paths['family_matrix'], index_col=0)
metadata = pd.read_csv(paths['genome_metadata'], usecols=['assembly_accession', 'genus', 'is_diazotroph'])
metadata = metadata.drop_duplicates('assembly_accession').set_index('assembly_accession')
metadata = metadata.loc[metadata.index.intersection(matrix.index)]
diazo = metadata[metadata.is_diazotroph.map(boolval)]
strict_genera = {'Candidatus', 'Crocosphaera', 'Euhalothece', 'Gloeothece', 'Halothece',
                'Rippkaea', 'Synechococcus', 'Synechocystis', 'cyanobacterium'}
broad_genera = strict_genera | {'Chroococcidiopsis', 'Cyanobacterium', 'Pleurocapsa'}
carrier_counts = {}
for tag, genera in [('strict', strict_genera), ('broad', broad_genera)]:
    for side, mask in [('unicellular', diazo.genus.isin(genera)), ('filamentous_proxy', ~diazo.genus.isin(genera))]:
        carrier_counts[f'{tag}_{side}_genomes'] = matrix.loc[diazo.index[mask]].astype(bool).sum()

# Reclassify every raw normalized protein row and collapse UNIQUE studies/family.
proteins = pd.read_csv(paths['normalized_proteomics'])
contexts = {'Panda et al. 2025': 'unicellular', 'Welkie et al. 2014': 'unicellular',
            'Sandh et al. 2014': 'heterocyst', 'Held et al. 2022': 'nonheterocystous filament'}
mapped_studies, up_studies, up_contexts = defaultdict(set), defaultdict(set), defaultdict(set)
up_proteins, study_rows = defaultdict(set), Counter()
class_records = []
welkie_ties = []
for row in proteins.to_dict('records'):
    is_up = bool(active_up(row))
    families = str(row['mapped_gene_families']).split(';') if row['map_status'] == 'mapped_to_family' else []
    families = [family for family in families if family and family != 'nan']
    for family in families:
        mapped_studies[family].add(row['paper_short'])
        study_rows[(family, row['paper_short'])] += 1
        if is_up:
            up_studies[family].add(row['paper_short'])
            up_contexts[family].add(contexts[row['paper_short']])
            up_proteins[family].add(row['source_protein_key'])
    class_records.append({'source_record_id': row['source_record_id'], 'study': row['paper_short'],
                          'family_ids': families, 'active_up_reconstructed': is_up})
    if row['paper_short'] == 'Welkie et al. 2014':
        m = metric_dict(row['evidence_metric'])
        values = {k: val(m, k) for k in ['CT_D0','CT_D3','CT_L0','CT_L3'] if math.isfinite(val(m, k))}
        if values:
            maxima = [k for k,v in values.items() if v == max(values.values())]
            if len(maxima) > 1 and max(values.values())-min(values.values()) >= 1:
                welkie_ties.append({'record': row['source_record_id'], 'maxima': maxima, 'families': families})

# Rebuild condensate mapping and ranks from protein rows. The first accession
# mapping is the released bridge choice; ambiguity is reported separately.
pg = pd.read_csv(paths['protein_family_map'], sep='\t', dtype=str)
acc_family_sets = pg.groupby('protein_accession').gene_family.agg(lambda v: sorted(set(v)))
ambiguous_accessions = acc_family_sets[acc_family_sets.map(len) > 1]
protein_map = {}
for acc, family in zip(pg.protein_accession, pg.gene_family):
    protein_map.setdefault(acc, family)
cond_rows = pd.read_csv(paths['raw_condensate_ranking'])
cond_by_family = defaultdict(list)
mapping_conflicts = []
for row in cond_rows.to_dict('records'):
    direct, cluster = protein_map.get(row['uniprot_id']), protein_map.get(row['cluster_id'])
    if direct and cluster and direct != cluster:
        mapping_conflicts.append([row['uniprot_id'], row['cluster_id'], direct, cluster])
    family = direct or cluster
    if family:
        cond_by_family[family].append(row)
cond_records = []
for family, rows in cond_by_family.items():
    cond_records.append({'gene_family': family, 'nrows': len(rows),
                         'best_score': max(float(r['driver_score']) for r in rows),
                         'best_rank': min(float(r['driver_rank_cyano']) for r in rows)})
cond_records.sort(key=lambda r: (-r['best_score'], r['best_rank'], r['gene_family']))
cond_lookup = {}
for rank, row in enumerate(cond_records, 1):
    row.update(rank=rank, percentile=rank / len(cond_records))
    cond_lookup[row['gene_family']] = row

# Reconstruct the actual related-atlas source: a shared IDENTIFIER join. This is
# deliberately separate from the normalized-product crosswalk described in text.
atlas = pd.read_csv(paths['related_atlas_cache'])
assert atlas.gene_family.is_unique
atlas_by_id = atlas.set_index('gene_family').to_dict('index')
atlas_by_product = defaultdict(list)
for row in atlas.to_dict('records'):
    atlas_by_product[normalize(row['product'])].append(row['gene_family'])

housekeeping_patterns = ['lineage marker', 'housekeeping', 'ribosomal protein', '30s ribosomal',
                         '50s ribosomal', 'ribosome-binding factor', 'elongation factor',
                         'translation initiation factor', 'translation elongation',
                         'aminoacyl-trna synthetase', 'trna synthetase']

def assign_tier(candidate, units, core, housekeeping, up_count):
    if core:
        return 'Core/diagnostic control'
    if housekeeping:
        return 'De-emphasize: housekeeping/lineage'
    if candidate == 'Model-Supported' and units >= 120:
        return 'Tier A: story-leading accessory'
    if candidate == 'Model-Supported' and units >= 90:
        return 'Tier B: strong cross-evidence accessory'
    if candidate == 'Highly Pure' and units >= 90:
        return 'Tier HP: evidence-supported high-purity marker'
    if units >= 90:
        return 'Tier C: cross-evidence diagnostic/support'
    if up_count >= 2:
        return 'Tier D: active-phase-up proteomics background'
    return 'Tier E: lower current evidence'

reconstructed = []
for family, source in inventory.iterrows():
    product, category = str(source['product']), str(source['module_bin'])
    core = category == 'Nitrogenase / nif machinery' or re.search(r'\bnif|nitrogenase|molybdenum-iron', product, re.I) is not None
    housekeeping = any(word in f'{product} {category}'.lower() for word in housekeeping_patterns)
    broad_flag = carrier_counts['broad_unicellular_genomes'][family] >= 1 and carrier_counts['broad_filamentous_proxy_genomes'][family] >= 1
    strict_flag = carrier_counts['strict_unicellular_genomes'][family] >= 3 and carrier_counts['strict_filamentous_proxy_genomes'][family] >= 1
    u, n, contexts_count = len(up_studies[family]), len(mapped_studies[family]), len(up_contexts[family])
    match = atlas_by_id.get(family)
    exact = bool(match) and str(match['product']).strip().lower() == product.strip().lower() and bool(product.strip())
    atlas_pct = float(match['consensus_rank_pct_mean']) if match else float('nan')
    hgt = bool(match) and boolval(match['hgt_passenger_flag'])
    cond = cond_lookup.get(family)
    percentile = cond['percentile'] if cond else float('nan')
    # Integer twentieths remove floating-point/rounding ambiguity.
    component_units = {
        'model_score': 40 if source['candidate_set']=='Model-Supported' else 20,
        'literature_proteomics_score': 20*min(u,3) + (10 if contexts_count==2 else 20 if contexts_count>=3 else 0),
        'related_atlas_score': 10*bool(match) + 5*exact + 5*(atlas_pct <= .25),
        'condensate_score': next((units for cutoff, units in [(.01,10),(.05,8),(.1,6),(.2,4)] if percentile <= cutoff),0),
        'morphotype_breadth_score': 25*int(broad_flag)+30*int(strict_flag),
    }
    total = sum(component_units.values())
    penalty = 20*int(core)+20*int(housekeeping)+15*int(hgt)
    adjustment = max(20-component_units['related_atlas_score'],0) if source['candidate_set']=='Highly Pure' else 0
    adjusted = total-penalty+adjustment
    row = {'gene_family': family, 'candidate_set': source['candidate_set'], 'product': product,
           **{k: v/20 for k,v in component_units.items()}, 'composite_evidence_score': total/20,
           'story_penalty': penalty/20, 'accessory_story_score': (total-penalty)/20,
           'related_atlas_scope_adjustment_score': adjustment/20, 'scope_adjusted_story_score': adjusted/20,
           'priority_tier': assign_tier(source['candidate_set'], adjusted, core, housekeeping, u),
           'n_literature_active_phase_up_studies': u, 'n_literature_mapped_studies': n,
           'active_phase_up_context_count': contexts_count, 'related_atlas_match': bool(match),
           'exact_product_match': exact, 'hgt_passenger_flag': hgt,
           'nif_or_nitrogenase_audit': core, 'housekeeping_lineage_flag': housekeeping,
           'primary_bridge_broad_ge1u_ge1f': bool(broad_flag), 'strict_unicellular_breadth_ge3u_ge1f': bool(strict_flag),
           'condensate_family_percentile': percentile, 'condensate_mapped_rows': cond['nrows'] if cond else None,
           'atlas_product_by_same_id': match['product'] if match else None,
           'any_normalized_product_match': normalize(product) in atlas_by_product,
           'normalized_product_matching_atlas_ids': atlas_by_product.get(normalize(product), []),
           'atlas_consensus_percentile': atlas_pct,
           'adjusted_score_if_omit_atlas_component_only':
               (adjusted - (component_units['related_atlas_score'] if source['candidate_set']=='Model-Supported' else 0))/20,
           'tier_if_omit_atlas_component_only': assign_tier(source['candidate_set'],
               adjusted - (component_units['related_atlas_score'] if source['candidate_set']=='Model-Supported' else 0),
               core, housekeeping, u),
           'tier_if_omit_atlas_component_and_imported_hgt_flag': assign_tier(source['candidate_set'],
               adjusted - (component_units['related_atlas_score'] if source['candidate_set']=='Model-Supported' else 0) + 15*int(hgt),
               core, housekeeping, u),
    }
    reconstructed.append(row)
rebuilt = pd.DataFrame(reconstructed).set_index('gene_family')

# Validation targets are loaded only here, after every component is calculated.
actual = pd.read_excel(paths['workbook'], sheet_name='AllComposite').set_index('gene_family')
top = pd.read_excel(paths['workbook'], sheet_name='TopFamilies')
shared_cols = [c for c in rebuilt if c in actual]
discrepancies = {}
for col in shared_cols:
    left, right = rebuilt[col], actual.loc[rebuilt.index,col]
    if col == 'active_phase_up_context_count':
        # Release leaves this joined count blank when no active-up row exists.
        # An independently recounted empty set has size0; normalize only here.
        assert left.loc[right.isna()].eq(0).all()
        right = right.fillna(0)
    if pd.api.types.is_numeric_dtype(left):
        equal = np.isclose(pd.to_numeric(left),pd.to_numeric(right),equal_nan=True,atol=1e-10)
    else:
        equal = left.fillna('').astype(str).eq(right.fillna('').astype(str)).to_numpy()
    if not equal.all():
        discrepancies[col] = [{'gene_family': f, 'reconstructed': str(left.loc[f]), 'released': str(right.loc[f])}
                              for f in left.index[~equal]]

assert not discrepancies, f'Independent reconstruction differs: {list(discrepancies)}'
chosen = rebuilt[rebuilt.priority_tier.str.startswith(('Tier A:','Tier B:','Tier HP:','Tier C:'))]
same_id = rebuilt[rebuilt.related_atlas_match]
wrong_product = same_id[~same_id.exact_product_match]
false_annotation = same_id[~same_id.any_normalized_product_match]
missed_annotation = rebuilt[(~rebuilt.related_atlas_match)&rebuilt.any_normalized_product_match]
score_ties = chosen.groupby(['priority_tier','scope_adjusted_story_score','accessory_story_score','composite_evidence_score']).size()
score_ties = score_ties[score_ties>1]
change = rebuilt[rebuilt.priority_tier != rebuilt.tier_if_omit_atlas_component_and_imported_hgt_flag]
atlas_only_change = rebuilt[rebuilt.priority_tier != rebuilt.tier_if_omit_atlas_component_only]
atlas_only_top = rebuilt[rebuilt.tier_if_omit_atlas_component_only.str.startswith(('Tier A:','Tier B:','Tier HP:','Tier C:'))]
top_order_keys = [
    (0 if r.priority_tier.startswith('Tier A:') else 1 if r.priority_tier.startswith('Tier B:') else 2,
     -r.scope_adjusted_story_score, -r.accessory_story_score, -r.composite_evidence_score)
    for r in top.itertuples(index=False)
]
quantile_bounds = {str(cut):{'last_in_rank':math.floor(cut*len(cond_records)),
                           'denominator_all_mapped_families':len(cond_records)} for cut in [.01,.05,.1,.2]}

report = {
 'method':'Independent input reconstruction; no production score functions imported; workbook only validation target.',
 'sources':hashes,
 'candidate_count':len(rebuilt), 'reconstructed_tier_counts':rebuilt.priority_tier.value_counts().to_dict(),
 'unique_selected_count':len(chosen), 'top_sheet_count':len(top),
 'top_sheet_sorted_by_documented_keys':top_order_keys==sorted(top_order_keys),
 'selected_id_symmetric_difference':sorted(set(chosen.index)^set(top.gene_family)),
 'column_discrepancies':discrepancies,
 'comparison_missingness_note':'For active_phase_up_context_count only, released blanks correspond to independently verified empty context sets and are compared as0. All other fields retain their original missingness semantics.',
 'context_count_blank_to_verified_zero':int(actual.active_phase_up_context_count.isna().sum()),
 'raw_source_counts':{'normalized_protein_rows':len(proteins),'matrix_rows':len(matrix),'diazotroph_metadata':len(diazo),
                      'raw_condensate_rows':len(cond_rows),'mapped_condensate_rows':sum(len(x) for x in cond_by_family.values()),
                      'condensate_families':len(cond_records),'atlas_cache_rows':len(atlas)},
 'proteomics':{'positive_candidate_families':int(rebuilt.n_literature_active_phase_up_studies.gt(0).sum()),
               'positive_selected_families':int(chosen.n_literature_active_phase_up_studies.gt(0).sum()),
               'max_mapped_studies':int(rebuilt.n_literature_mapped_studies.max()),
               'max_active_up_studies':int(rebuilt.n_literature_active_phase_up_studies.max()),
               'count_violations_up_exceeds_mapped':int((rebuilt.n_literature_active_phase_up_studies>rebuilt.n_literature_mapped_studies).sum()),
               'welkie_dynamic_max_ties':welkie_ties},
 'atlas_join':{'id_matches':len(same_id),'exact_product_same_id_matches':int(same_id.exact_product_match.sum()),
               'same_id_different_product':len(wrong_product),'id_matches_without_any_normalized_product_match':len(false_annotation),
               'product_matches_without_id_match':len(missed_annotation),
               'tierA_id_matches':int(chosen[chosen.priority_tier.str.startswith('Tier A:')].related_atlas_match.sum()),
               'tierB_id_matches':int(chosen[chosen.priority_tier.str.startswith('Tier B:')].related_atlas_match.sum()),
               'tierA_id_matches_without_any_product_match':list(false_annotation[false_annotation.priority_tier.str.startswith('Tier A:')].index),
               'hgt_flags_imported_on_id':int(rebuilt.hgt_passenger_flag.sum()),
               'hgt_flags_with_nonmatching_same_id_product':int(wrong_product.hgt_passenger_flag.sum()),
               'nonmatching_products_examples':wrong_product[['product','atlas_product_by_same_id','priority_tier','related_atlas_score','hgt_passenger_flag']].head(20).reset_index().to_dict('records')},
 'atlas_only_omission_sensitivity':{
     'description':'Set related-atlas additive component A=0 and recompute the existing Highly Pure max(1-A,0) scope adjustment. Preserve all other released rules, including the imported HGT penalty. Diagnostic only; not substituted for released tiers.',
     'changed_tier_count':len(atlas_only_change),
     'tier_counts':rebuilt.tier_if_omit_atlas_component_only.value_counts().to_dict(),
     'shortlist_count':len(atlas_only_top),
     'removed_from_shortlist':sorted(set(chosen.index)-set(atlas_only_top.index)),
     'added_to_shortlist':sorted(set(atlas_only_top.index)-set(chosen.index)),
     'changes':atlas_only_change[['product','priority_tier','scope_adjusted_story_score','related_atlas_score','adjusted_score_if_omit_atlas_component_only','tier_if_omit_atlas_component_only']].reset_index().to_dict('records')},
 'omission_sensitivity':{'description':'Diagnostic only: omit identifier-based atlas score for MS and imported HGT penalty; HP scope total remains1.0. Not an alternative recommended score.',
                          'changed_tier_count':len(change), 'changes':change[['priority_tier','tier_if_omit_atlas_component_and_imported_hgt_flag']].reset_index().to_dict('records')},
 'ties':{'selected_score_tuple_tie_groups':len(score_ties),'selected_families_in_ties':int(score_ties.sum()),
          'groups':[{'tier':key[0],'adjusted':key[1],'unadjusted':key[2],'sum':key[3],'n':int(n)} for key,n in score_ties.items()],
          'condensate_best_score_ties':int(pd.DataFrame(cond_records).duplicated('best_score',keep=False).sum()),
          'condensate_both_sortkey_ties':int(pd.DataFrame(cond_records).duplicated(['best_score','best_rank'],keep=False).sum())},
 'condensate':{'quantile_boundaries':quantile_bounds,'mapping_conflicts':mapping_conflicts,
               'ambiguous_accession_count':len(ambiguous_accessions)},
 'thresholds':{'exactly_6_nonexcluded_MS':rebuilt[(rebuilt.candidate_set=='Model-Supported')&(rebuilt.scope_adjusted_story_score==6)&(~rebuilt.nif_or_nitrogenase_audit)&(~rebuilt.housekeeping_lineage_flag)].index.tolist(),
               'exactly_4_5_nonexcluded':rebuilt[(rebuilt.scope_adjusted_story_score==4.5)&(~rebuilt.nif_or_nitrogenase_audit)&(~rebuilt.housekeeping_lineage_flag)].index.tolist(),
               'excluded_with_adjusted_ge4_5':rebuilt[(rebuilt.scope_adjusted_story_score>=4.5)&(rebuilt.nif_or_nitrogenase_audit|rebuilt.housekeeping_lineage_flag)][['priority_tier','scope_adjusted_story_score']].reset_index().to_dict('records')},
}
rebuilt.reset_index().to_json(OUT/'score_family_reconstruction.json',orient='records',indent=2)
(OUT/'score_audit_results.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps({'tier_counts':report['reconstructed_tier_counts'],
                  'shortlist_symmetric_difference':report['selected_id_symmetric_difference'],
                  'discrepancy_counts':{k:len(v) for k,v in discrepancies.items()},
                  'atlas_join_summary':{k:v for k,v in report['atlas_join'].items() if k!='nonmatching_products_examples'},
                  'atlas_only_sensitivity':report['atlas_only_omission_sensitivity']},indent=2))
