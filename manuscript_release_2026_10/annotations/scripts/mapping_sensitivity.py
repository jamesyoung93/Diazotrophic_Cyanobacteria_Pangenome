"""Read-only mapping sensitivity with a frozen 1,311-family rank reference.

Neither source scores nor the current corrected shortlist are modified. For each
row-exclusion policy, the primary conservative analysis disables the condensate
annotation when no retained row supplies the released family's best score/rank.
A secondary analysis places the best retained value against the original fixed
family reference, leaving all other families' ranks and the denominator intact.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd

E=Path(__file__).resolve().parent
ROOT=E.parent
TIERS=['Tier A: story-leading accessory','Tier B: strong cross-evidence accessory',
       'Tier HP: evidence-supported high-purity marker','Tier C: cross-evidence diagnostic/support']

def flag(s):
    return s.astype(str).str.lower().isin(['true','1','1.0'])

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def tier(r):
    if r.ubiquitous_control_flag: return 'Ubiquitous control: no presence/absence contrast'
    if r.nif_or_nitrogenase_audit: return 'Core/diagnostic control'
    if r.housekeeping_lineage_flag: return 'De-emphasize: housekeeping/lineage'
    units=round(r.scope_adjusted_story_score*20)
    if r.candidate_set=='Model-Supported' and units>=120: return TIERS[0]
    if r.candidate_set=='Model-Supported' and units>=90: return TIERS[1]
    if r.candidate_set=='Highly Pure' and units>=90: return TIERS[2]
    if units>=90:return TIERS[3]
    if r.n_literature_active_phase_up_studies>=2:return 'Tier D: active-phase-up proteomics background'
    return 'Tier E: lower current evidence'

def bonus(pct):
    return np.select([pct.le(.01),pct.le(.05),pct.le(.1),pct.le(.2)],[.5,.4,.3,.2],default=0.)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ranking',type=Path,default=ROOT/'inputs/raw/condensate_protein_ranking.csv')
    ap.add_argument('--map',type=Path,default=ROOT/'inputs/raw/genome_protein_family_map.tsv')
    ap.add_argument('--family-overlay',type=Path,default=ROOT/'inputs/raw/condensate_ranking_family_overlay.tsv')
    ap.add_argument('--corrected-csv',type=Path,default=ROOT/'results/current167/corrected_AllComposite.csv')
    ap.add_argument('--output-dir',type=Path,default=ROOT/'results/mapping_sensitivity')
    args=ap.parse_args(); out=args.output_dir;out.mkdir(exist_ok=True,parents=True)
    inputs=[args.ranking,args.map,args.family_overlay,args.corrected_csv]
    source_hashes={p.name:digest(p) for p in inputs}
    pg=pd.read_csv(args.map,sep='\t',dtype=str)
    rows=pd.read_csv(args.ranking)
    memberships=pg.groupby('protein_accession').gene_family.agg(set)
    ambiguous=set(memberships[memberships.map(len).gt(1)].index)
    # Reconstruct the released first-membership mapping in source row order.
    mapper={}
    for r in pg.itertuples(index=False):mapper.setdefault(r.protein_accession,r.gene_family)
    rows['direct_family']=rows.uniprot_id.map(mapper)
    rows['cluster_family']=rows.cluster_id.map(mapper)
    rows['gene_family']=rows.direct_family.fillna(rows.cluster_family)
    rows['cluster_only']=rows.direct_family.isna()&rows.cluster_family.notna()
    rows['selected_accession']=rows.uniprot_id.where(rows.direct_family.notna(),rows.cluster_id)
    rows['selected_ambiguous']=rows.selected_accession.isin(ambiguous)
    rows=rows[rows.gene_family.notna()].copy()
    ref=pd.read_csv(args.family_overlay,sep='\t').set_index('gene_family')
    rebuilt=rows.groupby('gene_family').agg(n=('driver_score','size'),best=('driver_score','max'),bestrank=('driver_rank_cyano','min')).sort_values(['best','bestrank'],ascending=[False,True])
    rebuilt['rank']=np.arange(1,len(rebuilt)+1)
    assert len(rebuilt)==len(ref)==1311
    assert np.allclose(rebuilt.best,ref.loc[rebuilt.index,'condensate_best_score'],rtol=0,atol=1e-14)
    assert rebuilt['rank'].eq(ref.loc[rebuilt.index,'condensate_family_rank']).all()
    assert rebuilt.n.eq(ref.loc[rebuilt.index,'condensate_mapped_rows']).all()
    d=pd.read_csv(args.corrected_csv).set_index('gene_family')
    for c in ['ubiquitous_control_flag','nif_or_nitrogenase_audit','housekeeping_lineage_flag']:d[c]=flag(d[c])
    assert d.apply(tier,axis=1).eq(d.priority_tier).all()
    assert np.allclose(d.condensate_score,bonus(d.condensate_family_percentile),atol=1e-9)
    base_short=set(d.index[d.priority_tier.isin(TIERS)])
    base_A=set(d.index[d.priority_tier.eq(TIERS[0])])
    assert len(base_short)==167 and len(base_A)==22
    summaries={}
    policies={'exclude_selected_ambiguous':~rows.selected_ambiguous,
              'exclude_cluster_only_fallback':~rows.cluster_only,
              'exclude_both':~(rows.selected_ambiguous|rows.cluster_only)}
    for policy,keep in policies.items():
        retained=rows[keep]
        remaining=retained.groupby('gene_family').agg(retained_rows=('driver_score','size'),retained_best_score=('driver_score','max'),retained_best_rank=('driver_rank_cyano','min'))
        f=ref[['condensate_mapped_rows','condensate_best_score','condensate_best_rank_cyano','condensate_family_rank','condensate_family_percentile']].join(remaining)
        f['lost_all_rows']=f.retained_rows.isna()
        f['released_best_supported']=f.retained_best_score.eq(f.condensate_best_score)&f.retained_best_rank.eq(f.condensate_best_rank_cyano)
        # Float values are round-tripped CSV decimals; tolerate only roundoff, not rank differences.
        f['released_best_supported']=(np.isclose(f.retained_best_score,f.condensate_best_score,rtol=0,atol=1e-14)&f.retained_best_rank.eq(f.condensate_best_rank_cyano))
        f['conservative_frozen_percentile']=f.condensate_family_percentile.where(f.released_best_supported)
        f['fixed_reference_percentile']=f.condensate_family_percentile
        f.loc[f.lost_all_rows,'fixed_reference_percentile']=np.nan
        for gf,r in f.loc[~f.lost_all_rows&~f.released_best_supported].iterrows():
            other=ref.drop(index=gf)
            better=other.condensate_best_score.gt(r.retained_best_score)
            equal=other.condensate_best_score.eq(r.retained_best_score)
            tie_before=other.condensate_best_rank_cyano.lt(r.retained_best_rank)|(other.condensate_best_rank_cyano.eq(r.retained_best_rank)&other.index.to_series().lt(gf))
            fixed_rank=1+int((better|(equal&tie_before)).sum())
            f.loc[gf,'fixed_reference_percentile']=fixed_rank/1311
        f['conservative_component']=bonus(f.conservative_frozen_percentile)
        f['fixed_reference_component']=bonus(f.fixed_reference_percentile)
        f['baseline_component']=bonus(f.condensate_family_percentile)
        assert (f.fixed_reference_component<=f.baseline_component+1e-9).all()
        assert (f.conservative_component<=f.fixed_reference_component+1e-9).all()
        f.to_csv(out/(policy+'_family_annotation.csv'))
        rows[~keep].to_csv(out/(policy+'_excluded_protein_rows.csv'),index=False)
        summary={'excluded_protein_rows':int((~keep).sum()),'retained_protein_rows':int(keep.sum()),
                 'families_with_removed_rows':int(rows.loc[~keep,'gene_family'].nunique()),
                 'families_losing_all_rows':int(f.lost_all_rows.sum()),
                 'families_losing_released_best_annotation':int((~f.released_best_supported).sum()),
                 'retained_mapped_families':len(remaining),'analyses':{}}
        for method in ['conservative','fixed_reference']:
            revised=d.copy()
            revised['prior_tier']=d.priority_tier
            revised['prior_condensate_score']=d.condensate_score
            revised['prior_adjusted_score']=d.scope_adjusted_story_score
            revised['condensate_score']=f[method+'_component'].reindex(d.index).fillna(0)
            delta_units=((revised.condensate_score-d.condensate_score)*20).round().astype(int)
            for c in ['composite_evidence_score','accessory_story_score','scope_adjusted_story_score']:
                revised[c]=((d[c]*20).round().astype(int)+delta_units)/20
            revised['priority_tier']=revised.apply(tier,axis=1)
            top=revised[revised.priority_tier.isin(TIERS)]
            A=set(revised.index[revised.priority_tier.eq(TIERS[0])])
            migrations=revised[revised.priority_tier.ne(revised.prior_tier)]
            detailcols=['product','candidate_set','prior_tier','priority_tier','prior_condensate_score','condensate_score','prior_adjusted_score','scope_adjusted_story_score','n_literature_active_phase_up_studies']
            migrations[detailcols].to_csv(out/(policy+'_'+method+'_tier_migrations.csv'))
            revised[detailcols].to_csv(out/(policy+'_'+method+'_all_scores.csv'))
            summary['analyses'][method]={'families_numeric_score_changes':int(delta_units.ne(0).sum()),
                'tier_counts':top.priority_tier.value_counts().to_dict(),'shortlist_count':len(top),
                'activeup_shortlist':int(top.n_literature_active_phase_up_studies.gt(0).sum()),
                'shortlist_removed_ids':sorted(base_short-set(top.index)),
                'shortlist_added_ids':sorted(set(top.index)-base_short),
                'TierA_count':len(A),'TierA_lost_ids':sorted(base_A-A),'TierA_added_ids':sorted(A-base_A),
                'tier_migrations':migrations.reset_index()[['gene_family']+detailcols].to_dict('records')}
        summaries[policy]=summary
    result={'source_hashes':source_hashes,'original_files_unchanged':all(digest(p)==source_hashes[p.name] for p in inputs),
        'baseline_reconstruction':'All 1311 family ranks, best scores and mapped-row counts agree with released overlay. All 1457 v2 tiers and condensate components reproduced independently.',
        'reference_denominator':1311,'mapped_protein_rows':len(rows),'unique_best_score_ties':int(ref.condensate_best_score.duplicated(keep=False).sum()),
        'methods':{
            'conservative':'Remove specified protein mappings. Keep the released family annotation and rank only if a retained row still supplies its released best score/rank; otherwise set the condensate annotation missing and its heuristic component to zero. Never rerank other families.',
            'fixed_reference':'If released best support is lost but another row remains, position its best retained score (then best cyano rank, then family ID) among the original 1311-family reference, excluding its own previous entry and retaining denominator1311. This is a sensitivity placement, not a newly estimated rank distribution. All unaffected ranks remain fixed.',
            'ambiguous':'Exclude only rows whose selected lookup accession has multiple GF memberships. An unused ambiguous cluster accession does not invalidate a direct unambiguous mapping.',
            'cluster_only':'Exclude only rows lacking a direct accession map and assigned through cluster_id. Direct mappings, including the22 direct-vs-cluster conflicts, retain the source direct-first policy.'},
        'sensitivities':summaries}
    (out/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
