"""Independent mask/metric/scaler checks plus declared official-ID partition test."""
from pathlib import Path
E=Path(__file__).resolve().parents[1]
source=(E/'analysis/verify_sensitivity.py').read_text()
source=source.replace("D=OUT/'model_sample_identity_sensitivity'","D=OUT/'model_official_taxonomy_sensitivity'")
a=source.index("    if scenario!='connected_group_repartition':")
b=source.index("    for fold,bucket in part.groupby('fold'):",a)
source=source[:a]+'''    taxonomy=pd.read_csv(OUT/'inputs/taxonomy_genus_audit.csv').set_index('assembly_accession')
    resolved=taxonomy.current_genus_tax_id.notna()
    keep=resolved if scenario=='official_genus_resolved422' else resolved&(taxonomy.public_assembly_status=='current')
    assert set(ids)==set(taxonomy.index[keep])
    numeric_ids=sorted(taxonomy.loc[ids,'current_genus_tax_id'].astype(int).unique())
    shuffled=np.array(numeric_ids,dtype=int);np.random.default_rng(42).shuffle(shuffled)
    expected_group_fold={int(g):1+i%5 for i,g in enumerate(shuffled)}
    for accession in ids:
        group=int(taxonomy.loc[accession,'current_genus_tax_id'])
        assert lookup.loc[accession,'fold']==expected_group_fold[group]
        assert int(lookup.loc[accession,'audit_group'])==group
    assert lookup.audit_group.nunique()==len(numeric_ids)
'''+source[b:]
source=source.replace('original_fold_ids_retained_for_all_local_exclusion_arms=True','official_genus_numeric_sort_seed42_splits_verified=True')
source=source.replace('independently_verified_masks=20','independently_verified_masks=10')
exec(compile(source,'independent_official_taxonomy_verifier','exec'),{'__file__':str(E/'analysis/verify_official_taxonomy.py')})
