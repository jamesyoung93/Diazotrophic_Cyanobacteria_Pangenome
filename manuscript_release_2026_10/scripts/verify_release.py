"""Check the frozen release and manuscript-facing summaries without dependencies."""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import math
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_manifest(path, base, nested=True):
    data = json.loads(path.read_text(encoding='utf-8'))
    hashes = data['files'] if nested else data
    for name, expected in hashes.items():
        target = base / name
        assert target.is_file(), f'Missing: {target}'
        assert digest(target) == expected, f'Hash mismatch: {target}'
    return len(hashes)


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def check_figure3():
    """Check the display against frozen candidate rows and model-input labels."""
    rows = read_csv(ROOT / 'annotations/results/figures/Figure3_source.csv')
    tier = {r['gene_family']: r for r in read_csv(
        ROOT / 'annotations/expected/corrected_TopFamilies.csv')
        if r['priority_tier'].split(':')[0] == 'Tier A'}
    assert len(rows) == len(tier) == 22
    assert {r['gene_family'] for r in rows} == set(tier)

    # These are display assignments, separate from all score-bearing bins.
    groups = [
        ('Amino-acid biosynthesis',
         {'GF_00051', 'GF_00143', 'GF_00741', 'GF_01376', 'GF_01387'}),
        ('Protein quality control',
         {'GF_00290', 'GF_00614', 'GF_00658', 'GF_01689', 'GF_01799', 'GF_02262'}),
        ('Translation/RNA-associated products',
         {'GF_00623', 'GF_00771', 'GF_01479', 'GF_01542'}),
        ('Tetrapyrrole-associated products', {'GF_00300', 'GF_01233'}),
        ('Other products',
         {'GF_00122', 'GF_00303', 'GF_00536', 'GF_00916', 'GF_01449'}),
    ]
    group_for = {family: (index, label)
                 for index, (label, families) in enumerate(groups, 1)
                 for family in families}
    assert set(group_for) == set(tier)

    matrix = read_csv(ROOT / 'models/inputs/gene_family_matrix.csv')
    labels = read_csv(ROOT / 'models/inputs/public_status/cohort_status_reconciled.csv')
    assert len(matrix) == len(labels) == 426
    assert len(matrix[0]) - 1 == 2286
    assert len({r[''] for r in matrix}) == 426
    assert len({r['assembly_accession'] for r in labels}) == 426
    assert {r[''] for r in matrix} == {r['assembly_accession'] for r in labels}
    assert {r['is_diazotroph'] for r in labels} == {'True', 'False'}
    positive = {r['assembly_accession'] for r in labels if r['is_diazotroph'] == 'True'}
    assert len(positive) == 112
    prevalence = {}
    for family in tier:
        assert all(r[family] in ('0', '1') for r in matrix), family
        carriers = {r[''] for r in matrix if r[family] == '1'}
        prevalence[family] = (len(carriers), len(carriers & positive))

    organism_codes = [
        ('Crocosphaera subtropica ATCC 51142', 'Cr'),
        ('Trichodesmium erythraeum IMS101', 'Tr'),
        ('Nostoc punctiforme PCC 73102', 'No'),
        ('Cyanothece sp. PCC 7822', 'Cy'),
    ]
    gene_mapping = read_csv(ROOT / 'figures/Figure4_sources/gene_to_protein_top100_agreement.csv')
    mapped_families = {r['mapped_gene_family'] for r in gene_mapping}
    overlap = set(tier) & mapped_families
    assert overlap == {'GF_01387'}

    def equal_number(left, right, field, family):
        if left == '' or right == '':
            assert left == right, (family, field, left, right)
        else:
            assert math.isclose(float(left), float(right), rel_tol=0, abs_tol=1e-12), (
                family, field, left, right)

    for row_index, row in enumerate(rows, 1):
        family = row['gene_family']
        source = tier[family]
        carriers, nif_positive = prevalence[family]
        group_index, group_label = group_for[family]
        assert int(row['display_row']) == row_index
        assert int(row['display_group_index']) == group_index
        assert row['display_group'] == group_label
        assert int(row['panel_genomes']) == 426
        assert int(row['panel_nifHDK_positive']) == 112
        assert int(row['carrier_count']) == carriers == int(source['verified_carrier_count'])
        assert int(row['nifHDK_positive_carrier_count']) == nif_positive
        purity = nif_positive / carriers
        equal_number(row['carrier_purity_fraction'], purity, 'carrier_purity_fraction', family)
        equal_number(source['diazotroph_pct_mean'], purity, 'diazotroph_pct_mean', family)
        assert row['carrier_purity_percent_display'] == f'{100 * purity:.1f}'
        for field in ('product', 'priority_tier', 'candidate_set', 'active_phase_up_organisms',
                      'primary_bridge_broad_ge1u_ge1f', 'strict_unicellular_breadth_ge3u_ge1f'):
            assert row[field] == source[field], (family, field)
        for field in ('scope_adjusted_story_score', 'consensus_rank_pct_mean',
                      'n_literature_active_phase_up_studies', 'n_literature_mapped_studies',
                      'condensate_family_percentile'):
            equal_number(row[field], source[field], field, family)
        assert row['corrected_score_display'] == f"{float(source['scope_adjusted_story_score']):.2f}"
        for raw, display in (('consensus_rank_pct_mean', 'model_percentile_display'),
                             ('condensate_family_percentile', 'condensate_percentile_display')):
            expected = 100 * float(source[raw]) if source[raw] else ''
            equal_number(row[display], expected, display, family)
        organisms = {x.strip() for x in source['active_phase_up_organisms'].split(';') if x.strip()}
        assert organisms <= {name for name, code in organism_codes}
        assert row['active_up_organism_codes'] == '+'.join(
            code for name, code in organism_codes if name in organisms)
        assert int(row['supplied_gene_top100_mapping_overlap']) == int(source['gene_top100_overlap'])
        assert bool(int(row['supplied_gene_top100_mapping_overlap'])) == (family in overlap)
        assert row['marker_display'] == ('*' if family in overlap else '')
        if family in overlap:
            hits = [r for r in gene_mapping if r['mapped_gene_family'] == family]
            assert len(hits) == 1
            hit = hits[0]
            assert row['marker_provenance'] == (
                'figures/Figure4_sources/gene_to_protein_top100_agreement.csv; '
                f"supplied feature {hit['collab_feature']}; rank {hit['collab_rank']}; mapped to {family}")
        else:
            assert row['marker_provenance'] == ''

    ordered = sorted(tier, key=lambda family: (
        group_for[family][0], -prevalence[family][1] / prevalence[family][0], family))
    assert [r['gene_family'] for r in rows] == ordered
    for extension in ('png', 'pdf'):
        assert digest(ROOT / f'figures/Figure3.{extension}') == digest(
            ROOT / f'annotations/results/figures/Figure3_TierA_evidence.{extension}')
    return {'figure3_rows': len(rows), 'figure3_display_groups': len(groups),
            'figure3_matrix_carrier_checks': len(prevalence)}


def main():
    counts = {
        'release_files': check_manifest(ROOT / 'RELEASE_SHA256.json', ROOT),
        'frozen_scientific_files': check_manifest(
            ROOT / 'validation/SCIENTIFIC_REFERENCE_SHA256.json', ROOT),
        'annotation_package_files': check_manifest(
            ROOT / 'annotations/PACKAGE_MANIFEST.json', ROOT / 'annotations'),
        'model_package_files': check_manifest(
            ROOT / 'models/bundle_sha256.json', ROOT / 'models', nested=False),
    }
    all_rows = read_csv(ROOT / 'annotations/results/current167/corrected_AllComposite.csv')
    shortlist = read_csv(ROOT / 'annotations/results/current167/corrected_TopFamilies.csv')
    assert len(all_rows) == len({r['gene_family'] for r in all_rows}) == 1457
    assert Counter(r['candidate_set'] for r in all_rows) == {'Model-Supported': 476, 'Highly Pure': 981}
    assert len(shortlist) == len({r['gene_family'] for r in shortlist}) == 167
    tiers = Counter(r['priority_tier'].split(':')[0] for r in shortlist)
    assert tiers == {'Tier A': 22, 'Tier B': 140, 'Tier HP': 5}, tiers
    assert sum(float(r['n_literature_active_phase_up_studies']) > 0 for r in shortlist) == 65
    assert all(r['ubiquitous_control_flag'].lower() not in ('true', '1') for r in shortlist)

    # Recount the historical displayed census from the actual score-bearing field.
    observed = Counter((r['module_bin'], r['candidate_set']) for r in all_rows)
    archived = read_csv(ROOT / 'docs/legacy_table2_keyword_counts.csv')
    expected = {(r['Keyword bin'], group): int(r[group + ' families'])
                for r in archived for group in ('Model-Supported', 'Highly Pure')}
    assert dict(observed) == expected
    counts['inventory_rows'] = len(all_rows)
    counts['shortlist_rows'] = len(shortlist)

    # Figure 3 adds verified display fields; Figure 5 retains its numeric sources.
    for expected_path in sorted((ROOT / 'annotations/expected/figures').glob('*.csv')):
        actual_path = ROOT / 'annotations/results/figures' / expected_path.name
        assert read_csv(actual_path) == read_csv(expected_path), expected_path.name
    counts.update(check_figure3())
    figure5 = read_csv(ROOT / 'annotations/results/figures/Figure5_summary_source.csv')
    panel_c = [r for r in figure5 if r['panel'] == 'C']
    assert [(int(r['numerator']), int(r['denominator'])) for r in panel_c] == [
        (17, 203), (21, 115), (19, 75), (15, 49), (12, 34), (16, 22)]
    for extension in ('png', 'pdf'):
        assert digest(ROOT / f'figures/Figure5.{extension}') == digest(
            ROOT / f'annotations/results/figures/Figure5_annotation_diagnostics.{extension}')

    supporting = json.loads((ROOT / 'validation/SUPPORTING_TABLES_SOURCE.json').read_text(encoding='utf-8'))
    for name, source in supporting.items():
        assert digest(ROOT / name) == source['sha256'], f'Supporting table changed: {name}'
    counts['supporting_tables'] = len(supporting)

    # Guard the current-release layout and the user-facing navigation.
    repo = ROOT.parent
    for obsolete in ('manuscript_release_2026_06', 'analysis', '.zenodo.json',
                     'unified_pipeline_clean/unified_pipeline_run_public',
                     'unified_pipeline_clean/nif_hdk_scan_release_clean/logs'):
        assert not (repo / obsolete).exists(), f'Superseded material returned: {obsolete}'
    for obsolete in ('prior169_AllComposite.csv', 'prior169_TopFamilies.csv',
                     'released_AllComposite.csv', 'released_TopFamilies.csv', 'workbook_payload.json'):
        assert not (ROOT / 'annotations/results/current167' / obsolete).exists(), obsolete
    documents = [repo / 'README.md', ROOT / 'README.md',
                 ROOT / 'annotations/README.md', ROOT / 'models/README.md',
                 ROOT / 'models/validation/VALIDATION.md', ROOT / 'tables/README.md',
                 repo / 'unified_pipeline_clean/README.md',
                 *sorted((repo / 'unified_pipeline_clean/docs').glob('*.md')),
                 *sorted((ROOT / 'docs').glob('*.md'))]
    link_count = 0
    for document in documents:
        for target in re.findall(r'\]\(([^)]+)\)', document.read_text(encoding='utf-8')):
            if '://' in target or target.startswith(('mailto:', '#')):
                continue
            target = unquote(target.split('#')[0])
            assert (document.parent / target).exists(), f'Broken link in {document.name}: {target}'
            link_count += 1
    counts['local_documentation_links'] = link_count
    print(json.dumps({'status': 'PASS', 'checks': counts}, indent=2))


if __name__ == '__main__':
    main()
