"""Check the frozen release and manuscript-facing summaries without dependencies."""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
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

    # Figure relabeling must not alter any numeric source table.
    for expected_path in sorted((ROOT / 'annotations/expected/figures').glob('*.csv')):
        actual_path = ROOT / 'annotations/results/figures' / expected_path.name
        assert read_csv(actual_path) == read_csv(expected_path), expected_path.name
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
