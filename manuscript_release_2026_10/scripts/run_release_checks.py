"""Run release integrity, annotation replay and saved-model verification in isolation."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out', required=True, type=Path,
                    help='New output directory outside the immutable release.')
args = parser.parse_args()
out = args.out.resolve()
if out == ROOT or ROOT in out.parents:
    parser.error('--out must be outside the immutable release')
if out.exists():
    parser.error('--out must be a new directory to preserve earlier verification')
out.mkdir(parents=True)
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONIOENCODING='utf-8',
           MPLBACKEND='Agg', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
records = []


def run(label, arguments):
    print(f'Running {label}', flush=True)
    with (out / f'{label}.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, '-B', '-X', 'utf8', *map(str, arguments)],
                                env=env, stdout=log, stderr=subprocess.STDOUT)
    records.append({'check': label, 'returncode': result.returncode})
    (out / 'verification.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    if result.returncode:
        print((out / f'{label}.log').read_text(encoding='utf-8')[-6000:])
        raise SystemExit(result.returncode)
    print(f'PASS: {label}', flush=True)


run('release_integrity_before', [ROOT / 'scripts/verify_release.py'])
run('annotation_replay', [ROOT / 'annotations/run_all.py', '--out', out / 'annotations'])
# The legacy verifier writes metric-check reports under --out. Verify a copy so
# those reports cannot overwrite the frozen reference files.
shutil.copytree(ROOT / 'models/reference', out / 'model_reference_copy')
run('saved_model_metrics_and_sensitivities', [ROOT / 'models/run_audit.py', '--mode', 'verify',
                                       '--out', out / 'model_reference_copy'])
run('release_integrity_after', [ROOT / 'scripts/verify_release.py'])
print('PASS: release integrity, complete annotation replay and saved model checks including sensitivities.')
print('No model refits were performed; strict clean-backend refit limitations remain documented.')
