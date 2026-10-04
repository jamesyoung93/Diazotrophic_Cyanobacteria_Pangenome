"""Replay the annotation/scoring audit using only packaged inputs (no network)."""
from pathlib import Path
import argparse, hashlib, importlib, importlib.metadata, json, os, subprocess, sys, time

P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--out',type=Path,default=Path('results'),help='Output folder relative to this package, or an explicit absolute folder.')
args=ap.parse_args()
out=args.out if args.out.is_absolute() else P/args.out
out.mkdir(parents=True,exist_ok=True)
logs=out/'logs';logs.mkdir(exist_ok=True)
def arg(p):return os.path.relpath(p,P)
manifest=P/'PACKAGE_MANIFEST.json'
if not manifest.exists():raise SystemExit('Missing PACKAGE_MANIFEST.json')
hashes=json.loads(manifest.read_text(encoding='utf8'))['files']
for name,expected in hashes.items():
    assert hashlib.sha256((P/name).read_bytes()).hexdigest()==expected,f'Package file changed: {name}'
modules={'numpy':'numpy','pandas':'pandas','scipy':'scipy','statsmodels':'statsmodels','openpyxl':'openpyxl','matplotlib':'matplotlib','pillow':'PIL'}
runtime={'python':sys.version.split()[0],'packages':{p:{'module_version':importlib.import_module(m).__version__,'distribution_metadata_version':importlib.metadata.version(p)} for p,m in modules.items()}}
(out/'runtime_versions.json').write_text(json.dumps(runtime,indent=2),encoding='utf8')
workbook='inputs/released/proteomics_composite_family_evidence.xlsx'
steps=[
 ('01_reconstruction',['scripts/reconstruct_annotations.py','--inputs','inputs','--output-dir',arg(out/'reconstruction')]),
 ('02_current167',['scripts/score_current.py','--workbook',workbook,'--matrix','inputs/raw/gene_family_matrix.csv','--independent-reconstruction',arg(out/'reconstruction/score_family_reconstruction.json'),'--output-dir',arg(out/'current167')]),
 ('03_historical_diagnostics',['scripts/historical_diagnostics.py','--workbook',workbook,'--matrix','inputs/raw/gene_family_matrix.csv','--labels','inputs/raw/complete_genomes_labeled.csv']),
 ('04_independent_regressions',['scripts/independent_regressions.py','--workbook',workbook,'--output-json',arg(out/'diagnostics/independent_regressions.json')]),
 ('05_mapping_sensitivity',['scripts/mapping_sensitivity.py','--corrected-csv',arg(out/'current167/corrected_AllComposite.csv'),'--output-dir',arg(out/'mapping_sensitivity')]),
 ('06_model_rank_from_saved_outputs',['rank_audit/run_rank_audit.py','--out',arg(out/'model_rank'),'--verify']),
 ('07_current_shortlist_rank',['scripts/current_shortlist_rank.py','--all-composite-csv',arg(out/'current167/corrected_AllComposite.csv'),'--rank-dir',arg(out/'model_rank'),'--output-dir',arg(out/'rank_sensitivity')]),
 ('08_figures',['scripts/plot_figures.py','--all-composite-csv',arg(out/'current167/corrected_AllComposite.csv'),'--output-dir',arg(out/'figures')]),
 ('09_verification',['scripts/verify_reproduction.py','--results',arg(out)]),
]
env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1'
timings=[]
for name,command in steps:
    print(f'Running {name} ...',flush=True)
    start=time.monotonic()
    proc=subprocess.run([sys.executable,*command],cwd=P,env=env,capture_output=True,text=True,encoding='utf8')
    (logs/(name+'.txt')).write_text(proc.stdout+('\nSTDERR:\n'+proc.stderr if proc.stderr else ''),encoding='utf8')
    timings.append({'step':name,'seconds':round(time.monotonic()-start,3),'returncode':proc.returncode,'command':['python',*command]})
    if proc.returncode:
        (out/'run_manifest.json').write_text(json.dumps(timings,indent=2),encoding='utf8')
        print(proc.stdout[-2000:]+proc.stderr[-4000:],file=sys.stderr)
        raise SystemExit(f'{name} failed; see its log.')
(out/'run_manifest.json').write_text(json.dumps(timings,indent=2),encoding='utf8')
output_hashes={p.relative_to(out).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_HASHES.json'}
(out/'OUTPUT_HASHES.json').write_text(json.dumps({'algorithm':'SHA-256','files':output_hashes},indent=2),encoding='utf8')
print('PASS: all annotation, score, regression, mapping, rank and figure-source checks passed.')
