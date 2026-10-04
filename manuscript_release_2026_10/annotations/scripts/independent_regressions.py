"""Reproduce the independent sensitivity checks from the released workbook.

Usage: python independent_checks.py [--workbook PATH]
Requires numpy, pandas, scipy, statsmodels, and openpyxl.
The default workbook is resolved relative to this script, not the working directory.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact

script_dir = Path(__file__).resolve().parent
workbook_name = 'proteomics_composite_family_evidence.xlsx'
default_candidates = [
    script_dir.parent / 'inputs' / 'released' / workbook_name,
    script_dir.parent / 'source_documents' / workbook_name,
    script_dir.parent / 'workbooks' / workbook_name,
    script_dir / workbook_name,
    script_dir.parent / workbook_name,
]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workbook', type=Path,
                    help='Workbook containing the AllComposite sheet. By default, find the supplied workbook relative to this script.')
parser.add_argument('--output-json', type=Path)
args = parser.parse_args()
workbook = args.workbook or next((p for p in default_candidates if p.is_file()), None)
if workbook is None:
    parser.error('No packaged workbook found; supply --workbook PATH.')
if not workbook.is_file():
    parser.error(f'Workbook does not exist: {workbook}')
ac = pd.read_excel(workbook, sheet_name='AllComposite')
T=lambda s:s.astype(str).str.lower().isin(['true','1','1.0'])
ms=ac[ac.candidate_set=='Model-Supported'].copy()
mm=ms[ms.n_literature_mapped_studies>=1].copy()
mm['up']=(mm.n_literature_active_phase_up_studies>=1).astype(int)
mm['strict']=T(mm.strict_unicellular_breadth_ge3u_ge1f).astype(int)
mm['nmap']=mm.n_literature_mapped_studies
mm['rank']=mm.consensus_rank_pct_mean
model_results=[]
def report(model):
    ci=model.conf_int()
    model_results.append({'formula':model.model.formula,'nobs':int(model.nobs),'parameters':{str(k):{'OR':float(np.exp(model.params[k])),'CI_low':float(np.exp(ci.loc[k,0])),'CI_high':float(np.exp(ci.loc[k,1])),'p':float(model.pvalues[k])} for k in model.params.index}})
    print(pd.DataFrame({'OR':np.exp(model.params),'CI_low':np.exp(ci[0]),'CI_high':np.exp(ci[1]),'p':model.pvalues}).to_string())
print('D5 original',len(mm),'complete rank',mm['rank'].notna().sum())
report(smf.logit('up ~ strict + nmap + rank',data=mm).fit(disp=0))
mm['rate']=mm.n_literature_active_phase_up_studies/mm.nmap
print('family-weighted and study-weighted rates')
for s,g in mm.groupby('strict'):
    print(s,len(g),g['rate'].mean(),g.n_literature_active_phase_up_studies.sum()/g.nmap.sum(),'nmapmean',g.nmap.mean())
print(mm.groupby(['strict','nmap']).up.agg(['sum','size','mean']).to_string())
print('D5 categorical opportunities sensitivity')
report(smf.logit('up ~ strict + C(nmap) + rank',data=mm).fit(disp=0))
cc=ac[ac.condensate_mapped_rows.fillna(0)>=1].copy()
cc['ms']=(cc.candidate_set=='Model-Supported').astype(int)
cc['top10']=T(cc.condensate_top10).astype(int)
cc['nrows']=cc.condensate_mapped_rows
print('D6 mapped sample and ranges',cc.groupby('candidate_set').n_member_genomes.agg(['size','min','max','median']).to_string())
print('D6 original without prevalence')
report(smf.logit('top10 ~ ms + np.log(nrows)',data=cc).fit(disp=0))
print('D6 original with prevalence')
report(smf.logit('top10 ~ ms + np.log(nrows) + np.log(n_member_genomes)',data=cc).fit(disp=0))
lo=max(cc.groupby('ms').n_member_genomes.min());hi=min(cc.groupby('ms').n_member_genomes.max())
ov=cc[cc.n_member_genomes.between(lo,hi)]
print('D6 common support',lo,hi,ov.groupby('ms').top10.agg(['sum','size','mean']).to_string())
print('D6 common support regression')
report(smf.logit('top10 ~ ms + np.log(nrows) + np.log(n_member_genomes)',data=ov).fit(disp=0))
print('D6 original strata')
for lo,hi in [(40,60),(60,80),(80,130)]:
    print(lo,hi,cc[cc.n_member_genomes.between(lo,hi,inclusive='left')].groupby('ms').top10.agg(['sum','size','mean']).to_string())

ABUNDANT=('clp|chaperon|dnak|groe|grol|dnaj|grpe|release factor|translation|elongation factor|typa|trna|nyn|ribosom|isopropylmalate|threonine synthase|ketol-acid|dihydroxy-acid|aminomutase|amino acid')
ms['abund']=ms['product'].str.lower().fillna('').str.contains(ABUNDANT).astype(int)
ms['tierA']=ms.priority_tier.fillna('').str.startswith('Tier A').astype(int)
ms['nmap']=ms.n_literature_mapped_studies.fillna(0).astype(int)
ms['strict']=T(ms.strict_unicellular_breadth_ge3u_ge1f).astype(int)
print('D7 cross-tab',ms.groupby(['tierA','nmap']).abund.agg(['sum','size','mean']).to_string())
print('D7 exploratory compositional association conditional on mapped studies')
report(smf.logit('abund ~ tierA + C(nmap)',data=ms).fit(disp=0))
strict=ms[ms.strict==1]
print('D7 strict only sample',len(strict))
report(smf.logit('abund ~ tierA + C(nmap)',data=strict).fit(disp=0))

if args.output_json:
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(model_results,indent=2),encoding="utf8")
