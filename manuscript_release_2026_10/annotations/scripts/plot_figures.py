"""Figures 3 and 5 from corrected scores with ubiquitous controls excluded from tiers.

Usage: python plot_corrected_score_figures_v2.py [--all-composite-csv PATH] [--output-dir DIR]
Requires: pandas, numpy, openpyxl, matplotlib. No inference uses a new taxonomy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

BLUE, ORANGE, INK, GREY = '#246AA4', '#BC5628', '#172536', '#5C6673'
LIGHT, GRID = '#F2F5F8', '#DBE1E7'
ABUNDANT = ('clp|chaperon|dnak|groe|grol|dnaj|grpe|release factor|translation|elongation factor|'
            'typa|trna|nyn|ribosom|isopropylmalate|threonine synthase|ketol-acid|dihydroxy-acid|'
            'aminomutase|amino acid')

SHORT = {
 'GF_00658': 'ClpX protease ATPase',
 'GF_00143': 'Ketol-acid reductoisomerase',
 'GF_00623': 'Peptide release factor 3',
 'GF_01387': '2-Isopropylmalate synthase',
 'GF_00051': 'Threonine synthase',
 'GF_00122': 'Glutathione / class III alcohol DH',
 'GF_00303': 'Family 2 glycosyltransferase',
 'GF_00771': 'NYN domain protein',
 'GF_01233': 'Glutamate-1-semialdehyde aminomutase',
 'GF_01542': 'tRNA methylthiotransferase MiaB',
 'GF_01689': 'Clp protease subunit',
 'GF_00614': 'Chaperone DnaK',
 'GF_00916': 'NAD(P)H-quinone oxidoreductase H',
 'GF_00741': '3-Isopropylmalate dehydrogenase',
 'GF_00360': 'Inorganic diphosphatase',
 'GF_00079': 'Clp protease subunit',
 'GF_00104': '15-cis-Phytoene desaturase',
 'GF_00136': 'Protochlorophyllide reductase B',
 'GF_00622': 'Amino acid ABC transporter ATPase',
 'GF_00624': 'NADP-isocitrate dehydrogenase',
 'GF_01799': 'Chaperonin GroEL',
 'GF_01816': 'NAD(P)H-quinone oxidoreductase K',
 'GF_01818': 'Glucose-6-phosphate dehydrogenase',
 'GF_01860': 'Hypothetical protein',
 'GF_00290': 'Clp protease ATPase',
 'GF_00300': 'Mg-protoporphyrin ester cyclase',
 'GF_00536': 'Gas vesicle protein GvpA',
 'GF_02262': 'Co-chaperone GroES',
 'GF_01376': '3-Isopropylmalate dehydratase, large',
 'GF_01449': 'Hydroxy-methylbutenyl-PP synthase',
 'GF_01479': 'Translational GTPase TypA',
}

def truth(s):
    return s.astype(str).str.lower().isin(['true', '1', '1.0'])

def wilson(k, n):
    if not n:
        return float('nan'), float('nan')
    z = 1.959963984540054
    p = k / n
    centre = (p + z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return max(0,centre-half), min(1,centre+half)

def save(fig, dest, stem):
    # Fixed physical dimensions: no tight bounding-box crop to alter print size.
    fig.savefig(dest / f'{stem}.png', dpi=300, facecolor='white')
    fig.savefig(dest / f'{stem}_printwidth96.png', dpi=96, facecolor='white')
    fig.savefig(dest / f'{stem}.pdf', facecolor='white',
                metadata={'Title': stem, 'Author': 'Beyond nif manuscript revision'})
    plt.close(fig)

def figure3(a, dest):
    t = a[a.priority_tier.str.startswith('Tier A', na=False)].copy()
    t = t.sort_values(['scope_adjusted_story_score','gene_family'], ascending=[False,True])
    assert len(t)==22 and t.gene_family.is_unique
    assert truth(t.primary_bridge_broad_ge1u_ge1f).all()
    assert truth(t.strict_unicellular_breadth_ge3u_ge1f).all()
    t['short_product'] = t.gene_family.map(SHORT)
    assert t.short_product.notna().all()
    t['model_percentile_display'] = 100*t.consensus_rank_pct_mean
    t['condensate_percentile_display'] = 100*t.condensate_family_percentile
    cols=['gene_family','product','short_product','priority_tier','candidate_set',
          'scope_adjusted_story_score','consensus_rank_pct_mean','model_percentile_display',
          'n_literature_active_phase_up_studies','n_literature_mapped_studies',
          'condensate_family_percentile','condensate_percentile_display',
          'primary_bridge_broad_ge1u_ge1f','strict_unicellular_breadth_ge3u_ge1f']
    t[cols].to_csv(dest/'Figure3_source.csv',index=False)
    fig = plt.figure(figsize=(6.2,4.85))
    fig.text(.025,.975,'All 22 Tier A Model-Supported families',ha='left',va='top',
             fontsize=10.6,fontweight='bold',color=INK)
    fig.text(.025,.940,'Corrected score; rows ordered by score, then family ID',ha='left',va='top',
             fontsize=7.4,color=GREY)
    ax=fig.add_axes([.018,.110,.966,.777]); ax.set_xlim(0,1); ax.set_ylim(22,-3); ax.axis('off')
    edges=[0,.434,.555,.645,.79,.89,1]
    # Labels occupy the first two columns; values use four independent scales.
    headers=['Product (abbreviated)','Family','Corrected\nscore','Released model\nimportance\npercentile','Active-up /\nmapped','Condensate\npercentile']
    for j,label in enumerate(headers):
        x=(edges[j]+edges[j+1])/2
        ax.text(edges[j]+.006 if j==0 else x,-1.50,label,ha='left' if j==0 else 'center',
                va='center',fontsize=6.5,fontweight='bold',color=INK,linespacing=1.1)
    ax.plot([0,1],[-.46,-.46],color=INK,lw=.65)
    cmap=plt.get_cmap('Blues')
    for i,(_,r) in enumerate(t.iterrows()):
        if i%2==0:
            ax.add_patch(Rectangle((0,i-.45),1,.94,facecolor=LIGHT,edgecolor='none',zorder=-5))
        ax.text(.006,i,r.short_product,ha='left',va='center',fontsize=6.75,color=INK)
        ax.text((edges[1]+edges[2])/2,i,r.gene_family,ha='center',va='center',fontsize=6.35,color=GREY)
        nr=int(r.n_literature_mapped_studies); nu=int(r.n_literature_active_phase_up_studies)
        p=r.condensate_family_percentile
        vals=[(f'{r.scope_adjusted_story_score:.2f}',r.scope_adjusted_story_score/9),
              (f'{100*r.consensus_rank_pct_mean:.1f}',r.consensus_rank_pct_mean),
              (f'{nu} / {nr}',nu/nr if nr else None),
              (f'{100*p:.1f}' if pd.notna(p) else 'NA',1-p if pd.notna(p) else None)]
        for j,(value,level) in enumerate(vals,start=2):
            x0,x1=edges[j]+.006,edges[j+1]-.006
            if level is None:
                fill='#F0F0F0'
            else:
                fill=cmap(.02+.60*min(max(float(level),0),1))
            ax.add_patch(Rectangle((x0,i-.405),x1-x0,.81,facecolor=fill,edgecolor='none'))
            ax.text((x0+x1)/2,i,value,ha='center',va='center',fontsize=6.8,color=INK)
    fig.text(.025,.070,'Model importance percentile: higher is better. Condensate rank percentile: lower is better.',
             fontsize=6.8,color=GREY,ha='left')
    fig.text(.025,.045,'Darker cells follow these directions, higher score, or up/mapped fraction. NA: no condensate rank.',
             fontsize=6.65,color=GREY,ha='left')
    fig.text(.025,.020,'All 22 pass both scored morphotype flags. Tier A is a heuristic score band.',
             fontsize=6.75,color=INK,ha='left')
    save(fig,dest,'Figure3_TierA_evidence')
    return t

def axis_style(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['bottom','left']].set_color('#9CA7B2')
    ax.tick_params(colors=GREY,labelsize=7,width=.6,length=2.5)
    ax.yaxis.grid(True,color=GRID,lw=.6); ax.set_axisbelow(True)

def add_point(ax,x,k,n,color,marker='o'):
    if n==0:
        return
    low,high=wilson(k,n); y=100*k/n
    ax.errorbar(x,y,yerr=[[y-100*low],[100*high-y]],fmt=marker,color=color,
                markersize=4.3,lw=1,capsize=2,markeredgewidth=.7,markeredgecolor='white',zorder=4)

def summary_record(panel,group,binlabel,k,n):
    lo,hi=wilson(k,n)
    return dict(panel=panel,group=group,stratum=binlabel,numerator=int(k),denominator=int(n),
                fraction=k/n if n else np.nan,ci95_low=lo,ci95_high=hi)

def figure5(a,t,dest):
    ms=a[a.candidate_set.eq('Model-Supported')].copy()
    ms['strict_breadth']=truth(ms.strict_unicellular_breadth_ge3u_ge1f)
    ms['any_up']=ms.n_literature_active_phase_up_studies.ge(1)
    ms['product_pattern_group']=ms['product'].str.lower().fillna('').str.contains(ABUNDANT,regex=True)
    ms['mapped_studies']=ms.n_literature_mapped_studies.fillna(0).astype(int)
    mm=ms[ms.mapped_studies.ge(1)].copy()
    mm['up_fraction_per_family']=mm.n_literature_active_phase_up_studies/mm.mapped_studies
    cc=a[a.condensate_mapped_rows.fillna(0).ge(1)].copy()
    cc['top10']=truth(cc.condensate_top10)
    cc['prevalence_bin']=pd.cut(cc.n_member_genomes,[40,60,80,130,np.inf],right=False,
                                labels=['40-<60','60-<80','80-<130','130+'])
    assert cc.prevalence_bin.notna().all()
    fig=plt.figure(figsize=(6.2,5.95))
    fig.text(.025,.98,'Annotation summaries by detection and prevalence',ha='left',va='top',
             fontsize=10.6,fontweight='bold',color=INK)
    axa=fig.add_axes([.09,.605,.38,.275]); axb=fig.add_axes([.595,.605,.38,.275])
    axc=fig.add_axes([.09,.155,.885,.235])
    summary=[]
    for ax in [axa,axb,axc]: axis_style(ax)
    axa.set_title('A  Active-phase-up annotation',loc='left',fontsize=8.5,fontweight='bold',pad=12)
    axa.set_ylim(0,103); axa.set_yticks([0,25,50,75,100]); axa.set_xlim(.5,4.5)
    axa.set_xticks(range(1,5)); axa.set_xlabel('Mapped proteomics studies',fontsize=7.4,labelpad=2)
    axa.set_ylabel('Families with any active-up (%)',fontsize=7.2,labelpad=3)
    for flag,col,offset,label,row_y in [(True,BLUE,-.07,'Strict breadth',-.31),
                                       (False,ORANGE,.07,'Other',-.41)]:
        vals=[]
        for nmap in range(1,5):
            g=mm[(mm.strict_breadth==flag)&(mm.mapped_studies==nmap)]
            k,n=int(g.any_up.sum()),len(g)
            summary.append(summary_record('A',label,str(nmap),k,n))
            if n: add_point(axa,nmap+offset,k,n,col,'o' if flag else 's'); vals.append((nmap+offset,100*k/n))
            axa.text(nmap,row_y,f'{k}/{n}' if n else '0/0',transform=axa.get_xaxis_transform(),
                     ha='center',va='top',fontsize=6.35,color=col)
        if vals: axa.plot(*zip(*vals),color=col,lw=.8,alpha=.7,zorder=2)
        axa.plot([],[],color=col,marker='o' if flag else 's',ms=3,lw=.8,label=label)
    axa.legend(frameon=False,fontsize=6.4,loc='upper left',borderaxespad=0.15,handlelength=1.2,labelspacing=.3)
    axa.text(-.22,-.31,'k/n:',transform=axa.transAxes,ha='left',va='top',fontsize=6.2,color=GREY)
    broad=100*mm.loc[mm.strict_breadth,'up_fraction_per_family'].mean()
    other=100*mm.loc[~mm.strict_breadth,'up_fraction_per_family'].mean()
    axa.text(0,-.55,f'Mean per-family up/mapped: {broad:.1f}% vs {other:.1f}%',
             transform=axa.transAxes,fontsize=6.75,color=INK,ha='left')

    axb.set_title('B  Condensate rank by prevalence',loc='left',fontsize=8.5,fontweight='bold',pad=12)
    axb.set_ylim(0,63); axb.set_yticks([0,20,40,60]); axb.set_xlim(-.5,3.5)
    axb.set_xticks(range(4)); axb.set_xticklabels(['40-<60','60-<80','80-<130','130+'],fontsize=6.5)
    axb.set_xlabel('Carrier genomes per family',fontsize=7.4,labelpad=2)
    axb.set_ylabel('Families in top 10% rank (%)',fontsize=7.2,labelpad=3)
    for cs,col,offset,label,row_y in [('Model-Supported',BLUE,-.07,'Model-Supported',-.31),
                                     ('Highly Pure',ORANGE,.07,'Highly Pure',-.41)]:
        vals=[]
        for i,binlabel in enumerate(cc.prevalence_bin.cat.categories):
            g=cc[(cc.candidate_set==cs)&(cc.prevalence_bin==binlabel)]
            k,n=int(g.top10.sum()),len(g)
            summary.append(summary_record('B',label,str(binlabel),k,n))
            if n: add_point(axb,i+offset,k,n,col,'o' if cs=='Model-Supported' else 's');vals.append((i+offset,100*k/n))
            axb.text(i,row_y,f'{k}/{n}' if n else '0/0',transform=axb.get_xaxis_transform(),
                     ha='center',va='top',fontsize=6.35,color=col)
        if vals: axb.plot(*zip(*vals),color=col,lw=.8,alpha=.7,zorder=2)
        axb.plot([],[],color=col,marker='o' if cs=='Model-Supported' else 's',ms=3,lw=.8,label=label)
    axb.legend(frameon=False,fontsize=6.3,loc='upper left',borderaxespad=.15,handlelength=1.1,labelspacing=.3)
    axb.text(-.22,-.31,'k/n:',transform=axb.transAxes,ha='left',va='top',fontsize=6.2,color=GREY)
    axb.text(0,-.55,'Limited overlap: Highly Pure n = 11 at 80+ carriers',
             transform=axb.transAxes,fontsize=6.65,color=INK,ha='left')

    axc.set_title('C  Product-name matches across detection strata',loc='left',fontsize=8.5,fontweight='bold',pad=12)
    axc.set_xlim(-.45,5.5);axc.set_ylim(0,100);axc.set_yticks([0,20,40,60,80,100])
    axc.set_xticks(range(6));axc.set_xticklabels(['0','1','2','3','4','Tier A'],fontsize=7.2)
    axc.set_xlabel('Mapped proteomics studies (0-4), with Tier A shown separately',fontsize=7.4,labelpad=3)
    axc.set_ylabel('Product-pattern group (%)',fontsize=7.2,labelpad=3)
    vals=[]
    for i in range(6):
        g=ms[ms.mapped_studies.eq(i)] if i<5 else ms[ms.priority_tier.str.startswith('Tier A',na=False)]
        k,n=int(g.product_pattern_group.sum()),len(g)
        summary.append(summary_record('C','Model-Supported',str(i) if i<5 else 'Tier A',k,n))
        add_point(axc,i,k,n,BLUE if i<5 else INK,'o' if i<5 else 'D')
        if i<5: vals.append((i,100*k/n))
        axc.text(i,-.29,f'{k}/{n}',transform=axc.get_xaxis_transform(),ha='center',va='top',fontsize=6.6,color=GREY)
    axc.plot(*zip(*vals),color=BLUE,lw=.85,zorder=2)
    axc.axvline(4.5,color=GRID,lw=1,ls='--')
    axc.text(.02,.89,'Name pattern: proteostasis / translation / amino-acid metabolism',transform=axc.transAxes,fontsize=6.6,color=GREY)
    fig.text(.025,.042,'Points: family proportions. Bars: descriptive 95% Wilson intervals; k/n counts appear below each panel.',
             fontsize=6.7,color=GREY)
    fig.text(.025,.022,'Families are not independent replicates. Product-name groups are heuristic; Tier A overlaps detection strata.',
             fontsize=6.7,color=GREY)
    save(fig,dest,'Figure5_annotation_diagnostics')
    pd.DataFrame(summary).to_csv(dest/'Figure5_summary_source.csv',index=False)
    ms[['gene_family','product','priority_tier','mapped_studies','n_literature_active_phase_up_studies',
        'strict_breadth','product_pattern_group','ubiquitous_control_flag']].to_csv(dest/'Figure5_AC_family_source.csv',index=False)
    cc[['gene_family','candidate_set','n_member_genomes','condensate_mapped_rows','condensate_family_percentile',
        'top10','prevalence_bin']].to_csv(dest/'Figure5_B_family_source.csv',index=False)
    return {'A_mean_per_family_up_fraction':{'strict':broad/100,'other':other/100},
            'mapped_condensate_families':int(len(cc)), 'pattern':ABUNDANT}

CAPTIONS='''Figure 3. Annotation matrix for all 22 Tier A Model-Supported families under the corrected score.
Rows are all 22 Tier A families (corrected adjusted score >=6.0 after ubiquitous-control, nif-keyword and housekeeping exclusions), ordered by decreasing score then family ID. The correction removes unverified cross-atlas identifier evidence and its imported HGT penalty while retaining the Highly Pure inventory offset; all other weights and thresholds are unchanged. Columns show corrected score, released mean model-importance percentile across available model/run entries, active-phase-up/mapped-study counts, and condensate rank percentile. The score retains a flat weight for historical inventory membership and does not incorporate corrected-model importance. Percentiles use a 0-100 scale: higher released model importance is better; lower condensate rank is better. Darker cells follow these directions, higher score or up/mapped fraction; shading is scaled separately by column. NA means no mapped condensate rank. Product labels are abbreviated (DH, dehydrogenase; PP, diphosphate); complete labels are in the source data. All 22 satisfy both scored morphotype flags. Tier A is a heuristic score band. These annotations guide experiments without independently validating candidates or establishing genetic requirement.

Figure 5. Detection opportunity and family prevalence constrain annotation summaries.
(A) Among Model-Supported families mapped in at least one external proteomics study, the proportion with any active-phase-up annotation is shown by mapped-study count and the strict unicellular-breadth flag (>=3 strict-proxy unicellular and >=1 filamentous-proxy diazotroph carrier). The reported 22.3% values are means of each family's active-up-study count divided by its mapped-study count, not pooled per-study probabilities. (B) The proportion of condensate-mapped candidate families in the top 10% of mapped family ranks is shown by inventory and carrier-genome count. Bins are half-open [40,60), [60,80), [80,130), and [130,infinity), matching the diagnostic stratification through 130 genomes. All nonempty strata are shown, including the 11 Highly Pure families at 80-129 carriers; no Highly Pure family has 130 or more carriers. (C) Matches to a product-name pattern covering proteostasis, translation, and amino-acid metabolism are shown across mapped-study counts in all 476 Model-Supported families. Diagnostic strata retain the frozen 476 Model-Supported families, including three ubiquitous controls excluded from shortlist eligibility. Under the corrected score, the pattern matches 16 of 22 Tier A families (72.7%). Tier A is a selected subset overlapping those strata, not an independent comparison group. Under each panel, k/n gives the numerator and denominator in matching series order and color; 0/0 denotes an empty stratum with no plotted estimate. Error bars are descriptive 95% Wilson intervals for family proportions and do not account for phylogenetic or other dependencies among families. Product-name matching is a heuristic diagnostic, not a validated functional classification. These patterns motivate detection- and prevalence-aware interpretation; they do not establish absence of biological roles for individual candidates.
'''

def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument('--all-composite-csv',type=Path,default=here.parent/'results/current167/corrected_AllComposite.csv')
    ap.add_argument('--output-dir',type=Path,default=here.parent/'results/figures')
    args=ap.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'text.color':INK,
        'axes.labelcolor':INK,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.65})
    a=pd.read_csv(args.all_composite_csv)
    assert a.related_atlas_score.eq(0).all() and not truth(a.hgt_passenger_flag).any()
    t=figure3(a,args.output_dir); diag=figure5(a,t,args.output_dir)
    (args.output_dir/'figure_captions.txt').write_text(CAPTIONS,encoding='utf-8')
    caption3,caption5=CAPTIONS.strip().split('\n\n',1)
    (args.output_dir/'Figure3_caption.txt').write_text(caption3+'\n',encoding='utf-8')
    (args.output_dir/'Figure5_caption.txt').write_text(caption5+'\n',encoding='utf-8')
    shutil.copyfile(args.output_dir/'Figure3_TierA_evidence.png',args.output_dir/'Figure3_all_TierA.png')
    shutil.copyfile(args.output_dir/'Figure5_annotation_diagnostics.png',args.output_dir/'Figure5_diagnostics.png')
    manifest={'status':'Proposed revision for author review',
              'all_composite_csv':args.all_composite_csv.name,
              'source_sha256':hashlib.sha256(args.all_composite_csv.read_bytes()).hexdigest(),
              'n_inventory_families':len(a),'n_figure3_families':len(t),
              'figure3_size_inches':[6.2,4.85],'figure5_size_inches':[6.2,5.95],
              'png_dpi':300,'diagnostic_notes':diag,
              'model_percentile_definition':'Mean absolute-importance percentile across available model/run entries; higher is better; displayed 0-100.',
              'condensate_percentile_definition':'Condensate rank divided by mapped-family count; lower is better; displayed 0-100.',
              'model_direction_source':'unified_pipeline_clean/nif_downstream_code/11_build_fox_gene_report.py:194-195,389,438'}
    (args.output_dir/'figure_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2))

if __name__=='__main__': main()

