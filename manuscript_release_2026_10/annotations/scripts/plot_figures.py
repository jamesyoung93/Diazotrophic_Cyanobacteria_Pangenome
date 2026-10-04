"""Figures 3 and 5 from corrected scores with ubiquitous controls excluded from tiers.

Usage: python plot_corrected_score_figures_v2.py [--all-composite-csv PATH] [--output-dir DIR]
Requires: pandas, numpy, openpyxl, matplotlib. No inference uses a new taxonomy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import textwrap
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

FIGURE3_GROUPS = [
    ('Amino-acid biosynthesis', ('GF_00741','GF_01387','GF_00143','GF_01376','GF_00051')),
    ('Protein quality control', ('GF_00658','GF_00614','GF_00290','GF_02262','GF_01689','GF_01799')),
    ('Translation/RNA-associated products', ('GF_00623','GF_00771','GF_01542','GF_01479')),
    ('Tetrapyrrole-associated products', ('GF_00300','GF_01233')),
    ('Other products', ('GF_00303','GF_00122','GF_00536','GF_00916','GF_01449')),
]
FIGURE3_ORGANISMS = {
    'Crocosphaera subtropica ATCC 51142': 'Cr',
    'Trichodesmium erythraeum IMS101': 'Tr',
    'Nostoc punctiforme PCC 73102': 'No',
}


def figure3(a, dest):
    """Display the unchanged Tier A selection with matrix-checked carrier fractions."""
    t = a[a.priority_tier.str.startswith('Tier A', na=False)].copy()
    assert len(t)==22 and t.gene_family.is_unique
    assert truth(t.primary_bridge_broad_ge1u_ge1f).all()
    assert truth(t.strict_unicellular_breadth_ge3u_ge1f).all()
    grouping = {gf: (i, name) for i, (name, ids) in enumerate(FIGURE3_GROUPS, 1) for gf in ids}
    assert set(t.gene_family) == set(grouping)
    raw = Path(__file__).resolve().parent.parent/'inputs/raw'
    matrix = pd.read_csv(raw/'gene_family_matrix.csv', index_col=0)
    labels = pd.read_csv(raw/'complete_genomes_labeled.csv').set_index('assembly_accession')
    assert matrix.index.is_unique and labels.index.is_unique
    positive = truth(labels.loc[matrix.index, 'is_diazotroph'])
    assert matrix.shape == (426,2286) and int(positive.sum()) == 112
    present = matrix[t.gene_family].gt(0)
    counts, positive_counts = present.sum(), present.loc[positive].sum()
    t['short_product'] = t.gene_family.map(SHORT)
    t['display_group_index'] = t.gene_family.map(lambda gf: grouping[gf][0])
    t['display_group'] = t.gene_family.map(lambda gf: grouping[gf][1])
    t['panel_genomes'], t['panel_nifHDK_positive'] = len(matrix), int(positive.sum())
    t['carrier_count'] = t.gene_family.map(counts).astype(int)
    t['nifHDK_positive_carrier_count'] = t.gene_family.map(positive_counts).astype(int)
    t['carrier_purity_fraction'] = t.nifHDK_positive_carrier_count/t.carrier_count
    assert t.carrier_count.eq(t.verified_carrier_count).all()
    assert np.allclose(t.carrier_purity_fraction,t.diazotroph_pct_mean,rtol=0,atol=1e-15)
    assert t.short_product.notna().all()
    t['active_up_organism_codes'] = t.active_phase_up_organisms.map(
        lambda value: '+'.join(FIGURE3_ORGANISMS[name] for name in value.split('; ')))
    t['supplied_gene_top100_mapping_overlap'] = t.gene_top100_overlap.astype(int)
    assert t.loc[t.supplied_gene_top100_mapping_overlap.ne(0),'gene_family'].tolist() == ['GF_01387']
    t['marker_display'] = np.where(t.supplied_gene_top100_mapping_overlap.ne(0),'*','')
    t['marker_provenance'] = np.where(t.supplied_gene_top100_mapping_overlap.ne(0),
        'figures/Figure4_sources/gene_to_protein_top100_agreement.csv; '
        'supplied feature genome_GCF_000015665.1_CDS_1351; rank 71; mapped to GF_01387','')
    t['model_percentile_display'] = 100*t.consensus_rank_pct_mean
    t['condensate_percentile_display'] = 100*t.condensate_family_percentile
    t['carrier_purity_percent_display'] = t.carrier_purity_fraction.map(lambda p:f'{100*p:.1f}')
    t['corrected_score_display'] = t.scope_adjusted_story_score.map(lambda s:f'{s:.2f}')
    t = t.sort_values(['display_group_index','carrier_purity_fraction','gene_family'],ascending=[True,False,True])
    t['display_row'] = np.arange(1,len(t)+1)
    cols=['display_row','display_group_index','display_group','gene_family','product','short_product',
          'priority_tier','candidate_set','panel_genomes','panel_nifHDK_positive','carrier_count',
          'nifHDK_positive_carrier_count','carrier_purity_fraction','carrier_purity_percent_display',
          'n_literature_active_phase_up_studies','n_literature_mapped_studies',
          'active_up_organism_codes','active_phase_up_organisms',
          'condensate_family_percentile','condensate_percentile_display',
          'scope_adjusted_story_score','corrected_score_display',
          'supplied_gene_top100_mapping_overlap','marker_display','marker_provenance',
          'consensus_rank_pct_mean','model_percentile_display',
          'primary_bridge_broad_ge1u_ge1f','strict_unicellular_breadth_ge3u_ge1f']
    t[cols].to_csv(dest/'Figure3_source.csv',index=False,lineterminator='\n')
    fig = plt.figure(figsize=(6.2,7.2))
    fig.text(.025,.980,'Family prevalence and annotations',ha='left',va='top',
             fontsize=11,fontweight='bold',color=INK)
    fig.text(.025,.950,'All 22 Tier A families; product-annotation groups',ha='left',va='top',
             fontsize=8,color=GREY)
    ax=fig.add_axes([.020,.147,.960,.778]); ax.set_xlim(0,1); ax.set_ylim(28.0,-2.3); ax.axis('off')
    edges=[0,.31,.425,.515,.645,.80,.91,1]
    headers=['Product (abbreviated)','Family','Carriers\n(n/426)',
             'nifHDK+\namong carriers\n(%)','Active-up /\nmapped;\norganism(s)',
             'Condensate\npercentile','Score']
    for j,label in enumerate(headers):
        x=(edges[j]+edges[j+1])/2
        ax.text(edges[j]+.005 if j==0 else x,-1.40,label,
                ha='left' if j==0 else 'center',va='center',fontsize=6.5 if j==5 else 6.8,
                fontweight='bold',color=INK,linespacing=1.13)
    ax.plot([0,1],[-.45,-.45],color=INK,lw=.7)
    cmap=plt.get_cmap('Blues')
    y=0
    for name, _ in FIGURE3_GROUPS:
        ax.text(.005,y,name,ha='left',va='center',fontsize=7.5,fontweight='bold',color=BLUE)
        ax.plot([0,1],[y+.30,y+.30],color=GRID,lw=.55)
        y+=1
        group_rows=list(t[t.display_group.eq(name)].iterrows())
        labels=[textwrap.fill(r.short_product,width=28,break_long_words=False,break_on_hyphens=False)
                for _,r in group_rows]
        for row_index,(_,r) in enumerate(group_rows):
            ax.text(.005,y,labels[row_index],
                    ha='left',va='center',fontsize=7.1,color=INK,linespacing=1.05)
            family_x=(edges[1]+edges[2])/2
            ax.text(family_x,y,r.gene_family,ha='center',va='center',fontsize=6.8,color=GREY)
            if r.marker_display:
                ax.text(edges[2]-.008,y-.13,'*',ha='right',va='center',fontsize=7.4,color=INK)
            p=float(r.carrier_purity_fraction)
            ax.add_patch(Rectangle((edges[3]+.012,y-.39),edges[4]-edges[3]-.024,.78,
                                   facecolor=cmap(p),edgecolor='none'))
            nu,nr=int(r.n_literature_active_phase_up_studies),int(r.n_literature_mapped_studies)
            rank=r.condensate_percentile_display
            vals=[f'{r.carrier_count}/426',r.carrier_purity_percent_display,
                  f'{nu}/{nr}  {r.active_up_organism_codes}',
                  f'{rank:.1f}' if pd.notna(rank) else 'NA',r.corrected_score_display]
            for j,value in enumerate(vals,2):
                color='white' if j==3 and p>=.60 else INK
                ax.text((edges[j]+edges[j+1])/2,y,value,ha='center',va='center',
                        fontsize=7.1 if j!=6 else 6.9,color=color)
            consecutive_wrapped=(row_index+1<len(labels) and '\n' in labels[row_index]
                                 and '\n' in labels[row_index+1])
            y+=1.2 if consecutive_wrapped else 1
    assert abs(y-27.4)<1e-10
    ax.plot([0,1],[y-.50,y-.50],color=GRID,lw=.65)
    # The only color scale is the carrier percentage, fixed at 0-100%.
    key=fig.add_axes([.025,.134,.25,.012])
    gradient=np.linspace(0,1,256)[None,:]
    key.imshow(gradient,aspect='auto',cmap=cmap,extent=[0,100,0,1],vmin=0,vmax=1)
    key.set_yticks([]);key.set_xticks([0,50,100]);key.tick_params(axis='x',labelsize=6.5,length=2,pad=1)
    for spine in key.spines.values():spine.set_visible(False)
    fig.text(.295,.139,'nifHDK+ among carriers (%) | Panel reference: 26.3%',
             ha='left',va='center',fontsize=6.8,color=GREY)
    fig.text(.025,.095,'Cr: Crocosphaera; Tr: Trichodesmium; No: Nostoc (strains in caption).',fontsize=6.9,color=GREY)
    fig.text(.025,.073,'Condensate percentile: lower is a higher rank; NA: no mapped rank.',fontsize=6.9,color=GREY)
    fig.text(.025,.051,'* Also mapped from the supplied gene-level top-100 list.',fontsize=6.9,color=GREY)
    fig.text(.025,.025,'Display groups and carrier percentages are descriptive; tier selection is unchanged.',
             fontsize=6.9,color=INK)
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
        axb.plot([],[],color=col,marker='o' if cs=='Model-Supported' else 's',ms=3,lw=.8,label={'Model-Supported':'Residual association','Highly Pure':'High purity'}[label])
    axb.legend(frameon=False,fontsize=6.3,loc='upper left',borderaxespad=.15,handlelength=1.1,labelspacing=.3)
    axb.text(-.22,-.31,'k/n:',transform=axb.transAxes,ha='left',va='top',fontsize=6.2,color=GREY)
    axb.text(0,-.55,'Limited overlap: high purity n = 11 at 80+ carriers',
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

CAPTIONS='''Figure 3. Family prevalence and annotations for all 22 Tier A residual-association families. Rows use heuristic product-annotation groups and, within each group, decreasing carrier purity followed by family ID. These display groups do not change scores or tier membership. Carriers gives the number of panel genomes containing the family. Shading shows the percentage of those carriers labeled nifHDK-positive, not coverage of positive genomes or evidence strength; 112 of 426 panel genomes are positive (26.3%). Active-up gives studies with an active-phase-up annotation over studies with a mapped protein, followed by the contributing organisms: Cr, Crocosphaera subtropica ATCC 51142; Tr, Trichodesmium erythraeum IMS101; No, Nostoc punctiforme PCC 73102. Study coverage and response rules differ, so these are descriptive counts. Condensate percentile summarizes the best mapped protein per family; lower values indicate higher ranks. NA means no mapped rank. The prioritization score selected Tier A at ≥6.0 after the stated exclusions and is shown for reference. An asterisk marks a match in the supplied gene-level top-100 mapping. Product labels are abbreviated (DH, dehydrogenase; PP, diphosphate); full labels accompany the source data. These annotations guide experiments without independently validating candidates or establishing genetic requirement.

Figure 5. Detection opportunity and family prevalence constrain annotation summaries.

(A) Among residual-association families mapped in at least one external proteomics study, the proportion with any active-phase-up annotation is shown by mapped-study count and the strict unicellular-breadth flag (>=3 strict-proxy unicellular and >=1 filamentous-proxy diazotroph carrier). The reported 22.3% values are means of each family's active-up-study count divided by its mapped-study count, not pooled per-study probabilities. (B) The proportion of condensate-mapped candidate families in the top 10% of mapped family ranks is shown by inventory and carrier-genome count. Bins are half-open [40,60), [60,80), [80,130), and [130,infinity), matching the diagnostic stratification through 130 genomes. All nonempty strata are shown, including the 11 high-purity families at 80-129 carriers; no high-purity family has 130 or more carriers. (C) Matches to a product-name pattern covering proteostasis, translation, and amino-acid metabolism are shown across mapped-study counts in all 476 residual-association families. Diagnostic strata retain all 476 residual-association families, including three ubiquitous controls excluded from shortlist eligibility. Under the prioritization score, the pattern matches 16 of 22 Tier A families (72.7%). Tier A is a selected subset overlapping those strata, not an independent comparison group. Under each panel, k/n gives the numerator and denominator in matching series order and color; 0/0 denotes an empty stratum with no plotted estimate. Error bars are descriptive 95% Wilson intervals for family proportions and do not account for phylogenetic or other dependencies among families. Product-name matching is a heuristic diagnostic, not a validated functional classification. These patterns motivate detection- and prevalence-aware interpretation; they do not establish absence of biological roles for individual candidates.'''

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
    (args.output_dir/'figure_captions.txt').write_text(CAPTIONS,encoding='utf-8',newline='\n')
    caption3,caption5=CAPTIONS.strip().split('\n\n',1)
    (args.output_dir/'Figure3_caption.txt').write_text(caption3+'\n',encoding='utf-8',newline='\n')
    (args.output_dir/'Figure5_caption.txt').write_text(caption5+'\n',encoding='utf-8')
    shutil.copyfile(args.output_dir/'Figure3_TierA_evidence.png',args.output_dir/'Figure3_all_TierA.png')
    shutil.copyfile(args.output_dir/'Figure5_annotation_diagnostics.png',args.output_dir/'Figure5_diagnostics.png')
    manifest={'status':'Proposed revision for author review',
              'all_composite_csv':args.all_composite_csv.name,
              'source_sha256':hashlib.sha256(args.all_composite_csv.read_bytes()).hexdigest(),
              'n_inventory_families':len(a),'n_figure3_families':len(t),
              'figure3_size_inches':[6.2,7.2],'figure5_size_inches':[6.2,5.95],
              'png_dpi':300,'diagnostic_notes':diag,
              'figure3_display_order':'Fixed heuristic product-annotation groups; decreasing carrier purity, then family ID within each group. Does not alter selection.',
              'figure3_shading':'nifHDK-positive carriers / all carriers; fixed 0-100 percent scale. Counts recomputed from the frozen matrix and labels.',
              'figure3_marker':'Also mapped from the supplied gene-level top-100 list; not independent validation.',
              'model_percentile_definition':'Mean absolute-importance percentile across available model/run entries; retained in source data, not displayed in Figure 3.',
              'condensate_percentile_definition':'Condensate rank divided by mapped-family count; lower is better; displayed 0-100.',
              'model_direction_source':'unified_pipeline_clean/nif_downstream_code/11_build_fox_gene_report.py:194-195,389,438'}
    (args.output_dir/'figure_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8',newline='\n')
    print(json.dumps(manifest,indent=2))

if __name__=='__main__': main()

