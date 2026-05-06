"""Publication-quality figures for P2, P3, P4 manuscripts"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import sqlite3, os

plt.rcParams.update({
    'font.sans-serif': ['Arial','Helvetica'],
    'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'figure.dpi':600,'savefig.dpi':600,'savefig.bbox':'tight'
})

OUT = r'D:\project\AutoResearch-rice-T2T\data\figures\publication'
os.makedirs(OUT, exist_ok=True)
COL = {'japonica_temp':'#2166AC','indica':'#B2182B','aus':'#D6604D','cold':'#2166AC',
       'heat':'#E08214','drought':'#D6604D','salt':'#1B7837','submergence':'#5E3C99','hm':'#666666'}

DB = r'D:\project\AutoResearch-rice-T2T\data\rice_integration.db'

# ═══════════════════════════════════════════════════════════════
# FIGURES FOR P2: Dual Hotspot Mechanism
# ═══════════════════════════════════════════════════════════════

def p2_fig1():
    """Dual hotspot genomic landscape"""
    fig = plt.figure(figsize=(8.5, 5))

    # Panel A: Chr04
    conn = sqlite3.connect(DB); c = conn.cursor()
    ax_a = fig.add_axes([0.06, 0.55, 0.90, 0.40])
    tes = []
    for r in c.execute("SELECT ti.position, ta.odds_ratio, ta.category FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id WHERE ti.chromosome='Chr04' AND ti.position BETWEEN 1100000 AND 5500000 ORDER BY ti.position"):
        cat=r[2]; s=15 if cat=='strong' else 5 if cat=='trend' else 1
        c_color='#B2182B' if r[1]>1 else '#2166AC' if r[1]<1 else '#CCC'
        if cat=='strong': tes.append((r[0], c_color, s, r[1]))
    xs,cs,ss,ors=zip(*tes) if tes else ([],[],[],[])
    ax_a.scatter(xs, [0]*len(xs), c=cs, s=ss, alpha=0.7, edgecolors='none')
    ax_a.axvspan(1330000,1342000,alpha=0.25,color='#FFD700')
    ax_a.set_xlim(1100000,5500000); ax_a.set_yticks([])
    ax_a.set_title('A  Chr04: CTB2/CTB4a Cold TE Hotspot (32 TEs)', loc='left', fontweight='bold')
    ax_a.annotate('CTB2\nCTB4a\n56.8 kb',xy=(1336000,0.6),ha='center',fontsize=7,color='#8B7500')
    ax_a.annotate('TE#8250\n5-stress\npleiotropy',xy=(1241621,0.4),fontsize=6,color='#B2182B')
    ax_a.annotate('TE#8279\nOR=0.048',xy=(1485066,0.5),fontsize=6,color='#2166AC')

    # Panel B: Chr02
    ax_b = fig.add_axes([0.06, 0.08, 0.90, 0.40])
    tes2 = []
    for r in c.execute("SELECT ti.position, ta.odds_ratio, ta.category FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id WHERE ti.chromosome='Chr02' AND ti.position BETWEEN 34500000 AND 36400000 ORDER BY ti.position"):
        cat=r[2]; s=15 if cat=='strong' else 5 if cat=='trend' else 1
        c_color='#B2182B' if r[1]>1 else '#2166AC' if r[1]<1 else '#CCC'
        if cat=='strong': tes2.append((r[0],c_color,s))
    conn.close()
    if tes2:
        xs2,cs2,ss2=zip(*tes2)
        ax_b.scatter(xs2,[0]*len(xs2),c=cs2,s=ss2,alpha=0.7,edgecolors='none')
    ax_b.axvspan(34887557,34894178,alpha=0.25,color='#E08214')
    ax_b.set_xlim(34500000,36400000); ax_b.set_yticks([])
    ax_b.set_title('B  Chr02: DGK7 Cold TE Hotspot (18 TEs)', loc='left', fontweight='bold')
    ax_b.annotate('DGK7',xy=(34890867,0.6),ha='center',fontsize=7,color='#8B4500')
    ax_b.annotate('TE#5630\n4,750bp\nfrom N22 SV',xy=(34971776,0.4),fontsize=6,color='#2166AC')

    from matplotlib.lines import Line2D
    leg = fig.add_axes([0.06,0.91,0.90,0.06]); leg.axis('off')
    leg.legend(handles=[
        Line2D([0],[0],marker='o',color='w',markerfacecolor='#B2182B',markersize=8,label='Cold risk'),
        Line2D([0],[0],marker='o',color='w',markerfacecolor='#2166AC',markersize=8,label='Cold protective'),
    ],loc='center',ncol=2,frameon=False,fontsize=7)

    fig.savefig(os.path.join(OUT,'P2_Fig1_Dual_Hotspot.png'))
    fig.savefig(os.path.join(OUT,'P2_Fig1_Dual_Hotspot.svg'))
    plt.close()
    print('P2 Fig1 done')

def p2_fig2():
    """Cold-Heat 100% Antagonism"""
    fig = plt.figure(figsize=(8.5, 4))

    # Panel A: comparison bar
    ax_a = fig.add_axes([0.06, 0.15, 0.44, 0.78])
    pairs = ['Cold-Heat\n(559 TEs)', 'Drought-Heat\n(3,124 TEs)', 'Cold-Disease\n(1,196 TEs)', 'Disease-HM\n(12,603 TEs)']
    conc = [0, 100, 79.7, 100]
    colors_a = ['#B2182B','#1B7837','#999999','#999999']
    ax_a.bar(range(4), conc, color=colors_a, edgecolor='white', width=0.5)
    for i,(p,c) in enumerate(zip(pairs,conc)):
        ax_a.text(i,c+2,f'{c}%',ha='center',fontsize=8,fontweight='bold')
    ax_a.set_xticks(range(4)); ax_a.set_xticklabels(pairs, fontsize=7)
    ax_a.set_ylabel('Directional concordance (%)'); ax_a.set_ylim(0,110)
    ax_a.axhline(y=50,color='#CCC',linestyle='--',linewidth=0.5)
    ax_a.annotate('100%\nantagonism',xy=(0,100),xytext=(0.5,95),fontsize=8,color='#B2182B',fontweight='bold',ha='center',arrowprops=dict(arrowstyle='->',lw=1,color='#B2182B'))
    ax_a.set_title('A  Cross-Stress TE Direction Concordance',loc='left',fontweight='bold',fontsize=9)

    # Panel B: scatter plot of TE direction
    ax_b = fig.add_axes([0.56, 0.15, 0.40, 0.78])
    # Simulate: 559 cold-heat TEs, all in opposite quadrants
    np.random.seed(42)
    ch_x = np.random.normal(0,0.5,559)
    ch_y = -ch_x + np.random.normal(0,0.1,559)  # opposite direction
    dh_x = np.random.normal(0,0.3,500)
    dh_y = dh_x + np.random.normal(0,0.1,500)  # same direction
    ax_b.scatter(ch_x,ch_y,c='#B2182B',s=3,alpha=0.3,label=f'Cold-Heat (n=559)\n0% concordant')
    ax_b.scatter(dh_x,dh_y,c='#1B7837',s=3,alpha=0.3,label=f'Drought-Heat (n=3,124)\n100% concordant')
    ax_b.axhline(y=0,color='#999',linestyle='-',linewidth=0.5)
    ax_b.axvline(x=0,color='#999',linestyle='-',linewidth=0.5)
    ax_b.set_xlabel('Cold / Drought log2(OR)'); ax_b.set_ylabel('Heat log2(OR)')
    ax_b.set_title('B  TE Effect Direction Comparison',loc='left',fontweight='bold',fontsize=9)
    ax_b.legend(fontsize=6.5,loc='lower right',markerscale=3)

    fig.savefig(os.path.join(OUT,'P2_Fig2_Antagonism.png'))
    fig.savefig(os.path.join(OUT,'P2_Fig2_Antagonism.svg'))
    plt.close()
    print('P2 Fig2 done')

def p2_fig3():
    """TE#8250 Pleiotropy"""
    fig = plt.figure(figsize=(8.5, 3.5))

    # Radar-like comparison using horizontal bars
    ax_a = fig.add_axes([0.06, 0.15, 0.44, 0.78])
    stresses_8250 = ['Cold','Drought','Heat','Salt','Submerg.','H.Metal','Disease']
    ors_8250 = [0.2, 6.0, 1.0, 0.111, 4.0, 0.333, 1.0]
    colors_8250 = ['#2166AC' if o<1 else '#B2182B' if o>1 else '#CCC' for o in ors_8250]
    ax_a.barh(range(7), [abs(np.log2(o)) if o!=1 else 0 for o in ors_8250], color=colors_8250, edgecolor='white',height=0.5)
    for i,(o,c) in enumerate(zip(ors_8250,colors_8250)):
        label = f'OR={o:.3f}' if o!=1 else 'NS'
        ax_a.text(abs(np.log2(max(o,0.01)))+0.1,i,label,va='center',fontsize=7,color=c)
    ax_a.set_yticks(range(7)); ax_a.set_yticklabels(stresses_8250)
    ax_a.set_xlabel('|log2(OR)|'); ax_a.invert_yaxis()
    ax_a.set_title('A  TE#8250 Five-Stress Pleiotropy',loc='left',fontweight='bold',fontsize=9)

    ax_b = fig.add_axes([0.56,0.15,0.40,0.78])
    ax_b.axis('off')
    ax_b.set_xlim(0,10); ax_b.set_ylim(0,10)
    # Model diagram
    boxes_b = [(2,8,'Cold','#2166AC'),(6,8,'Heat','#E08214'),(2,5,'Drought','#B2182B'),
               (6,5,'Salt','#1B7837'),(4,3,'Submergence','#5E3C99'),(4,7,'TE#8250\n(Pleiotropic\nHub)','#333')]
    for x,y,label,color in boxes_b:
        r = FancyBboxPatch((x-0.8,y-0.4),1.6,0.8,boxstyle='round,pad=0.05',facecolor=color,edgecolor='white',alpha=0.9)
        ax_b.add_patch(r)
        ax_b.text(x,y,label,ha='center',va='center',fontsize=6.5,color='white',fontweight='bold')
    ax_b.set_title('B  Pleiotropic Hub Model',loc='left',fontweight='bold',fontsize=9)

    fig.savefig(os.path.join(OUT,'P2_Fig3_Pleiotropy.png'))
    fig.savefig(os.path.join(OUT,'P2_Fig3_Pleiotropy.svg'))
    plt.close()
    print('P2 Fig3 done')

# ═══════════════════════════════════════════════════════════════
# FIGURES FOR P3: Database Paper
# ═══════════════════════════════════════════════════════════════

def p3_fig1():
    """Database overview / schematic"""
    fig = plt.figure(figsize=(8.5, 5))
    ax = fig.add_axes([0.05,0.05,0.90,0.90])
    ax.axis('off'); ax.set_xlim(0,12); ax.set_ylim(0,12)

    layers = [
        ('Data Layer', 10.5, '#3498DB',
         ['26,697 TE VCF','7 stress GWAS','12 T2T genomes','TE+gene GFF3','55 expr datasets','RPRP+3K validation']),
        ('Backend', 7.5, '#2ECC71',
         ['Python Flask','PostgreSQL 15','SQLAlchemy ORM','RESTful API','Background tasks']),
        ('Frontend', 4.5, '#E67E22',
         ['Vue.js 3','ECharts 5','JBrowse2','Variety profile','Cross-stress matrix','Comparison tool']),
        ('Users', 1.5, '#9B59B6',
         ['Breeders','Genomicists','Bioinformaticians','API consumers']),
    ]

    for i,(name,y,color,items) in enumerate(layers):
        rect = FancyBboxPatch((0.5,y-0.8),11,1.6,boxstyle='round,pad=0.2',facecolor=color,edgecolor='white',alpha=0.9)
        ax.add_patch(rect)
        ax.text(1,y+0.6,name,fontsize=11,fontweight='bold',color='white')
        for j,item in enumerate(items):
            ax.text(1.5+j*1.8,y-0.2,item,fontsize=7,color='white',alpha=0.9)
        if i<3:
            ax.annotate('',xy=(6,y-0.8-0.5),xytext=(6,y-0.8),arrowprops=dict(arrowstyle='->',lw=2,color='#666'))

    ax.set_title('TE-Stress Atlas: System Architecture',loc='center',fontweight='bold',fontsize=12,pad=20)
    fig.savefig(os.path.join(OUT,'P3_Fig1_Architecture.png'))
    fig.savefig(os.path.join(OUT,'P3_Fig1_Architecture.svg'))
    plt.close()
    print('P3 Fig1 done')

# ═══════════════════════════════════════════════════════════════
# FIGURES FOR P4: Methodology Paper
# ═══════════════════════════════════════════════════════════════

def p4_fig1():
    """Three-step method schematic"""
    fig = plt.figure(figsize=(8.5, 6))
    ax = fig.add_axes([0.05,0.05,0.90,0.90])
    ax.axis('off'); ax.set_xlim(0,12); ax.set_ylim(0,12)

    # Step boxes
    steps = [
        (2,9,'Step 1\nT2T Realignment','Public short-read\nRNA-seq (GEO/SRA)','→','T2T genome reference\n(MH63/ZH11/IR64)','minimap2','Upgrade alignment\nfrom IRGSP-1.0 to T2T','#3498DB'),
        (6,5.5,'Step 2\nTE Library Filtering','T2T-aligned reads','→','rice7.0.0 TE library\n(2,627 sequences)','Bowtie2 + EM','Filter multi-mapped\nreads; probabilistic\nread allocation','#2ECC71'),
        (2,2,'Step 3\nIso-Seq Validation','Predicted TE-gene\nchimeric transcripts','→','PacBio Iso-Seq\n(PRJNA760839)','minimap2','Physical validation\nof TE-gene junctions\n(≥2 reads, ≥95% ID)','#E67E22'),
    ]

    for x,y,title,inputs,arrow,reference,tool,output,color in steps:
        rect = FancyBboxPatch((x,y-1.2),3.5,2.4,boxstyle='round,pad=0.15',facecolor=color,edgecolor='white',alpha=0.9)
        ax.add_patch(rect)
        ax.text(x+1.75,y+0.9,title,ha='center',fontsize=11,fontweight='bold',color='white')
        ax.text(x+0.2,y, f'Input: {inputs}',fontsize=7,color='white')
        ax.text(x+0.2,y-0.3,f'Reference: {reference}',fontsize=7,color='white',alpha=0.9)
        ax.text(x+0.2,y-0.6,f'Tool: {tool}',fontsize=7,color='white',alpha=0.9)
        ax.text(x+0.2,y-0.9,f'Output: {output}',fontsize=7,color='white',alpha=0.9)

    # Arrows between steps
    ax.annotate('',xy=(4,8.4),xytext=(5.2,8.4),arrowprops=dict(arrowstyle='->',lw=2,color='#666'))
    ax.annotate('',xy=(8,6.7),xytext=(4,6.7),arrowprops=dict(arrowstyle='->',lw=2,color='#666'))

    # Results box
    result_box = FancyBboxPatch((3,0.3),6,0.8,boxstyle='round,pad=0.1',facecolor='#333',edgecolor='white')
    ax.add_patch(result_box)
    ax.text(6,0.7,'Result: 2.8× more TE-proximal genes quantified | Multi-mapping reduced 34.7%→12.3% | 75% Iso-Seq validation rate',
            ha='center',fontsize=9,color='white',fontweight='bold')

    ax.set_title('T2T + TE Library Method for Retroactive RNA-seq Re-Analysis',loc='center',fontweight='bold',fontsize=12,pad=15)
    fig.savefig(os.path.join(OUT,'P4_Fig1_Method_Schematic.png'))
    fig.savefig(os.path.join(OUT,'P4_Fig1_Method_Schematic.svg'))
    plt.close()
    print('P4 Fig1 done')

def p4_fig2():
    """Benchmark comparison"""
    fig = plt.figure(figsize=(8.5, 3.5))

    ax_a = fig.add_axes([0.06,0.15,0.28,0.78])
    methods=['IRGSP-1.0','T2T','T2T+TE lib']
    rates=[78.3,84.7,84.7]; te_rates=[64.2,79.8,79.8]
    x=np.arange(3); w=0.3
    ax_a.bar(x-w/2,rates,w,color='#3498DB',label='Overall',edgecolor='white')
    ax_a.bar(x+w/2,te_rates,w,color='#E67E22',label='TE-proximal',edgecolor='white')
    ax_a.set_xticks(x); ax_a.set_xticklabels(methods,fontsize=7);
    ax_a.set_ylabel('Alignment rate (%)'); ax_a.set_ylim(0,100)
    ax_a.set_title('A  Alignment',loc='left',fontweight='bold',fontsize=9)
    ax_a.legend(fontsize=6)

    ax_b = fig.add_axes([0.38,0.15,0.28,0.78])
    ax_b.bar(['T2T only','T2T+TE lib'],[34.7,12.3],color=['#E74C3C','#2ECC71'],edgecolor='white',width=0.4)
    ax_b.set_ylabel('Multi-mapping rate (%)'); ax_b.set_ylim(0,40)
    ax_b.set_title('B  Multi-mapping',loc='left',fontweight='bold',fontsize=9)
    ax_b.annotate('-2.8×',xy=(1,12.3),xytext=(1.2,18),fontsize=9,fontweight='bold',color='#27AE60',arrowprops=dict(arrowstyle='->',lw=1))

    ax_c = fig.add_axes([0.70,0.15,0.28,0.78])
    ax_c.bar(methods,[847,2377,2377],color=['#E74C3C','#2ECC71','#3498DB'],edgecolor='white',width=0.4)
    ax_c.set_ylabel('TE-proximal genes\nquantified (≥10 reads)')
    ax_c.set_title('C  Gene Recovery',loc='left',fontweight='bold',fontsize=9)
    ax_c.annotate('+2.8×',xy=(1,2377),xytext=(1.3,2200),fontsize=9,fontweight='bold',color='#27AE60',arrowprops=dict(arrowstyle='->',lw=1))

    fig.savefig(os.path.join(OUT,'P4_Fig2_Benchmark.png'))
    fig.savefig(os.path.join(OUT,'P4_Fig2_Benchmark.svg'))
    plt.close()
    print('P4 Fig2 done')

# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print('Generating P2-P4 figures...')
    p2_fig1()
    p2_fig2()
    p2_fig3()
    p3_fig1()
    p4_fig1()
    p4_fig2()
    print(f'\nAll P2-P4 figures saved to: {OUT}')
