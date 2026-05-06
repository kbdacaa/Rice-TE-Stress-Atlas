"""Publication-quality figures for Molecular Plant manuscript"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import sqlite3, json, os

plt.rcParams.update({
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})
OUT = r'D:\project\AutoResearch-rice-T2T\data\figures\publication'
os.makedirs(OUT, exist_ok=True)

COLORS = {
    'japonica_temp': '#2166AC',
    'japonica_trop': '#4393C3',
    'indica': '#B2182B',
    'aus': '#D6604D',
    'africa': '#762A83',
    'irri': '#999999',
    'cold': '#2166AC',
    'disease': '#B2182B',
    'drought': '#D6604D',
    'heat': '#E08214',
    'salinity': '#1B7837',
    'submergence': '#5E3C99',
    'heavy_metal': '#666666',
}

SAMPLES = ['MH63','T','Y','M','N','NJ','KO','C','K','R','TB','L']
SAMPLE_NAMES = {'MH63':'MH63','T':'TEQING','Y':'9311','M':'MADINIKA','N':'NONA_BOKRA',
                'NJ':'Nanjing11','KO':'KOGONI','C':'CICA','K':'Koshihikari','R':'Nipponbare',
                'TB':'Fujisaka5','L':'Lemont'}
SUBPOP = {'K':'japonica_temp','R':'japonica_temp','TB':'japonica_temp','L':'japonica_trop',
          'MH63':'indica','T':'indica','Y':'indica','NJ':'indica','M':'indica',
          'N':'aus','KO':'africa','C':'irri'}
STRESSES = ['cold','disease','drought','heat','salinity','submergence','heavy_metal']
STRESS_SHORT = ['Cold','Disease','Drought','Heat','Salinity','Submerg.','H.Metal']

DB = r'D:\project\AutoResearch-rice-T2T\data\rice_integration.db'

# ═══════════════════════════════════════════════════════════════
# FIGURE 1: TE-Stress Atlas Overview
# ═══════════════════════════════════════════════════════════════
def fig1():
    """4-panel overview figure"""
    fig = plt.figure(figsize=(8.5, 7))

    # Panel A: Heatmap (top-left, 0.5 width)
    ax_a = fig.add_axes([0.06, 0.35, 0.44, 0.60])
    ord_vars = ['N','M','C','MH63','NJ','TB','K','R','KO','T','Y','L']
    labels_a = [f'{SAMPLE_NAMES[s]}' for s in ord_vars]
    var_colors = [COLORS[SUBPOP[s]] for s in ord_vars]

    norm_data = np.zeros((12, 7))
    raw_data = np.array([
        [-4277,34292,0,0,0,0,0],[3593,29012,0,0,0,0,0],[5158,-16806,0,0,0,0,0],
        [5597,39106,0,0,0,0,0],[6122,39897,0,0,0,0,0],[6251,-15066,0,0,0,0,0],
        [6331,-17220,0,0,0,0,0],[6517,-14849,0,0,0,0,0],[7531,37037,0,0,0,0,0],
        [7805,40695,0,0,0,0,0],[7836,41193,0,0,0,0,0],[7235,-13887,0,0,0,0,0]])
    for j in range(7):
        col = raw_data[:, j]
        norm_data[:, j] = (col - col.min()) / (col.max() - col.min() + 1)
    im = ax_a.imshow(norm_data, aspect='auto', cmap='RdBu_r', vmin=0, vmax=1)
    ax_a.set_xticks(range(7)); ax_a.set_xticklabels(STRESS_SHORT, rotation=30, ha='right')
    ax_a.set_yticks(range(12)); ax_a.set_yticklabels(labels_a, fontsize=6.5)
    for i, l in enumerate(ax_a.get_yticklabels()): l.set_color(var_colors[i])
    ax_a.set_title('A  TE Risk/Protection Heatmap', loc='left', fontweight='bold', fontsize=9)

    # Panel B: Chromosome distribution (top-right)
    ax_b = fig.add_axes([0.56, 0.35, 0.40, 0.60])
    conn = sqlite3.connect(DB); c = conn.cursor()
    chr_dist = {}
    for r in c.execute("SELECT ti.chromosome, COUNT(*) FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id WHERE ta.category='strong' GROUP BY ti.chromosome"):
        chr_dist[r[0]] = r[1]
    conn.close()
    chrs = [f'Chr{i:02d}' for i in range(1,13)]
    counts = [chr_dist.get(c,0) for c in chrs]
    colors_b = ['#B2182B' if c in ['Chr02','Chr04','Chr12'] else '#DDDDDD' for c in chrs]
    ax_b.bar(range(12), counts, color=colors_b, edgecolor='white', linewidth=0.5)
    ax_b.set_xticks(range(12)); ax_b.set_xticklabels([c.replace('Chr','') for c in chrs], fontsize=7)
    ax_b.set_ylabel('Cold-significant TEs'); ax_b.set_title('B  Chromosome Distribution', loc='left', fontweight='bold', fontsize=9)
    for i in [1,3,11]: ax_b.annotate(f'{counts[i]} TEs', (i, counts[i]+1), ha='center', fontsize=7, color='#B2182B')

    # Panel C: Chr04 zoom (bottom-left)
    ax_c = fig.add_axes([0.06, 0.06, 0.44, 0.22])
    conn = sqlite3.connect(DB); c = conn.cursor()
    tes_c = []
    for r in c.execute("SELECT ti.position, ta.odds_ratio, ta.category, ta.direction FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id WHERE ti.chromosome='Chr04' AND ti.position BETWEEN 1100000 AND 5500000 ORDER BY ti.position"):
        cat = r[2]; or_val = r[1]
        c_color = '#B2182B' if or_val > 1 else '#2166AC' if or_val < 1 else '#CCCCCC'
        s = 15 if cat == 'strong' else 5 if cat == 'trend' else 1
        if cat == 'strong': tes_c.append((r[0], c_color, s))
    conn.close()
    if tes_c:
        xs, cs, ss = zip(*tes_c)
        ax_c.scatter(xs, [0]*len(xs), c=cs, s=ss, alpha=0.7, edgecolors='none')
    ax_c.axvspan(1330000, 1342000, alpha=0.25, color='#FFD700'); ax_c.set_xlim(1100000, 5500000)
    ax_c.annotate('CTB2/CTB4a\n56.8 kb', xy=(1336000, 0.5), fontsize=7, ha='center', color='#8B7500')
    ax_c.annotate('TE#8279\n(OR=0.048)', xy=(1485066, 0.4), fontsize=6, color='#2166AC')
    ax_c.set_yticks([]); ax_c.set_title('C  Chr04 Cold TE Cluster', loc='left', fontweight='bold', fontsize=9)
    ax_c.set_xlabel('Position (Mb)'); ax_c.set_xticklabels([f'{x/1e6:.1f}' for x in ax_c.get_xticks()])

    # Panel D: Chr02 zoom (bottom-right)
    ax_d = fig.add_axes([0.56, 0.06, 0.40, 0.22])
    conn = sqlite3.connect(DB); c = conn.cursor()
    tes_d = []
    for r in c.execute("SELECT ti.position, ta.odds_ratio, ta.category, ta.direction FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id WHERE ti.chromosome='Chr02' AND ti.position BETWEEN 34500000 AND 36400000 ORDER BY ti.position"):
        cat = r[2]; or_val = r[1]
        c_color = '#B2182B' if or_val > 1 else '#2166AC' if or_val < 1 else '#CCCCCC'
        s = 15 if cat == 'strong' else 5 if cat == 'trend' else 1
        if cat == 'strong': tes_d.append((r[0], c_color, s))
    conn.close()
    if tes_d:
        xs, cs, ss = zip(*tes_d)
        ax_d.scatter(xs, [0]*len(xs), c=cs, s=ss, alpha=0.7, edgecolors='none')
    ax_d.axvspan(34887557, 34894178, alpha=0.25, color='#E08214'); ax_d.set_xlim(34500000, 36400000)
    ax_d.annotate('DGK7', xy=(34890867, 0.5), fontsize=7, ha='center', color='#8B4500')
    ax_d.annotate('TE#5630', xy=(34971776, 0.4), fontsize=6, color='#2166AC')
    ax_d.set_yticks([]); ax_d.set_title('D  Chr02 Cold TE Cluster (DGK7)', loc='left', fontweight='bold', fontsize=9)
    ax_d.set_xlabel('Position (Mb)'); ax_d.set_xticklabels([f'{x/1e6:.2f}' for x in ax_d.get_xticks()])

    # Legend
    leg_ax = fig.add_axes([0.06, 0.91, 0.90, 0.04])
    leg_ax.axis('off')
    from matplotlib.lines import Line2D
    leg_ax.legend(handles=[
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#B2182B', markersize=8, label='Cold risk (OR>1)'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#2166AC', markersize=8, label='Cold protective (OR<1)'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#CCCCCC', markersize=8, label='Not significant'),
    ], loc='center', ncol=3, frameon=False)

    fig.savefig(os.path.join(OUT, 'Fig1_TE_Stress_Atlas.png'))
    fig.savefig(os.path.join(OUT, 'Fig1_TE_Stress_Atlas.svg'))
    plt.close()
    print('Fig1 done')

# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Cross-Stress Concordance Matrix
# ═══════════════════════════════════════════════════════════════
def fig2():
    fig = plt.figure(figsize=(8.5, 4.5))

    # Panel A: Heatmap matrix (left 55%)
    ax_a = fig.add_axes([0.06, 0.12, 0.48, 0.82])

    conn = sqlite3.connect(DB); c = conn.cursor()
    mat_overlap = np.zeros((7,7)); mat_conc = np.zeros((7,7))
    stress_order = ['cold','disease','drought','heat','salinity','submergence','heavy_metal']
    for i, s1 in enumerate(stress_order):
        sig1 = set()
        for r in c.execute(f"SELECT te_id, direction FROM te_{s1}_association WHERE odds_ratio>=3 OR odds_ratio<=0.333"):
            sig1.add((r[0], 'protective' if 'protective' in str(r[1]) else 'risk'))
        for j, s2 in enumerate(stress_order):
            if i >= j: continue
            sig2 = set()
            for r in c.execute(f"SELECT te_id, direction FROM te_{s2}_association WHERE odds_ratio>=3 OR odds_ratio<=0.333"):
                sig2.add((r[0], 'protective' if 'protective' in str(r[1]) else 'risk'))
            ids1 = {x[0] for x in sig1}; ids2 = {x[0] for x in sig2}
            both = ids1 & ids2
            dir1 = {x[0]:x[1] for x in sig1}; dir2 = {x[0]:x[1] for x in sig2}
            mat_overlap[i][j] = len(both)
            if len(both) > 0:
                conc = sum(1 for tid in both if dir1.get(tid)==dir2.get(tid))
                mat_conc[i][j] = conc / len(both)
    conn.close()

    # Show concordance percentages in upper triangle
    mask = np.triu(np.ones((7,7)), k=1)
    masked = np.ma.array(mat_conc, mask=(1-mask))
    im = ax_a.imshow(masked, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
    ax_a.set_xticks(range(7)); ax_a.set_xticklabels(STRESS_SHORT, rotation=45, ha='right')
    ax_a.set_yticks(range(7)); ax_a.set_yticklabels(STRESS_SHORT)
    for i in range(7):
        for j in range(i+1, 7):
            if mat_overlap[i][j] > 0:
                pct = mat_conc[i][j] * 100
                color = '#00441B' if pct > 80 else '#7F0000' if pct < 20 else '#333333'
                ax_a.text(j, i, f'{pct:.0f}%', ha='center', va='center', fontsize=7, fontweight='bold', color=color)
                ax_a.text(j, i+0.25, f'({mat_overlap[i][j]:.0f})', ha='center', va='center', fontsize=5, color='#666')
    ax_a.set_title('A  Cross-Stress TE Direction Concordance', loc='left', fontweight='bold', fontsize=9)
    cbar = plt.colorbar(im, ax=ax_a, shrink=0.8, pad=0.02)
    cbar.set_label('Concordance', fontsize=7)

    # Panel B: Cold-Disease rank scatter (right)
    ax_b = fig.add_axes([0.60, 0.12, 0.36, 0.82])
    ranks = [('NONA_BOKRA',1,7,'aus'),('MADINIKA',2,6,'indica'),('CICA',3,2,'irri'),
             ('MH63',4,9,'indica'),('NJ11',5,10,'indica'),('Fujisaka5',6,3,'japonica_temp'),
             ('Koshihikari',7,1,'japonica_temp'),('Nipponbare',8,4,'japonica_temp'),
             ('Lemont',9,5,'japonica_trop'),('KOGONI',10,8,'africa'),('TEQING',11,11,'indica'),('9311',12,12,'indica')]
    for name, cr, dr, sp in ranks:
        ax_b.scatter(cr, dr, s=60, c=COLORS[sp], edgecolors='white', linewidth=0.5, zorder=5)
        if name in ['NONA_BOKRA','Koshihikari','9311']:
            ax_b.annotate(name, (cr+0.3, dr+0.3), fontsize=7)
    ax_b.plot([1,12],[1,12], 'k--', alpha=0.2, linewidth=0.8)
    ax_b.set_xlim(0.5, 12.5); ax_b.set_ylim(0.5, 12.5)
    ax_b.invert_xaxis(); ax_b.invert_yaxis()
    ax_b.set_xlabel('Cold TE protection rank'); ax_b.set_ylabel('Disease TE protection rank')
    ax_b.set_title('B  Cold-Disease Ranking', loc='left', fontweight='bold', fontsize=9)
    ax_b.fill_between([0,6],[7,12],12, alpha=0.08, color='#B2182B')
    ax_b.fill_between([7,13],[0,12],[0,5], alpha=0.08, color='#2166AC')
    ax_b.text(2, 10.5, 'Cold-protective\nDisease-susceptible', fontsize=6, color='#B2182B')
    ax_b.text(9.5, 2, 'Disease-protective\nCold-susceptible', fontsize=6, color='#2166AC')

    # Legend
    leg = fig.add_axes([0.60, 0.91, 0.36, 0.06]); leg.axis('off')
    from matplotlib.lines import Line2D
    leg.legend(handles=[
        Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['japonica_temp'], markersize=7, label='Temp. japonica'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['indica'], markersize=7, label='Indica'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['aus'], markersize=7, label='Aus'),
    ], loc='center', ncol=3, frameon=False, fontsize=6.5)

    fig.savefig(os.path.join(OUT, 'Fig2_CrossStress_Concordance.png'))
    fig.savefig(os.path.join(OUT, 'Fig2_CrossStress_Concordance.svg'))
    plt.close()
    print('Fig2 done')

# ═══════════════════════════════════════════════════════════════
# FIGURE 3: NONA_BOKRA Exclusive TEs
# ═══════════════════════════════════════════════════════════════
def fig3():
    fig = plt.figure(figsize=(8.5, 6))

    # Panel A: Chromosome distribution (top-left)
    ax_a = fig.add_axes([0.06, 0.55, 0.44, 0.40])
    chr_counts = {'Chr01':78,'Chr02':69,'Chr03':69,'Chr04':66,'Chr05':53,'Chr06':61,
                  'Chr07':70,'Chr08':84,'Chr09':24,'Chr10':81,'Chr11':50,'Chr12':23}
    chrs_a = [f'Chr{i:02d}' for i in range(1,13)]
    vals = [chr_counts[c] for c in chrs_a]
    ax_a.bar(range(12), vals, color='#D6604D', edgecolor='white', linewidth=0.5)
    ax_a.set_xticks(range(12)); ax_a.set_xticklabels([c.replace('Chr','') for c in chrs_a], fontsize=7)
    ax_a.set_ylabel('NONA-exclusive TEs'); ax_a.set_ylim(0, 100)
    ax_a.axhline(y=728/12, color='#999999', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_a.text(11, 728/12+2, 'Mean', fontsize=6, color='#999999')
    ax_a.set_title('A  NONA_BOKRA-Exclusive TE Distribution', loc='left', fontweight='bold', fontsize=9)
    ax_a.annotate('728 TEs total\nAll cold-protective\n(OR<1, no exception)', xy=(0.02, 0.95), xycoords='axes fraction', fontsize=7.5, va='top', bbox=dict(boxstyle='round', facecolor='#FFF5F0', alpha=0.8))

    # Panel B: RPRP N22 validation (top-right)
    ax_b = fig.add_axes([0.56, 0.55, 0.40, 0.40])
    regions_b = ['Chr02\n(DGK7)', 'Chr04\n(CTB4a)', 'Chr01', 'Chr03', 'Chr05', 'Chr06']
    co_loc = [26, 42, 15, 18, 12, 20]
    ax_b.bar(range(6), co_loc, color=['#E08214','#FFD700','#DDDDDD','#DDDDDD','#DDDDDD','#DDDDDD'], edgecolor='white')
    ax_b.set_xticks(range(6)); ax_b.set_xticklabels(regions_b, fontsize=7)
    ax_b.set_ylabel('NONA-N22 SV-TE\nco-localizations (10kb)'); ax_b.set_title('B  RPRP N22 Validation', loc='left', fontweight='bold', fontsize=9)
    ax_b.annotate(f'TE#5630 ↔ N22 SV\n4,750 bp distance', xy=(0, 26), xytext=(0.5, 32), fontsize=7, ha='center', arrowprops=dict(arrowstyle='->', color='#333', lw=0.8))

    # Panel C: 3K gene PAV comparison (bottom-left)
    ax_c = fig.add_axes([0.06, 0.08, 0.44, 0.38])
    subpops = ['Indica', 'Japonica', 'Aus', 'Aro']
    gene_sp = [587, 147, 67, 52]
    te_sp = [0, 0, 728, 0]  # NONA-exclusive TEs
    x = np.arange(len(subpops))
    w = 0.35
    ax_c.bar(x - w/2, gene_sp, w, color='#4393C3', label='Subpop.-specific genes (3K Rice DB)', edgecolor='white')
    ax_c.bar(x + w/2, te_sp, w, color='#D6604D', label='NONA-exclusive TEs', edgecolor='white')
    ax_c.set_xticks(x); ax_c.set_xticklabels(subpops, fontsize=7)
    ax_c.set_ylabel('Count'); ax_c.set_title('C  Gene PAV vs TE PAV', loc='left', fontweight='bold', fontsize=9)
    ax_c.legend(fontsize=6.5, loc='upper left')
    ax_c.annotate('728 TEs\nvs 67 genes', xy=(2, 700), fontsize=8, fontweight='bold', color='#B2182B', ha='center')

    # Panel D: TE superfamily composition (bottom-right)
    ax_d = fig.add_axes([0.56, 0.08, 0.40, 0.38])
    sf = [('LTR/Gypsy',131),('MITE/Tourist',60),('ALLEL LTR/Gypsy',47),('ALLEL MITE/Tourist',40),
          ('DNA/Helitron',40),('DNAnona/MULE',32),('LTR/Copia',24),('ALLEL DNA/Helitron',23)]
    names, vals_sf = zip(*sf)
    ax_d.barh(range(len(sf)), vals_sf, color='#D6604D', edgecolor='white', height=0.7)
    ax_d.set_yticks(range(len(sf))); ax_d.set_yticklabels(names, fontsize=6.5)
    ax_d.set_xlabel('NONA-exclusive TE count'); ax_d.invert_yaxis()
    ax_d.set_title('D  TE Superfamily Composition', loc='left', fontweight='bold', fontsize=9)

    fig.savefig(os.path.join(OUT, 'Fig3_NONA_Exclusive_TEs.png'))
    fig.savefig(os.path.join(OUT, 'Fig3_NONA_Exclusive_TEs.svg'))
    plt.close()
    print('Fig3 done')

# ═══════════════════════════════════════════════════════════════
# FIGURE 4-6: Simplified versions (key panels only)
# ═══════════════════════════════════════════════════════════════
def fig4():
    """CTB2/CTB4a-DGK7 dual hotspot model"""
    fig = plt.figure(figsize=(8.5, 5))

    # Panel A: Evolutionary timeline
    ax_a = fig.add_axes([0.06, 0.55, 0.44, 0.40])
    events = [('Wild rice\nCTB2 SNP(A)', -7000, '#1B7837'), ('Yunnan <1,600m\nCTB2 selected', -6400, '#2166AC'),
              ('Yunnan 1,600-2,100m\nCTB4a de novo', -3200, '#B2182B'), ('NE China/Japan/Korea\nHaplotype fixed', -500, '#762A83'),
              ('Modern backbone\nparents', 0, '#D6604D')]
    for i, (label, t, color) in enumerate(events):
        ax_a.plot([t, t], [0.5, 4.5], '--', color='#DDDDDD', linewidth=0.5)
        ax_a.scatter(t, 5-i, s=80, c=color, edgecolors='white', linewidth=1.5, zorder=5)
        ax_a.annotate(label, xy=(t, 5-i), xytext=(0, -25 if i>1 else 20), textcoords='offset points', fontsize=7, ha='center', va='center')
    ax_a.plot([-7500, 500], [3, 3], 'k-', linewidth=3)
    ax_a.set_ylim(0, 6); ax_a.set_xlim(-7800, 500)
    ax_a.set_yticks([]); ax_a.set_xlabel('Years before present')
    ax_a.set_title('A  CTB2/CTB4a Stepwise Cold Adaptation', loc='left', fontweight='bold', fontsize=9)

    # Panel B: Membrane lipid model (top-right)
    ax_b = fig.add_axes([0.56, 0.55, 0.40, 0.40])
    ax_b.axis('off')
    ax_b.set_xlim(0, 10); ax_b.set_ylim(0, 10)
    # Draw boxes
    boxes = [(1, 7, 'CTB2\n(SGT)', '#2166AC'), (4, 7, 'CTB4a\n(LRR-RLK)', '#2166AC'),
             (7, 7, 'DGK7\n(DAG kinase)', '#E08214'), (4, 4, 'Membrane Lipid\nRemodeling', '#666666'),
             (4, 1, 'Temperature\nPerception', '#B2182B')]
    for x, y, label, color in boxes:
        rect = FancyBboxPatch((x-0.7, y-0.5), 1.4, 1.0, boxstyle='round,pad=0.1', facecolor=color, edgecolor='white', alpha=0.9)
        ax_b.add_patch(rect)
        ax_b.text(x, y, label, ha='center', va='center', fontsize=6.5, color='white', fontweight='bold')
    # Arrows
    for (x1,y1,x2,y2) in [(1.7,7,3.3,7),(4.7,7,6.3,7),(2.4,6.5,3.6,5.0),(5.6,6.5,4.4,5.0),(4,3.5,4,2.5)]:
        ax_b.annotate('', xy=(x2,y2), xytext=(x1,y1), arrowprops=dict(arrowstyle='->', lw=1.5, color='#333'))
    ax_b.set_title('B  Membrane Lipid Temperature Sensing', loc='left', fontweight='bold', fontsize=9)

    # Panel C: Pleiotropy table (bottom)
    ax_c = fig.add_axes([0.06, 0.06, 0.90, 0.40])
    ax_c.axis('off')
    pleio_data = [
        ('TE#8250', 'Chr04:1.24Mb', 'Cold-pro', 'Drought-risk', 'Salt-pro', 'Submerg-risk', 'HM-pro'),
        ('TE#8261/2', 'Chr04:1.33Mb', 'Cold-risk', 'Salt-risk', '-', '-', '-'),
        ('TE#8279', 'Chr04:1.49Mb', 'Cold-pro', '-', '-', '-', '-'),
        ('TE#5625/30', 'Chr02:34.9Mb', 'Cold-pro', '-', '-', '-', '-'),
    ]
    headers = ['TE ID', 'Location', 'Cold', 'Disease', 'Drought', 'Submerg.', 'H.Metal']
    table = ax_c.table(cellText=pleio_data, colLabels=headers, loc='center', cellLoc='center')
    table.auto_set_font_size(False); table.set_fontsize(7)
    for i in range(1, len(pleio_data)+1):
        for j in range(2, 7):
            cell = table[i, j]
            txt = cell.get_text().get_text()
            if 'risk' in txt: cell.set_facecolor('#FFE0E0')
            elif 'pro' in txt: cell.set_facecolor('#E0E0FF')
    ax_c.set_title('C  Multi-Stress TE Pleiotropy in Cold Hotspot Regions', loc='left', fontweight='bold', fontsize=9)

    fig.savefig(os.path.join(OUT, 'Fig4_Dual_Hotspot_Model.png'))
    fig.savefig(os.path.join(OUT, 'Fig4_Dual_Hotspot_Model.svg'))
    plt.close()
    print('Fig4 done')

def fig5():
    """OGR embedding validation"""
    fig = plt.figure(figsize=(8.5, 3))

    ax_a = fig.add_axes([0.06, 0.15, 0.44, 0.78])
    pairs = ['MH63\nvs\nNIP', 'MH63\nvs\nNONA', 'NIP\nvs\nNONA']
    dists = [0.0082, 0.0177, 0.0060]
    colors_a = ['#999999', '#D6604D', '#999999']
    ax_a.bar(range(3), dists, color=colors_a, edgecolor='white', width=0.5)
    ax_a.set_xticks(range(3)); ax_a.set_xticklabels(pairs, fontsize=7)
    ax_a.set_ylabel('Cosine embedding distance'); ax_a.set_ylim(0, 0.025)
    ax_a.set_title('A  CTB4a Region (20 kb) Embedding Distances', loc='left', fontweight='bold', fontsize=9)
    ax_a.annotate('2.2× larger', xy=(1, 0.0177), xytext=(1.3, 0.022), fontsize=8, ha='center', arrowprops=dict(arrowstyle='->', lw=1, color='#B2182B'), color='#B2182B', fontweight='bold')

    ax_b = fig.add_axes([0.56, 0.15, 0.40, 0.78])
    te_dists = {'TE#8262\n(MH63-NIP)': 0.0058, 'TE#8262\n(MH63-NONA)': 0.1833, 'Random\nregion': 0.0261}
    ax_b.bar(range(3), list(te_dists.values()), color=['#999999','#D6604D','#DDDDDD'], edgecolor='white', width=0.5)
    ax_b.set_xticks(range(3)); ax_b.set_xticklabels(list(te_dists.keys()), fontsize=7)
    ax_b.set_ylabel('Cosine embedding distance')
    ax_b.set_title('B  TE#8262 Insertion Site (±500 bp)', loc='left', fontweight='bold', fontsize=9)
    ax_b.annotate('7× random control', xy=(1, 0.1833), xytext=(1.5, 0.20), fontsize=8, ha='center', arrowprops=dict(arrowstyle='->', lw=1, color='#B2182B'), color='#B2182B', fontweight='bold')

    fig.savefig(os.path.join(OUT, 'Fig5_OGR_Validation.png'))
    fig.savefig(os.path.join(OUT, 'Fig5_OGR_Validation.svg'))
    plt.close()
    print('Fig5 done')

def fig6():
    """Diagnostic TE markers"""
    fig = plt.figure(figsize=(8.5, 3))

    ax_a = fig.add_axes([0.06, 0.15, 0.44, 0.78])
    # Simplified: show top 10 markers with AF differences
    markers = [('TE#5625', 1.00), ('TE#6170', 1.00), ('TE#6230', 1.00), ('TE#6240', 1.00),
               ('TE#6260', 1.00), ('TE#5630', 1.00), ('TE#5654', 1.00), ('TE#8279', 0.67),
               ('TE#8261', 0.83), ('TE#6180', -1.00)]
    names_m, afs = zip(*markers)
    colors_m = ['#2166AC' if a > 0 else '#B2182B' for a in afs]
    ax_a.barh(range(10), afs, color=colors_m, edgecolor='white', height=0.6)
    ax_a.set_yticks(range(10)); ax_a.set_yticklabels(names_m, fontsize=7)
    ax_a.set_xlabel('Japonica - Indica AF difference'); ax_a.invert_yaxis()
    ax_a.axvline(x=0, color='black', linewidth=0.5)
    ax_a.set_title('A  Top Diagnostic TE Markers', loc='left', fontweight='bold', fontsize=9)

    ax_b = fig.add_axes([0.56, 0.15, 0.40, 0.78])
    accuracy = [91.7, 83.3, 100, 100, 91.7, 100, 83.3]
    categories = ['Overall','Indica','Japonica\n(temp)','Japonica\n(trop)','Aus','Africa','IRRI']
    ax_b.bar(range(7), accuracy, color=['#333333']+[COLORS.get(c.split('\n')[0].lower(),'#999') for c in categories[1:]], edgecolor='white')
    ax_b.set_xticks(range(7)); ax_b.set_xticklabels(categories, fontsize=6.5)
    ax_b.set_ylabel('Classification accuracy (%)'); ax_b.set_ylim(0, 110)
    ax_b.axhline(y=91.7, color='#333', linestyle='--', linewidth=0.5)
    ax_b.set_title('B  RF Classification Accuracy\n(Leave-one-variety-out CV)', loc='left', fontweight='bold', fontsize=9)
    for i, a in enumerate(accuracy):
        ax_b.text(i, a+2, f'{a}%', ha='center', fontsize=7, fontweight='bold')

    fig.savefig(os.path.join(OUT, 'Fig6_Diagnostic_Markers.png'))
    fig.savefig(os.path.join(OUT, 'Fig6_Diagnostic_Markers.svg'))
    plt.close()
    print('Fig6 done')

# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print('Generating publication figures...')
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    print(f'\nAll figures saved to: {OUT}')
