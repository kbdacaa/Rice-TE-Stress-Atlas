"""OGR滑动窗口嵌入分析可视化"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json, sqlite3, os

plt.rcParams.update({
    'font.sans-serif': ['Arial','Helvetica'],
    'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'figure.dpi':600,'savefig.dpi':600,'savefig.bbox':'tight'
})

OUT = r'D:\project\AutoResearch-rice-T2T\data\figures\publication'

# Load data
with open('D:/project/AutoResearch-rice-T2T/data/ogr_sliding_window.json','r') as f:
    data = json.load(f)

sw = data['sliding_window']
fine = data['fine_grained']

# Get TE positions
DB = 'D:/project/AutoResearch-rice-T2T/data/rice_integration.db'
conn = sqlite3.connect(DB); c = conn.cursor()
te_positions = []
for r in c.execute("""
    SELECT ti.position, ta.odds_ratio, ta.category, ta.direction
    FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id
    WHERE ti.chromosome='Chr04' AND ti.position BETWEEN 1100000 AND 5500000
    ORDER BY ti.position
"""):
    cat = r[2]; or_val = r[1]
    if cat == 'strong':
        c_color = '#B2182B' if or_val > 1 else '#2166AC'
        te_positions.append({'pos': r[0], 'color': c_color, 'size': 30, 'or': or_val})
    elif cat == 'trend':
        c_color = '#E8A0A0' if or_val > 1 else '#A0C8E8'
        te_positions.append({'pos': r[0], 'color': c_color, 'size': 10, 'or': or_val})
conn.close()

# ═══ Main Figure ═══
fig = plt.figure(figsize=(10, 8))

# Panel A: Full sliding window overview
ax_a = fig.add_axes([0.08, 0.55, 0.88, 0.40])

# Plot embedding distances as lines
xs = [w['pos']/1e6 for w in sw]
mn = [w['d_MH63_NIP'] for w in sw]
ma = [w['d_MH63_NONA'] for w in sw]
na = [w['d_NIP_NONA'] for w in sw]

ax_a.plot(xs, mn, 'o-', color='#888888', markersize=6, linewidth=1.5, label='MH63(Ind) vs NIP(Jap)', zorder=2)
ax_a.plot(xs, ma, 's-', color='#D6604D', markersize=6, linewidth=1.5, label='MH63(Ind) vs NONA(Aus)', zorder=2)
ax_a.plot(xs, na, '^-', color='#2166AC', markersize=6, linewidth=1.5, label='NIP(Jap) vs NONA(Aus)', zorder=2)

# TE positions as scatter on top
for te in te_positions:
    ax_a.scatter(te['pos']/1e6, 0.21, s=te['size'], c=te['color'], alpha=0.6, edgecolors='none', zorder=1)

# CTB4a region highlight
ax_a.axvspan(1.33, 1.34, alpha=0.2, color='#FFD700', zorder=0)
ax_a.annotate('CTB4a gene\n(conserved)', xy=(1.335, 0.15), fontsize=8, ha='center', color='#8B7500',
              bbox=dict(boxstyle='round', facecolor='#FFD700', alpha=0.2))

# Top peak annotations
ax_a.annotate('★ Peak 1\n0.1946', xy=(2.1, 0.1946), xytext=(2.3, 0.22), fontsize=8, ha='center',
              arrowprops=dict(arrowstyle='->', color='#B2182B', lw=1.5), color='#B2182B', fontweight='bold')
ax_a.annotate('Peak 2\n0.1626', xy=(4.9, 0.1626), xytext=(4.7, 0.19), fontsize=7, ha='center',
              arrowprops=dict(arrowstyle='->', color='#666', lw=1), color='#666')
ax_a.annotate('Peak 3\n0.1122', xy=(4.3, 0.1122), xytext=(4.1, 0.14), fontsize=7, ha='center',
              arrowprops=dict(arrowstyle='->', color='#666', lw=1), color='#666')

ax_a.set_xlim(1.0, 5.6); ax_a.set_ylim(0, 0.25)
ax_a.set_ylabel('Cosine embedding distance'); ax_a.set_xlabel('Chr04 position (Mb)')
ax_a.set_title('A  OGR Embedding Divergence Across Chr04 Cold TE Hotspot (20 kb windows)', loc='left', fontweight='bold')
ax_a.legend(fontsize=7, loc='upper right', ncol=1)

# Gene annotation markers
ax_a.annotate('CTB2', xy=(0.18, 0.22), fontsize=7, color='#FFD700', ha='center')

# Panel B: Zoomed TE cluster divergence (bar chart comparison)
ax_b = fig.add_axes([0.08, 0.06, 0.42, 0.42])

regions = ['CTB4a\ngene body\n(1.33Mb)', 'CTB4a\ndownstream\n(1.5Mb)', 'Peak 1\nTE cluster\n(2.1Mb)',
           'Peak 2\n(4.9Mb)', 'Peak 3\n(4.3Mb)', 'Random\ncontrol\n(3.1Mb)']
values_mn = [0.0076, 0.0996, 0.1946, 0.1626, 0.1122, 0.0366]
values_ma = [0.0071, 0.0106, 0.1385, 0.1165, 0.1159, 0.0711]

x = np.arange(len(regions))
w = 0.35
bars1 = ax_b.bar(x - w/2, values_mn, w, color='#888888', edgecolor='white', label='MH63 vs NIP')
bars2 = ax_b.bar(x + w/2, values_ma, w, color='#D6604D', edgecolor='white', label='MH63 vs NONA')

# Highlight CTB4a
ax_b.annotate('CONSERVED\n(<0.01)', xy=(0, 0.0076), xytext=(0.2, 0.10), fontsize=8, ha='center',
              arrowprops=dict(arrowstyle='->', color='#27AE60', lw=1.5), color='#27AE60', fontweight='bold')
# Highlight Peak 1
ax_b.annotate('25.6×\nhigher', xy=(2, 0.1946), xytext=(2.5, 0.22), fontsize=8, ha='center',
              arrowprops=dict(arrowstyle='->', color='#B2182B', lw=1.5), color='#B2182B', fontweight='bold')

ax_b.set_xticks(x); ax_b.set_xticklabels(regions, fontsize=6.5)
ax_b.set_ylabel('OGR embedding distance'); ax_b.set_ylim(0, 0.28)
ax_b.set_title('B  CTB4a Gene vs TE Clusters: Embedding Divergence Comparison', loc='left', fontweight='bold', fontsize=9)
ax_b.legend(fontsize=7)

# Panel C: TE effect size vs embedding distance
ax_c = fig.add_axes([0.56, 0.06, 0.40, 0.42])

# Per-window: compute correlation between TE count and embedding distance
te_counts = [w['te_count'] for w in sw]
max_dists = [w['d_max'] for w in sw]
ax_c.scatter(te_counts, max_dists, c=np.array(te_counts)>0, cmap='Reds', s=40, edgecolors='white', alpha=0.8)
ax_c.set_xlabel('Number of TEs in 20kb window'); ax_c.set_ylabel('Max embedding distance')
ax_c.set_title('C  TE Density vs Embedding Divergence', loc='left', fontweight='bold', fontsize=9)

# Highlight CTB4a point (2 TEs, low distance)
ax_c.annotate('CTB4a\n(2 TEs,\nlow div.)', xy=(2, 0.0076), xytext=(4, 0.03), fontsize=7, ha='center',
              arrowprops=dict(arrowstyle='->', color='#27AE60', lw=1), color='#27AE60')
# Highlight peak point
ax_c.annotate('Peak 1\n(3 TEs,\nhigh div.)', xy=(3, 0.1946), xytext=(5, 0.17), fontsize=7, ha='center',
              arrowprops=dict(arrowstyle='->', color='#B2182B', lw=1), color='#B2182B')

# Regression line
from scipy.stats import linregress
if len(te_counts) > 2:
    slope, intercept, r, p, _ = linregress(te_counts, max_dists)
    xs_line = np.array([0, 7])
    ax_c.plot(xs_line, slope*xs_line + intercept, '--', color='#999', linewidth=0.8)
    ax_c.text(5, 0.02, f'r={r:.2f}, p={p:.3f}', fontsize=7, color='#999')

# Panel D: Fine-grained detail around Peak 1 (2.1Mb ±50kb)
ax_d = fig.add_axes([0.56, 0.55, 0.40, 0.40])
if fine:
    # Filter for the 2.1Mb center
    peak1_fine = [f for f in fine if abs(f['center'] - 2100000) < 100000]
    if not peak1_fine:
        peak1_fine = [f for f in fine if abs(f['center'] - fine[0]['center']) < 100000]

    if peak1_fine:
        fx = [f['pos']/1e6 for f in peak1_fine]
        fmn = [f['d_MH63_NIP'] for f in peak1_fine]
        fma = [f['d_MH63_NONA'] for f in peak1_fine]
        ax_d.plot(fx, fmn, 'o-', color='#888888', markersize=4, linewidth=1.5)
        ax_d.plot(fx, fma, 's-', color='#D6604D', markersize=4, linewidth=1.5)
        ax_d.fill_between(fx, fmn, fma, alpha=0.1, color='#D6604D')
        ax_d.set_xlabel('Chr04 position (Mb)')
        ax_d.set_ylabel('Cosine distance')
        ax_d.set_title('D  Fine Structure Around Peak 1 (±50kb, 5kb step)', loc='left', fontweight='bold', fontsize=9)
        # Mark TE positions
        for te in te_positions:
            if 2.05 < te['pos']/1e6 < 2.20:
                ax_d.axvline(x=te['pos']/1e6, color=te['color'], alpha=0.4, linewidth=0.5)

fig.savefig(os.path.join(OUT, 'P5_OGR_SlidingWindow.png'))
fig.savefig(os.path.join(OUT, 'P5_OGR_SlidingWindow.svg'))
plt.close()
print('P5 Figure saved!')
