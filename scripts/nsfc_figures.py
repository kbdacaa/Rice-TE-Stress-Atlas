"""NSFC申报书5张核心Figure生成脚本
Fig1: 12品种×7胁迫TE风险热图
Fig2: Chr02+Chr04冷TE簇与基因共定位
Fig3: 冷-病品种排名反转散点图
Fig4: CTB2/CTB4a演化时间线+TE分布
Fig5: 7胁迫TE多效性UpSet图

输出: data/figures/ (PNG, 300dpi) + SVG矢量图
"""
import json, sqlite3, os
from collections import defaultdict, Counter
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300

DB = r"D:\project\AutoResearch-rice-T2T\data\rice_integration.db"
OUT_DIR = r"D:\project\AutoResearch-rice-T2T\data\figures"
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = ['MH63','T','Y','M','N','NJ','KO','C','K','R','TB','L']
SAMPLE_LABELS = {
    'MH63':'明恢63','T':'特青','Y':'9311','M':'MADINIKA','N':'NONA_BOKRA',
    'NJ':'南京11','KO':'KOGONI','C':'CICA','K':'越光','R':'日本晴','TB':'Fujisaka5','L':'Lemont'
}
SUBPOP_COLORS = {
    'indica': '#E74C3C',
    'japonica_temp': '#2980B9',
    'japonica_trop': '#1ABC9C',
    'aus': '#F39C12',
    'africa': '#8E44AD',
    'irri_breed': '#95A5A6',
}
SAMPLE_SUBPOP = {
    'MH63':'indica','T':'indica','Y':'indica','NJ':'indica',
    'K':'japonica_temp','R':'japonica_temp','TB':'japonica_temp',
    'L':'japonica_trop','N':'aus','M':'indica','KO':'africa','C':'irri_breed',
}
STRESSES = ['cold','disease','drought','heat','salinity','submergence','heavy_metal']
STRESS_LABELS = {
    'cold':'冷','disease':'病','drought':'旱','heat':'热',
    'salinity':'盐','submergence':'淹','heavy_metal':'重金属'
}
STRESS_COLORS = {
    'cold':'#3498DB','disease':'#E74C3C','drought':'#F39C12','heat':'#E67E22',
    'salinity':'#1ABC9C','submergence':'#9B59B6','heavy_metal':'#7F8C8D'
}

# ═══════════════════════════════════════════════════════════════
def load_stress_risk_data():
    """Load all-stress risk profile data"""
    fpath = r"D:\project\AutoResearch-rice-T2T\data\all_stress_risk_profile.json"
    if os.path.exists(fpath):
        with open(fpath, 'r', encoding='utf-8') as f:
            return json.load(f)
    return None

def load_haplotype_data():
    fpath = r"D:\project\AutoResearch-rice-T2T\data\te_haplotype_analysis.json"
    if os.path.exists(fpath):
        with open(fpath, 'r', encoding='utf-8') as f:
            return json.load(f)
    return None

def load_ctb4a_data():
    fpath = r"D:\project\AutoResearch-rice-T2T\data\ctb4a_te_landscape.json"
    if os.path.exists(fpath):
        with open(fpath, 'r', encoding='utf-8') as f:
            return json.load(f)
    return None


# ═══════════════════════════════════════════════════════════════
# Fig1: 12品种×7胁迫TE风险热图
# ═══════════════════════════════════════════════════════════════
def fig1_stress_risk_heatmap():
    """12-variety × 7-stress TE risk score heatmap"""
    data = load_stress_risk_data()
    if not data:
        print("[FIG1] 数据不可用，跳过")
        return

    # Build risk matrix
    risk_matrix = np.zeros((12, 7))
    for i, s in enumerate(SAMPLES):
        for j, st in enumerate(STRESSES):
            # Try to extract from data
            if 'variety_risk' in data:
                vr = data['variety_risk'].get(s, {})
                risk_matrix[i, j] = vr.get(st, 0)

    # Normalize per stress for visualization
    norm_matrix = np.zeros_like(risk_matrix)
    for j in range(7):
        col = risk_matrix[:, j]
        rng = np.ptp(col)
        if rng > 0:
            norm_matrix[:, j] = (col - np.min(col)) / rng
        else:
            norm_matrix[:, j] = 0.5

    # Sort rows by avg risk
    avg_risk = norm_matrix.mean(axis=1)
    order = np.argsort(avg_risk)
    risk_matrix_sorted = risk_matrix[order]
    samples_sorted = [SAMPLES[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, 8))

    cmap = plt.cm.RdBu_r
    im = ax.imshow(norm_matrix[order], aspect='auto', cmap=cmap, vmin=0, vmax=1)

    # Annotate cells with risk score
    for i in range(12):
        for j in range(7):
            val = risk_matrix_sorted[i, j]
            color = 'white' if 0.3 < norm_matrix[order][i, j] < 0.7 else 'black'
            ax.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=7,
                    color=color, fontweight='bold')

    ax.set_xticks(range(7))
    ax.set_xticklabels([STRESS_LABELS[s] for s in STRESSES], fontsize=11)
    ax.set_yticks(range(12))
    ax.set_yticklabels([f"{SAMPLE_LABELS.get(s,s)}({s})" for s in samples_sorted], fontsize=9)

    # Color y-tick labels by subpopulation
    for i, s in enumerate(samples_sorted):
        sp = SAMPLE_SUBPOP.get(s, 'indica')
        color = SUBPOP_COLORS.get(sp, 'black')
        ax.get_yticklabels()[i].set_color(color)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('标准化TE风险评分', fontsize=10)

    ax.set_title('12骨干亲本 × 7种胁迫 TE风险谱', fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('胁迫类型', fontsize=12)
    ax.set_ylabel('品种（按TE平均风险排序）', fontsize=12)

    # Legend for subpopulations
    legend_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['indica'], label='籼稻(indica)', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['japonica_temp'], label='温带粳稻(japonica)', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['japonica_trop'], label='热带粳稻', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['aus'], label='aus', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['africa'], label='非洲稻', markersize=10),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='lower right', ncol=2)

    # Annotations
    ax.annotate('NONA_BOKRA\n冷TE保护型', xy=(0, 0), xytext=(-1.8, -1.5),
                fontsize=8, arrowprops=dict(arrowstyle='->', color='red'), color='red')
    ax.annotate('越光·日本晴\n病TE保护型',
                xy=(1, list(samples_sorted).index('K')),
                xytext=(8.2, list(samples_sorted).index('K')-0.5),
                fontsize=8, arrowprops=dict(arrowstyle='->', color='blue'), color='blue')

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'Fig1_12varieties_7stress_risk_heatmap.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, 'Fig1_12varieties_7stress_risk_heatmap.svg'), bbox_inches='tight')
    plt.close()
    print("[FIG1] 完成: 12品种×7胁迫TE风险热图")


# ═══════════════════════════════════════════════════════════════
# Fig2: Chr02 + Chr04 Cold TE cluster co-localization
# ═══════════════════════════════════════════════════════════════
def fig2_cold_te_colocalization():
    """Two-panel: Chr04 and Chr02 cold TE clusters with gene annotations"""
    conn = sqlite3.connect(DB)
    c = conn.cursor()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios':[1,1]})

    # Panel A: Chr04 CTB2/CTB4a region
    te_data_chr04 = c.execute(
        "SELECT ti.position, ta.odds_ratio, ta.category, ta.direction "
        "FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id "
        "WHERE ti.chromosome='Chr04' AND ti.position BETWEEN 1100000 AND 5500000 "
        "ORDER BY ti.position"
    ).fetchall()

    positions_04 = [t[0] for t in te_data_chr04]
    ors_04 = [t[1] for t in te_data_chr04]
    cats_04 = [t[2] for t in te_data_chr04]

    colors_04 = []
    sizes_04 = []
    for cat, or_val in zip(cats_04, ors_04):
        if cat == 'strong':
            colors_04.append('#E74C3C' if or_val > 1 else '#2980B9')
            sizes_04.append(60)
        elif cat == 'trend':
            colors_04.append('#F1948A' if or_val > 1 else '#85C1E9')
            sizes_04.append(30)
        else:
            colors_04.append('#E0E0E0')
            sizes_04.append(10)

    ax1.scatter(positions_04, [0]*len(positions_04), c=colors_04, s=sizes_04, alpha=0.7, edgecolors='none')
    ax1.set_xlim(1100000, 5500000)
    ax1.set_ylim(-0.5, 0.5)

    # Mark CTB2/CTB4a 56.8kb region
    ctb_region_start = 1330000
    ctb_region_end = 1342000
    ax1.axvspan(ctb_region_start, ctb_region_end, alpha=0.3, color='#F1C40F', label='CTB2/CTB4a 56.8kb (Li 2021)')
    ax1.annotate('CTB2\nCTB4a\n56.8kb', xy=(1336000, 0.35), fontsize=9, ha='center',
                 bbox=dict(boxstyle='round', facecolor='#F1C40F', alpha=0.7))

    # Mark TE#8279 (strongest cold protection)
    ax1.annotate('TE#8279\n(OR=0.048)\n最强冷保护', xy=(1485066, 0.1), xytext=(1550000, 0.35),
                 fontsize=8, arrowprops=dict(arrowstyle='->', color='green'), color='green')

    ax1.set_title('Chr04: CTB2/CTB4a区域冷TE簇 (1.1-5.5Mb, 32个冷显著TE)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('', fontsize=10)
    ax1.set_yticks([])
    ax1.legend(fontsize=9, loc='upper right')

    # Panel B: Chr02 DGK7 region
    te_data_chr02 = c.execute(
        "SELECT ti.position, ta.odds_ratio, ta.category, ta.direction "
        "FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id "
        "WHERE ti.chromosome='Chr02' AND ti.position BETWEEN 34500000 AND 36400000 "
        "ORDER BY ti.position"
    ).fetchall()

    positions_02 = [t[0] for t in te_data_chr02]
    ors_02 = [t[1] for t in te_data_chr02]
    cats_02 = [t[2] for t in te_data_chr02]

    colors_02 = []
    sizes_02 = []
    for cat, or_val in zip(cats_02, ors_02):
        if cat == 'strong':
            colors_02.append('#E74C3C' if or_val > 1 else '#2980B9')
            sizes_02.append(60)
        elif cat == 'trend':
            colors_02.append('#F1948A' if or_val > 1 else '#85C1E9')
            sizes_02.append(30)
        else:
            colors_02.append('#E0E0E0')
            sizes_02.append(10)

    ax2.scatter(positions_02, [0]*len(positions_02), c=colors_02, s=sizes_02, alpha=0.7, edgecolors='none')
    ax2.set_xlim(34500000, 36400000)
    ax2.set_ylim(-0.5, 0.5)

    # Mark DGK7
    dgk7_start = 34887557
    dgk7_end = 34894178
    ax2.axvspan(dgk7_start, dgk7_end, alpha=0.3, color='#E67E22', label='DGK7 (Cell 2025)')
    ax2.annotate('DGK7\n热感知基因\n(Cell 2025)', xy=(34890867, 0.35), fontsize=9, ha='center',
                 bbox=dict(boxstyle='round', facecolor='#E67E22', alpha=0.7))

    # Mark key cold TEs near DGK7
    for te_pos, te_label in [(34927145, 'TE#5625'), (34971776, 'TE#5630')]:
        ax2.annotate(f'{te_label}\n冷保护', xy=(te_pos, 0.15), xytext=(te_pos+300000, 0.35),
                     fontsize=8, arrowprops=dict(arrowstyle='->', color='blue'), color='blue')

    ax2.set_title('Chr02: DGK7区域冷TE簇 (34.5-36.4Mb, 18个冷显著TE)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('基因组位置 (Mb)', fontsize=11)
    ax2.set_yticks([])
    ax2.legend(fontsize=9, loc='upper right')

    # Legend for colors
    legend_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#E74C3C', label='冷风险(OR>1)', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#2980B9', label='冷保护(OR<1)', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#E0E0E0', label='无显著关联', markersize=10),
    ]

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'Fig2_cold_TE_colocalization.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, 'Fig2_cold_TE_colocalization.svg'), bbox_inches='tight')
    plt.close()
    conn.close()
    print("[FIG2] 完成: Chr02+Chr04冷TE簇与基因共定位")


# ═══════════════════════════════════════════════════════════════
# Fig3: Cold-Disease Ranking Inversion Scatter
# ═══════════════════════════════════════════════════════════════
def fig3_cold_disease_ranking_inversion():
    """Scatter plot showing cold vs disease TE ranking inversion"""
    # Hardcoded from evidence_summary.txt section 2.10
    rankings = {
        'Lemont':       (9, 5),
        'CICA':         (3, 2),
        'MADINIKA':     (2, 6),
        'Koshihikari':  (7, 1),
        'Fujisaka5':    (6, 3),
        'Minghui63':    (4, 9),
        'Nipponbare':   (8, 4),
        '9311':         (12, 12),
        'NONA_BOKRA':   (1, 7),
        'KOGONI':       (10, 8),
        'Nanjing11':    (5, 10),
        'TEQING':       (11, 11),
    }

    fig, ax = plt.subplots(figsize=(8, 7))

    for sample_name, (cold_rank, disease_rank) in rankings.items():
        # Find sample code
        s_code = None
        for sc, sl in SAMPLE_LABELS.items():
            if sl in sample_name or sample_name in sl:
                s_code = sc
                break
        if not s_code:
            for sc in SAMPLES:
                if sample_name.lower().replace(' ','_') == sc.lower().replace(' ','_'):
                    s_code = sc
                    break

        sp = SAMPLE_SUBPOP.get(s_code, 'indica') if s_code else 'indica'
        color = SUBPOP_COLORS.get(sp, 'black')

        ax.scatter(cold_rank, disease_rank, s=200, c=color, edgecolors='black', linewidths=1.5, zorder=5)
        label = SAMPLE_LABELS.get(s_code, sample_name) if s_code else sample_name
        ax.annotate(label, (cold_rank + 0.3, disease_rank + 0.1), fontsize=9)

    # Diagonal line
    ax.plot([1, 12], [1, 12], 'k--', alpha=0.3, linewidth=1)
    ax.annotate('冷=病\n(完全一致)', xy=(3, 3.5), fontsize=8, color='gray', rotation=45)

    # Trade-off zones
    ax.fill_between([1, 12], [7, 12], [12, 12], alpha=0.1, color='red')
    ax.fill_between([7, 12], [1, 12], [1, 6], alpha=0.1, color='blue')
    ax.annotate('冷保护+病害风险\n(NONA_BOKRA型)', xy=(2, 10.5), fontsize=9, color='red')
    ax.annotate('病害保护+冷风险\n(粳稻型)', xy=(9, 2), fontsize=9, color='blue')

    ax.set_xlabel('冷胁迫TE风险排名 (1=最保护)', fontsize=12)
    ax.set_ylabel('病害TE风险排名 (1=最保护)', fontsize=12)
    ax.set_title('冷-病害TE保护排名反转', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 13.5)
    ax.invert_yaxis()
    ax.invert_xaxis()  # 1=best, upper-left

    # Legend
    legend_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['indica'], label='籼稻', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['japonica_temp'], label='温带粳稻', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['japonica_trop'], label='热带粳稻', markersize=10),
        Line2D([0],[0], marker='o', color='w', markerfacecolor=SUBPOP_COLORS['aus'], label='aus', markersize=10),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc='lower left')

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'Fig3_cold_disease_ranking_inversion.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, 'Fig3_cold_disease_ranking_inversion.svg'), bbox_inches='tight')
    plt.close()
    print("[FIG3] 完成: 冷-病排名反转散点图")


# ═══════════════════════════════════════════════════════════════
# Fig4: CTB2/CTB4a Evolutionary Timeline + TE Distribution
# ═══════════════════════════════════════════════════════════════
def fig4_ctb2_ctb4a_evolution():
    """Evolutionary timeline + TE distribution at CTB2/CTB4a 56.8kb region"""
    fig = plt.figure(figsize=(14, 10))

    # Panel A: Evolutionary timeline (top)
    ax_time = fig.add_axes([0.1, 0.72, 0.85, 0.25])
    ax_time.set_xlim(8000, 0)
    ax_time.set_ylim(-1, 5)

    # Timeline
    events = [
        (7000, 2.5, '普通野生稻\nCTB2 SNP(A)\n(standing variation)', '#27AE60'),
        (6400, 1.5, '云南低海拔(<1600m)\nCTB2首先被选择', '#2980B9'),
        (3200, 0.5, '云南中高海拔(1600-2100m)\nCTB4a de novo mutation\n(4个功能SNP)', '#E74C3C'),
        (500, -0.5, '扩散至东北/日韩\n高纬度冷区\n联合单倍型固定', '#8E44AD'),
    ]

    y_positions = []
    for yr, y, label, color in events:
        ax_time.plot([yr, yr], [-0.8, 4.5], '--', color='gray', alpha=0.3, linewidth=0.8)
        ax_time.scatter(yr, y, s=300, c=color, edgecolors='black', linewidths=2, zorder=5)
        ax_time.annotate(label, (yr, y), textcoords="offset points", xytext=(0, 25 if y < 2 else -35),
                         fontsize=9, ha='center', va='center',
                         bbox=dict(boxstyle='round', facecolor=color, alpha=0.15))
        y_positions.append(y)

    ax_time.set_xlabel('距今年数 (yr BP)', fontsize=11)
    ax_time.set_title('CTB2/CTB4a逐步选择促进粳稻冷适应演化 (Li et al. 2021, New Phytologist)', fontsize=13, fontweight='bold')
    ax_time.set_yticks([])
    ax_time.spines['left'].set_visible(False)

    # Arrow of time
    ax_time.annotate('', xy=(8000, 4.2), xytext=(0, 4.2),
                     arrowprops=dict(arrowstyle='->', lw=2, color='gray'))

    # Panel B: TE genotype matrix at CTB4a core region (bottom-left)
    conn = sqlite3.connect(DB)
    c = conn.cursor()

    tes = c.execute(
        "SELECT id, position, vcf_id, te_type FROM te_insertions "
        "WHERE chromosome='Chr04' AND position BETWEEN 1310000 AND 1400000 ORDER BY position"
    ).fetchall()

    if tes:
        ax_matrix = fig.add_axes([0.08, 0.08, 0.55, 0.55])
        n_tes = len(tes)

        # Build genotype matrix from VCF
        VCF = r"D:\fanjiyinzu\RICEPTEDB\rice.pTE.vcf"
        vcf_data = {}
        with open(VCF, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if parts[0] != 'Chr04': continue
                pos = int(parts[1])
                if 1310000 <= pos <= 1400000:
                    vcf_data[pos] = {
                        'id': parts[2],
                        'gts': parts[9:] if len(parts) > 9 else []
                    }

        matrix = np.zeros((n_tes, 12))
        for i, (te_id, pos, vcf_id, te_type) in enumerate(tes):
            vd = vcf_data.get(pos)
            if vd:
                for j, s in enumerate(SAMPLES):
                    if j < len(vd['gts']):
                        gt = vd['gts'][j]
                        if gt in ('0/0','0|0','0'): matrix[i,j] = 0
                        elif gt in ('0/1','1/0','0|1','1|0'): matrix[i,j] = 1
                        elif gt in ('1/1','1|1','1'): matrix[i,j] = 2
                        else: matrix[i,j] = -1

        im = ax_matrix.imshow(matrix, aspect='auto', cmap=plt.cm.RdYlBu_r, vmin=0, vmax=2)

        ax_matrix.set_yticks(range(n_tes))
        ax_matrix.set_yticklabels([t[2].replace('OsTIP04G00','') for t in tes], fontsize=6)
        ax_matrix.set_xticks(range(12))
        ax_matrix.set_xticklabels([f"{SAMPLE_LABELS.get(s,s)}" for s in SAMPLES], fontsize=7, rotation=45, ha='right')

        # Color x-tick labels
        for j, s in enumerate(SAMPLES):
            sp = SAMPLE_SUBPOP.get(s, 'indica')
            ax_matrix.get_xticklabels()[j].set_color(SUBPOP_COLORS.get(sp, 'black'))

        ax_matrix.set_title(f'CTB4a核心区域 TE基因型矩阵 ({n_tes}个TE)', fontsize=11, fontweight='bold')
        ax_matrix.set_xlabel('品种', fontsize=10)

        # Mark CTB4a gene position
        ctb4a_te_positions = [1334699, 1336356]  # TE#8261, TE#8262
        for pos in ctb4a_te_positions:
            te_idx = None
            for i, t in enumerate(tes):
                if abs(t[1] - pos) < 1000:
                    te_idx = i
                    break
            if te_idx is not None:
                ax_matrix.axhline(y=te_idx, color='yellow', linewidth=2, alpha=0.7)
                ax_matrix.annotate(f'CTB4a基因\n内含子TE', xy=(11.5, te_idx), fontsize=7, color='red',
                                   xytext=(13, te_idx), arrowprops=dict(arrowstyle='->', color='red'))

    # Panel C: Legend & annotations (bottom-right)
    ax_legend = fig.add_axes([0.68, 0.08, 0.28, 0.55])
    ax_legend.axis('off')

    annotations = [
        "CTB2 (甾醇糖基转移酶)",
        "  → 甾醇糖苷 → 膜脂重塑",
        "CTB4a (LRR-RLK)",
        "  → ATP合酶 → 能量代谢",
        "56.8kb共定位 (Li et al. 2021)",
        "  → 联合单倍型",
        "",
        "TE基因型图例:",
        "  深红 = AA (替代等位纯合)",
        "  白色 = RA (杂合)",
        "  深蓝 = RR (参考等位纯合)",
        "",
        "★ TE#8262: CtB4a基因内含子",
        "    粳稻固定插入(亚群分化)",
        "★ TE#8279: 最强冷保护信号",
        "    CICA/NONA_BOKRA独占",
    ]
    for i, ann in enumerate(annotations):
        ax_legend.text(0, 0.95 - i*0.07, ann, transform=ax_legend.transAxes, fontsize=8,
                       va='top')

    conn.close()
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'Fig4_CTB2_CTB4a_evolution_TE_distribution.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, 'Fig4_CTB2_CTB4a_evolution_TE_distribution.svg'), bbox_inches='tight')
    plt.close()
    print("[FIG4] 完成: CTB2/CTB4a演化时间线+TE分布")


# ═══════════════════════════════════════════════════════════════
# Fig5: Multi-stress TE pleiotropy UpSet plot
# ═══════════════════════════════════════════════════════════════
def fig5_multistress_te_upset():
    """UpSet-like plot showing TE overlap across 7 stresses"""
    conn = sqlite3.connect(DB)
    c = conn.cursor()

    # For each stress, get significant TE IDs (category='strong' or OR significant)
    sig_tes = {}
    for stress in STRESSES:
        tbl = f"te_{stress}_association"
        tes = c.execute(f"SELECT te_id FROM {tbl} WHERE category='strong'").fetchall()
        sig_tes[stress] = set(t[0] for t in tes)
        print(f"  {stress}: {len(sig_tes[stress])} significant TEs")

    # Count overlaps
    from itertools import combinations

    # Build UpSet-style matrix
    stress_list = STRESSES
    n_stress = len(stress_list)

    # For each TE, which stresses is it significant in?
    te_membership = defaultdict(int)
    for te_id in range(1, 26698):
        mask = 0
        for j, st in enumerate(stress_list):
            if te_id in sig_tes[st]:
                mask |= (1 << j)
        if mask > 0:
            te_membership[mask] += 1

    # Get top intersection sets
    intersections = []
    for r in range(1, n_stress + 1):
        for combo in combinations(range(n_stress), r):
            mask = sum(1 << i for i in combo)
            count = te_membership.get(mask, 0)
            if count > 0:
                stress_names = [stress_list[i] for i in combo]
                intersections.append((stress_names, count))

    intersections.sort(key=lambda x: -x[1])
    top_n = min(15, len(intersections))
    intersections = intersections[:top_n]

    if not intersections:
        print("[FIG5] 无TE重叠数据，跳过")
        conn.close()
        return

    fig, (ax_bar, ax_matrix) = plt.subplots(2, 1, figsize=(14, 7),
                                              gridspec_kw={'height_ratios': [3, 2]}, sharex=True)

    # Bar chart of intersection sizes
    x_pos = range(len(intersections))
    sizes = [isc[1] for isc in intersections]
    bars = ax_bar.bar(x_pos, sizes, color=['#34495E' if len(isc[0])>=3 else '#95A5A6' for isc in intersections])

    for i, (isc, bar) in enumerate(zip(intersections, bars)):
        if isc[1] > 0:
            ax_bar.text(i, bar.get_height() + max(sizes)*0.02, str(isc[1]),
                       ha='center', fontsize=8, fontweight='bold')

    ax_bar.set_ylabel('TE数量', fontsize=11)
    ax_bar.set_title('7种胁迫间TE多效性重叠', fontsize=14, fontweight='bold')
    ax_bar.spines['top'].set_visible(False)

    # Matrix: which stresses in each intersection
    matrix_data = np.zeros((len(intersections), n_stress))
    for i, (stress_names, count) in enumerate(intersections):
        for st in stress_names:
            j = stress_list.index(st)
            matrix_data[i, j] = 1

    ax_matrix.imshow(matrix_data, aspect='auto', cmap=plt.cm.Blues, alpha=0.7)
    ax_matrix.set_yticks(range(len(intersections)))
    ax_matrix.set_yticklabels([f"{'+'.join(STRESS_LABELS[s][:2] for s in st)}" for st, _ in intersections], fontsize=7)
    ax_matrix.set_xticks(range(n_stress))
    ax_matrix.set_xticklabels([STRESS_LABELS[s] for s in stress_list], fontsize=10)

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'Fig5_multistress_TE_upset.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, 'Fig5_multistress_TE_upset.svg'), bbox_inches='tight')
    plt.close()
    conn.close()
    print("[FIG5] 完成: 7胁迫TE多效性UpSet图")


# ═══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("=" * 60)
    print("NSFC申报书 Figure生成")
    print("=" * 60)

    fig1_stress_risk_heatmap()
    fig2_cold_te_colocalization()
    fig3_cold_disease_ranking_inversion()
    fig4_ctb2_ctb4a_evolution()
    fig5_multistress_te_upset()

    print(f"\n所有Figure已保存至: {OUT_DIR}")
