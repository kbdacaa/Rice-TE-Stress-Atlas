"""NSFC技术路线图 — 纯计算方向"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300

fig, ax = plt.subplots(1, 1, figsize=(16, 12))
ax.set_xlim(0, 16)
ax.set_ylim(0, 12)
ax.axis('off')

# Color scheme
colors = {
    'data': '#3498DB',
    'compute': '#2ECC71',
    'model': '#E67E22',
    'output': '#9B59B6',
    'arrow': '#7F8C8D',
}

def draw_box(ax, x, y, w, h, text, color, fontsize=9, bold=False):
    """Draw a rounded box with text"""
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3",
                                     facecolor=color, edgecolor='white', alpha=0.9, linewidth=2)
    ax.add_patch(rect)
    weight = 'bold' if bold else 'normal'
    lines = text.split('\n')
    line_h = h / max(len(lines), 1)
    for i, line in enumerate(lines):
        ax.text(x + w/2, y + h - line_h/2 - i*line_h, line, ha='center', va='center',
                fontsize=fontsize, fontweight=weight, color='white')

def draw_arrow(ax, x1, y1, x2, y2, color='#7F8C8D', lw=2):
    """Draw an arrow"""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw))

def draw_label(ax, x, y, text, fontsize=8, color='#2C3E50'):
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize, color=color)

# ═══ Title ═══
ax.text(8, 11.5, '水稻骨干亲本TE多态性与多胁迫适应谱的关联图谱构建及调控网络解析',
        ha='center', fontsize=16, fontweight='bold', color='#2C3E50')
ax.text(8, 11.0, '技术路线图', ha='center', fontsize=12, color='#7F8C8D')

# ═══ Layer 0: Data Sources (top) ═══
data_y = 9.5
draw_box(ax, 0.3, data_y, 3.2, 1.0, '数据层\nRICEPTEDB 12品种基因组+TE+基因注释\n26,697 TIPs × 7胁迫GWAS关联\nRPRP M16 泛基因组图+16品种注释', colors['data'])
draw_box(ax, 4.0, data_y, 3.2, 1.0, '公共多组学数据\n8个水稻胁迫表达数据集(GEO/SRA)\nENCODE水稻染色质状态数据\n文献中CTB2/CTB4a/DGK7功能验证', colors['data'])
draw_box(ax, 7.7, data_y, 3.2, 1.0, '品种资源\n12骨干亲本(IRRI MAGIC)\nRPRP 16品种独立验证群体\n165品种冷存活率表型', colors['data'])
draw_box(ax, 11.4, data_y, 4.2, 1.0, '前期积累\n138冷基因×TE交叉\n7×7胁迫TE方向矩阵\nNONA_BOKRA 728独占TE\n2个TE热点单倍型图谱', colors['data'])

# Layer 0 → Layer 1 arrows
for x_center in [2.0, 5.6, 9.3, 13.5]:
    draw_arrow(ax, x_center, data_y, x_center, data_y - 0.7, colors['arrow'])

# ═══ Layer 1: Computational Analysis (middle) ═══
comp_y1 = 5.5
comp_y2 = 3.2

# Year 1
draw_box(ax, 0.3, comp_y1, 3.5, 0.9, '第1年: TE胁迫关联图谱\n构建12品种TE泛基因组图谱\nTE×7胁迫关联效应全矩阵\n公共ChIP/ATAC/RNA调控注释', colors['compute'], bold=True)

# Year 2
draw_box(ax, 4.2, comp_y1, 3.5, 0.9, '第2年: 调控网络与多效性\nTE-基因-胁迫三层网络建模\nTE多效性系统分类(协同/权衡)\nCTB2/CTB4a-DGK7双热点精细解析', colors['compute'], bold=True)

# Year 3
draw_box(ax, 8.1, comp_y1, 3.5, 0.9, '第3年: 预测模型与群体验证\n广适性TE预测模型(ML)\nRPRP独立群体验证\n骨干亲本TE标记体系建立', colors['compute'], bold=True)

# Year 4
draw_box(ax, 12.0, comp_y1, 3.5, 0.9, '第4年: 数据库与成果固化\nTE-Stress Atlas公开数据库\n系列论文撰写\n育种推荐工具原型', colors['compute'], bold=True)

# Cross connections between years
draw_arrow(ax, 3.8, comp_y1 + 0.45, 4.2, comp_y1 + 0.45, colors['arrow'])
draw_arrow(ax, 7.7, comp_y1 + 0.45, 8.1, comp_y1 + 0.45, colors['arrow'])
draw_arrow(ax, 11.6, comp_y1 + 0.45, 12.0, comp_y1 + 0.45, colors['arrow'])

# Year 1 → 2 → 3 → 4 bidirectional data flow
for x1, x2 in [(3.8, 4.2), (7.7, 8.1), (11.6, 12.0)]:
    draw_arrow(ax, x2, comp_y1 + 0.2, x1, comp_y1 + 0.2, '#BDC3C7', lw=1)

# Layer 1 → Layer 2 arrows (each year has sub-outputs)
sub_y = comp_y1 - 1.5
sub_boxes = [
    (0.3, 'TE-Stress Atlas v1.0\nTE-PAV矩阵\nTE调控潜力注释'),
    (4.2, '三层调控网络\nTE多效性分类体系\nChr02/Chr04精细模型'),
    (8.1, '广适性ML预测模型\n诊断性TE标记组合\n跨群体验证报告'),
    (12.0, '在线数据库\n3-4篇SCI论文\nTE标记使用规范'),
]

for i, (x, text) in enumerate(sub_boxes):
    draw_box(ax, x, sub_y, 3.5, 0.7, text, colors['model'])
    draw_arrow(ax, x + 1.75, comp_y1, x + 1.75, sub_y + 0.7, colors['arrow'])

# All sub-boxes → final output
draw_arrow(ax, 8, sub_y - 0.1, 8, sub_y - 0.8, colors['arrow'])

# ═══ Layer 2: Final Output ═══
final_y = 0.5
draw_box(ax, 2.0, final_y, 12.0, 1.2,
         '最终产出\n'
         'TE-Stress Atlas公开数据库 | TE标记体系(10-15个诊断标记) | 广适性预测模型 | 4-5篇SCI论文 | 育种亲本选配推荐工具',
         colors['output'], bold=True, fontsize=11)

# ═══ Side annotations ═══
ax.text(15.6, 7.5, '数\n据\n驱\n动', fontsize=10, color='#3498DB', ha='center')
ax.text(15.6, 4.5, '计\n算\n建\n模', fontsize=10, color='#2ECC71', ha='center')
ax.text(15.6, 2.0, '应\n用\n产\n出', fontsize=10, color='#9B59B6', ha='center')

# ═══ Legend at bottom ═══
legend_y = 0.1
legend_items = [
    (colors['data'], '数据基础'),
    (colors['compute'], '计算分析'),
    (colors['model'], '阶段产出'),
    (colors['output'], '最终成果'),
]
for i, (color, label) in enumerate(legend_items):
    x = 2.0 + i * 3.5
    rect = mpatches.FancyBboxPatch((x, legend_y), 0.3, 0.2, boxstyle="round,pad=0.1",
                                     facecolor=color, edgecolor='white')
    ax.add_patch(rect)
    ax.text(x + 0.5, legend_y + 0.1, label, fontsize=8, va='center')

plt.tight_layout()
out_dir = r'D:\project\AutoResearch-rice-T2T\data\figures'
import os
os.makedirs(out_dir, exist_ok=True)
fig.savefig(os.path.join(out_dir, 'Fig6_Technical_Roadmap.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(out_dir, 'Fig6_Technical_Roadmap.svg'), bbox_inches='tight')
plt.close()
print('技术路线图已保存: Fig6_Technical_Roadmap.png / .svg')
