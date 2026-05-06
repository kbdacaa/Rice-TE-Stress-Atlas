"""研究假说支撑实验：三项独立验证
1. 7×7胁迫间TE方向一致性矩阵 → 胁迫协同vs权衡的分子基础
2. NONA_BOKRA全基因组独占TE系统鉴定 → aus亚群独特TE保护谱
3. TE超家族胁迫关联富集分析 → TE类型的功能特异性
"""
import sqlite3, json, os
from collections import defaultdict, Counter
from itertools import combinations
from scipy.stats import chi2_contingency, fisher_exact

DB = r"D:\project\AutoResearch-rice-T2T\data\rice_integration.db"
VCF = r"D:\fanjiyinzu\RICEPTEDB\rice.pTE.vcf"
OUT_DIR = r"D:\project\AutoResearch-rice-T2T\data"

STRESSES = ['cold','disease','drought','heat','salinity','submergence','heavy_metal']
STRESS_LABELS = {'cold':'冷','disease':'病','drought':'旱','heat':'热',
                 'salinity':'盐','submergence':'淹','heavy_metal':'重金属'}
SAMPLES = ['MH63','T','Y','M','N','NJ','KO','C','K','R','TB','L']

conn = sqlite3.connect(DB)
c = conn.cursor()

# ═══════════════════════════════════════════════════════════════
# 1. 7×7 Cross-Stress TE Direction Concordance Matrix
# ═══════════════════════════════════════════════════════════════
print("=" * 80)
print("分析1: 7×7胁迫间TE方向一致性矩阵")
print("=" * 80)

# For each stress pair, find TEs significant in both
concordance_matrix = {}
pair_details = {}

for s1, s2 in combinations(STRESSES, 2):
    # Get significant TEs for each stress
    sig1 = {}
    for r in c.execute(f"SELECT te_id, direction FROM te_{s1}_association WHERE odds_ratio >= 3 OR odds_ratio <= 0.333"):
        sig1[r[0]] = 'protective' if 'protective' in str(r[1]) else ('risk' if 'associated' in str(r[1]) else 'neutral')

    sig2 = {}
    for r in c.execute(f"SELECT te_id, direction FROM te_{s2}_association WHERE odds_ratio >= 3 OR odds_ratio <= 0.333"):
        sig2[r[0]] = 'protective' if 'protective' in str(r[1]) else ('risk' if 'associated' in str(r[1]) else 'neutral')

    # TEs significant in both
    both = set(sig1.keys()) & set(sig2.keys())

    concordant = 0
    conflicting = 0
    for te_id in both:
        if sig1[te_id] == sig2[te_id]:
            concordant += 1
        elif sig1[te_id] != 'neutral' and sig2[te_id] != 'neutral':
            conflicting += 1

    total_double = len(both)
    total_sig1 = len(sig1)
    total_sig2 = len(sig2)

    pair_key = f"{s1}-{s2}"
    concordance_matrix[pair_key] = {
        'n_sig_s1': total_sig1,
        'n_sig_s2': total_sig2,
        'n_both_sig': total_double,
        'concordant': concordant,
        'conflicting': conflicting,
        'concordant_pct': concordant / max(total_double, 1) * 100,
        'conflicting_pct': conflicting / max(total_double, 1) * 100,
        'overlap_pct': total_double / max(total_sig1, total_sig2, 1) * 100,
    }
    pair_details[pair_key] = {'sig1_only': list(sig1.keys())[:5], 'sig2_only': list(sig2.keys())[:5]}

# Print matrix
print(f"\n{'胁迫对':<16} {'胁迫1显著':<12} {'胁迫2显著':<12} {'共同显著':<10} {'一致%':<10} {'冲突%':<10} {'重叠率%':<10}")
print("-" * 80)
for pair_key in concordance_matrix:
    d = concordance_matrix[pair_key]
    s1, s2 = pair_key.split('-')
    print(f"{STRESS_LABELS[s1]}-{STRESS_LABELS[s2]}{'':<8} "
          f"{d['n_sig_s1']:<12} {d['n_sig_s2']:<12} {d['n_both_sig']:<10} "
          f"{d['concordant_pct']:<10.1f} {d['conflicting_pct']:<10.1f} {d['overlap_pct']:<10.1f}")

# Identify top concordant pairs (co-adaptation) and conflicting pairs (trade-off)
sorted_pairs = sorted(concordance_matrix.items(), key=lambda x: -x[1]['n_both_sig'])
print(f"\n★ 共同显著TE最多的胁迫对:")
for pk, d in sorted_pairs[:5]:
    s1, s2 = pk.split('-')
    print(f"  {STRESS_LABELS[s1]}-{STRESS_LABELS[s2]}: {d['n_both_sig']}个TE (一致{d['concordant_pct']:.0f}%, 冲突{d['conflicting_pct']:.0f}%)")

sorted_conflict = sorted(concordance_matrix.items(), key=lambda x: -x[1]['conflicting_pct'])
print(f"\n★ 冲突率最高的胁迫对 (权衡信号):")
for pk, d in sorted_conflict[:5]:
    if d['n_both_sig'] >= 5:
        s1, s2 = pk.split('-')
        print(f"  {STRESS_LABELS[s1]}-{STRESS_LABELS[s2]}: 冲突{d['conflicting_pct']:.0f}% ({d['conflicting']}/{d['n_both_sig']})")

# ═══════════════════════════════════════════════════════════════
# 2. NONA_BOKRA Genome-Wide Unique TE Systematic Identification
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*80}")
print(f"分析2: NONA_BOKRA全基因组独有TE系统鉴定")
print(f"{'='*80}")

n_idx = SAMPLES.index('N')

def gt_to_allele(gt):
    if gt in ('0/0','0|0','0'): return 0
    if gt in ('0/1','1/0','0|1','1|0'): return 1
    if gt in ('1/1','1|1','1'): return 2
    return None

nona_unique_tes = []
nona_het_tes = []
with open(VCF, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        chrom = parts[0]
        pos = int(parts[1])
        vcf_id = parts[2]
        gts = parts[9:] if len(parts) > 9 else []

        if len(gts) <= n_idx: continue

        n_gt = gts[n_idx]
        n_allele = gt_to_allele(n_gt)

        if n_allele is None: continue

        # Check other 11 varieties
        other_alleles = []
        for j, s in enumerate(SAMPLES):
            if j == n_idx or j >= len(gts): continue
            a = gt_to_allele(gts[j])
            if a is not None:
                other_alleles.append(a)

        if not other_alleles: continue

        is_nona_unique_alt = (n_allele > 0 and all(a == 0 for a in other_alleles))
        is_nona_het = (n_allele == 1 and all(a == 0 for a in other_alleles))
        is_nona_hom = (n_allele == 2 and all(a == 0 for a in other_alleles))
        is_nona_unique_ref = (n_allele == 0 and all(a > 0 for a in other_alleles))

        if is_nona_unique_alt:
            anno = [x.replace('ANNO=','') for x in parts[7].split(';') if x.startswith('ANNO=')]
            nona_unique_tes.append({
                'chrom': chrom, 'pos': pos, 'vcf_id': vcf_id,
                'n_genotype': 'het' if is_nona_het else 'hom_alt',
                'te_type': anno[0] if anno else 'NA',
                'gts': {s: gts[SAMPLES.index(s)] if SAMPLES.index(s) < len(gts) else './.' for s in SAMPLES},
            })
            if is_nona_het:
                nona_het_tes.append(vcf_id)

print(f"\nNONA_BOKRA独占替代等位TE: {len(nona_unique_tes)}个 (het={len(nona_het_tes)}, hom={len(nona_unique_tes)-len(nona_het_tes)})")
print(f"NONA_BOKRA独占RR(所有其他品种都有替代等位): 见下方分析")

# Chromosome distribution
chr_dist = Counter(t['chrom'] for t in nona_unique_tes)
print(f"\n染色体分布:")
for ch in sorted(chr_dist.keys(), key=lambda x: int(x.replace('Chr',''))):
    print(f"  {ch}: {chr_dist[ch]} TEs")

# Check stress associations for NONA unique TEs
print(f"\n胁迫关联分析 (NONA独占TE中的显著胁迫TE):")
nona_ids = set()
for t in nona_unique_tes:
    db_te = c.execute("SELECT id FROM te_insertions WHERE vcf_id=?", (t['vcf_id'],)).fetchone()
    if db_te:
        nona_ids.add(db_te[0])

nona_stress_counts = Counter()
nona_stress_details = defaultdict(list)
for te_id in nona_ids:
    for stress in STRESSES:
        r = c.execute(f"SELECT odds_ratio, direction, category FROM te_{stress}_association WHERE te_id=? AND (odds_ratio >= 3 OR odds_ratio <= 0.333)", (te_id,)).fetchone()
        if r:
            nona_stress_counts[stress] += 1
            nona_stress_details[stress].append({
                'te_id': te_id,
                'or_val': r[0],
                'direction': 'protective' if 'protective' in r[1] else 'risk'
            })

print(f"  胁迫显著TE总数: {sum(nona_stress_counts.values())}")
for stress in STRESSES:
    cnt = nona_stress_counts[stress]
    details = nona_stress_details[stress]
    n_prot = sum(1 for d in details if d['direction'] == 'protective')
    n_risk = sum(1 for d in details if d['direction'] == 'risk')
    if cnt > 0:
        print(f"  {STRESS_LABELS[stress]}: {cnt}个 (保护{n_prot}, 风险{n_risk})")

# TE type enrichment in NONA unique
nona_type_counts = Counter(t['te_type'] for t in nona_unique_tes)
print(f"\nNONA独占TE类型分布 (前10):")
for tt, cnt in nona_type_counts.most_common(10):
    print(f"  {tt}: {cnt}")

# Compare with genome-wide: NONA has more protective or risk TEs?
nona_protective = sum(1 for stress, details in nona_stress_details.items() for d in details if d['direction']=='protective')
nona_risk = sum(1 for stress, details in nona_stress_details.items() for d in details if d['direction']=='risk')
print(f"\nNONA独有TE的胁迫方向: 保护{nona_protective}, 风险{nona_risk}")

# ═══════════════════════════════════════════════════════════════
# 3. TE Superfamily Enrichment in Stress Associations
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*80}")
print(f"分析3: TE超家族胁迫关联富集分析")
print(f"{'='*80}")

# Genome-wide TE superfamily distribution (background)
all_sf = Counter()
for r in c.execute("SELECT te_superfamily FROM te_insertions"):
    sf = r[0] or 'Unclassified'
    # Simplify: take first component
    sf = sf.split(',')[0] if sf else 'Unclassified'
    all_sf[sf] += 1
total_all = sum(all_sf.values())

# For each stress, get superfamily distribution of significant TEs
enrichment_results = {}
for stress in STRESSES:
    sig_sf = Counter()
    for r in c.execute(f"""
        SELECT ti.te_superfamily
        FROM te_{stress}_association ta
        JOIN te_insertions ti ON ta.te_id=ti.id
        WHERE ta.odds_ratio >= 3 OR ta.odds_ratio <= 0.333
    """):
        sf = r[0] or 'Unclassified'
        sf = sf.split(',')[0] if sf else 'Unclassified'
        sig_sf[sf] += 1

    total_sig = sum(sig_sf.values())

    if total_sig == 0:
        print(f"\n{STRESS_LABELS[stress]}: 0个显著TE, 跳过")
        continue

    # Enrichment analysis for top superfamilies
    enrichments = []
    for sf in all_sf:
        obs = sig_sf.get(sf, 0)
        exp = all_sf[sf] / total_all * total_sig

        if obs >= 5 and exp >= 1:
            # 2x2 contingency table
            cell_a = obs
            cell_b = total_sig - obs
            cell_c = all_sf[sf] - obs
            cell_d = total_all - total_sig - cell_c

            if min(cell_a, cell_b, cell_c, cell_d) >= 0:
                or_val = (cell_a * cell_d) / max(cell_b * cell_c, 1)
                try:
                    _, p_val = fisher_exact([[cell_a, cell_b], [cell_c, cell_d]])
                except:
                    p_val = 1.0

                enrichments.append({
                    'superfamily': sf,
                    'observed': int(obs),
                    'expected': float(exp),
                    'fold_change': float(obs / max(exp, 0.01)),
                    'p_value': float(p_val),
                    'significant': bool(p_val < 0.05),
                })

    enrichments.sort(key=lambda x: -x['fold_change'])
    enrichment_results[stress] = enrichments

    print(f"\n{STRESS_LABELS[stress]}: {total_sig}显著TE")
    print(f"  {'超家族':<30} {'观察':<8} {'期望':<8} {'倍数':<8} {'p值':<10}")
    print(f"  {'-'*64}")

    for e in enrichments[:8]:
        sig_mark = '*' if e['significant'] else ''
        print(f"  {e['superfamily']:<30} {e['observed']:<8} {e['expected']:<8.1f} {e['fold_change']:<8.2f} {e['p_value']:<10.4f} {sig_mark}")

    # Also show depleted superfamilies
    enrichments.sort(key=lambda x: x['fold_change'])
    print(f"\n  贫化 (fold<1, p<0.05):")
    depleted = [e for e in enrichments if e['fold_change'] < 1 and e['significant']]
    if depleted:
        for e in depleted[:3]:
            print(f"    {e['superfamily']}: FC={e['fold_change']:.2f}, p={e['p_value']:.4f}")
    else:
        print(f"    无显著贫化")

# ═══════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════
output = {
    'cross_stress_concordance': {
        pk: {
            'n_sig_s1': d['n_sig_s1'],
            'n_sig_s2': d['n_sig_s2'],
            'n_both_sig': d['n_both_sig'],
            'concordant': d['concordant'],
            'conflicting': d['conflicting'],
            'concordant_pct': round(d['concordant_pct'], 1),
            'conflicting_pct': round(d['conflicting_pct'], 1),
            'overlap_pct': round(d['overlap_pct'], 1),
        } for pk, d in concordance_matrix.items()
    },
    'nona_bokra_unique_tes': {
        'total_count': len(nona_unique_tes),
        'het_count': len(nona_het_tes),
        'hom_alt_count': len(nona_unique_tes) - len(nona_het_tes),
        'chromosome_distribution': dict(chr_dist),
        'stress_associations': {
            stress: {
                'count': len(nona_stress_details[stress]),
                'protective': sum(1 for d in nona_stress_details[stress] if d['direction']=='protective'),
                'risk': sum(1 for d in nona_stress_details[stress] if d['direction']=='risk'),
            } for stress in STRESSES
        },
        'te_type_distribution': dict(nona_type_counts.most_common(20)),
        'top_20_unique_tes': nona_unique_tes[:20],
    },
    'te_superfamily_enrichment': {
        stress: [
            {'superfamily': e['superfamily'], 'fold_change': round(e['fold_change'], 2),
             'observed': e['observed'], 'expected': round(e['expected'], 1),
             'p_value': round(e['p_value'], 4), 'significant': e['significant']}
            for e in enrichments[:10]
        ] for stress, enrichments in enrichment_results.items() if enrichments
    },
}

out_path = os.path.join(OUT_DIR, 'hypothesis_support_analysis.json')
with open(out_path, 'w', encoding='utf-8') as f:
    json.dump(output, f, indent=2, ensure_ascii=False)

print(f"\n{'='*80}")
print(f"分析完成, 结果保存至: {out_path}")
conn.close()
