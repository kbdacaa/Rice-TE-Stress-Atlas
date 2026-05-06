"""OGR模型验证：亚群嵌入 + TE效应 + 文献闭环"""
import torch, gzip, json, os, pickle, sys
from transformers import AutoModelForCausalLM, AutoTokenizer
from scipy.spatial.distance import cosine
import numpy as np

CACHE = 'D:/project/AutoResearch-rice-T2T/onegenome_rice_weights/chr04_seqs.pkl'
MODEL_DIR = 'D:/project/AutoResearch-rice-T2T/onegenome_rice_weights'

def load_chr04_cached(path, label):
    """Load Chr04 from genome FASTA with caching"""
    if os.path.exists(CACHE):
        with open(CACHE, 'rb') as f:
            data = pickle.load(f)
        if label in data:
            print(f'  {label}: loaded from cache ({len(data[label]):,}bp)', flush=True)
            return data[label]
    return None

def load_chr04(path, label):
    """Load Chr04 from genome FASTA (fast: read all at once)"""
    result = load_chr04_cached(path, label)
    if result: return result

    print(f'  Loading {label} from FASTA (decompressing ~120MB)...', flush=True)
    with gzip.open(path, 'rb') as f:
        raw = f.read()
    text = raw.decode('ascii', errors='ignore')
    chroms = {}
    cur = None
    for line in text.split('\n'):
        if line.startswith('>'):
            cur = line.strip().split()[0][1:]
            chroms[cur] = ''
        elif cur:
            chroms[cur] += line.strip()
    for k in chroms:
        if '04' in k or 'Chr04' in k:
            print(f'  {label}: found Chr04 ({k}) {len(chroms[k]):,}bp', flush=True)
            return chroms[k]
    return ''

# ── Step 1: Cache genome sequences ──
print('='*60, flush=True)
print('Step 1: 提取Chr04序列...', flush=True)
genomes = {
    'MH63': 'D:/fanjiyinzu/RICEPTEDB/MH63RS3.genome.fasta.gz',
    'NIP': 'D:/fanjiyinzu/RICEPTEDB/R.genome.fasta.gz',
    'NONA': 'D:/fanjiyinzu/RICEPTEDB/N.genome.fasta.gz',
}

# Try cache first, fall back to loading
need_save = not os.path.exists(CACHE)
chr04_seqs = {}
if os.path.exists(CACHE):
    with open(CACHE, 'rb') as f:
        chr04_seqs = pickle.load(f)
    print(f'  Loaded from cache: {list(chr04_seqs.keys())}', flush=True)

for label, path in genomes.items():
    if label not in chr04_seqs:
        chr04_seqs[label] = load_chr04(path, label)
        need_save = True

if need_save:
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, 'wb') as f:
        pickle.dump(chr04_seqs, f)
    print(f'  Saved cache', flush=True)

for label, seq in chr04_seqs.items():
    print(f'  {label}: {len(seq):,}bp', flush=True)

# ── Step 2: Load OGR model ──
print(f'\nStep 2: 加载OGR模型...', flush=True)
tokenizer = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModelForCausalLM.from_pretrained(
    MODEL_DIR, torch_dtype=torch.bfloat16, device_map='cuda:0',
    trust_remote_code=True, low_cpu_mem_usage=True)
model.eval()
print(f'  Model on {next(model.parameters()).device}', flush=True)

def get_embedding(seq):
    inputs = tokenizer(seq[:4000], return_tensors='pt', truncation=True, max_length=2048).to('cuda:0')
    with torch.no_grad():
        out = model(**inputs, output_hidden_states=True)
        return out.hidden_states[-1].mean(dim=1).squeeze().cpu().numpy()

# ── Step 3: Test 1 - 亚群嵌入距离 ──
print(f'\n{"="*60}', flush=True)
print('Test 1: 亚群嵌入距离 (CTB4a/CTB2/DGK7/Random)', flush=True)
print('='*60, flush=True)

regions = [
    ('CTB4a_gene',    1332869, 1342717),
    ('CTB2_gene',      175000,  185000),
    ('DGK7_region', 29000000, 29010000),
    ('Random_ctrl_1', 5000000, 5010000),
    ('Random_ctrl_2',15000000,15010000),
]

results_t1 = []
for rname, s, e in regions:
    mh63_s = chr04_seqs['MH63'][s:min(e, len(chr04_seqs['MH63']))]
    nip_s = chr04_seqs['NIP'][s:min(e, len(chr04_seqs['NIP']))]
    nona_s = chr04_seqs['NONA'][s:min(e, len(chr04_seqs['NONA']))]

    if len(mh63_s) < 100: continue

    e_m = get_embedding(mh63_s)
    e_nip = get_embedding(nip_s)
    e_nona = get_embedding(nona_s)

    d_mn = cosine(e_m, e_nip)
    d_mnona = cosine(e_m, e_nona)
    results_t1.append({'region': rname, 'MH63_NIP': float(d_mn), 'MH63_NONA': float(d_mnona), 'bp': len(mh63_s)})
    print(f'  {rname:18s} MH63↔NIP={d_mn:.4f}  MH63↔NONA={d_mnona:.4f}  ({len(mh63_s)}bp)', flush=True)

# ── Step 4: Test 2 - TE#8262 效应 ──
print(f'\n{"="*60}', flush=True)
print('Test 2: TE#8262对OGR嵌入的影响', flush=True)
print('='*60, flush=True)

# CTB4a positon: MH63 Chr04 ~1,334,869-1,340,717
# TE#8262 is at 1,336,356 (ALLEL insertion inside CTB4a intron)
# MH63 = no TE (RR), NIP = has TE (AA)
pos = 1336356
flank = 1000

mh63_te_region = chr04_seqs['MH63'][pos-flank:pos+flank]
nip_te_region = chr04_seqs['NIP'][pos-flank:pos+flank]
nona_te_region = chr04_seqs['NONA'][pos-flank:pos+flank]

e_mh63_te = get_embedding(mh63_te_region)
e_nip_te = get_embedding(nip_te_region)
e_nona_te = get_embedding(nona_te_region)

d_te_mn = cosine(e_mh63_te, e_nip_te)
d_te_mnona = cosine(e_mh63_te, e_nona_te)
d_te_nipnona = cosine(e_nip_te, e_nona_te)

# Random control
rand_s = 8000000
mh63_rand = chr04_seqs['MH63'][rand_s:rand_s+2000]
nip_rand = chr04_seqs['NIP'][rand_s:rand_s+2000]
d_rand = cosine(get_embedding(mh63_rand), get_embedding(nip_rand))

print(f'  CTB4a TE#8262区域 (±{flank}bp):', flush=True)
print(f'    MH63(RR,noTE) ↔ NIP(AA,hasTE): {d_te_mn:.4f}', flush=True)
print(f'    MH63(RR,noTE) ↔ NONA(RR,noTE): {d_te_mnona:.4f}', flush=True)
print(f'    NIP(AA,hasTE) ↔ NONA(RR,noTE): {d_te_nipnona:.4f}', flush=True)
print(f'  Random control (2000bp): {d_rand:.4f}', flush=True)
print(f'  TE± embedding ratio: {d_te_mn/d_rand:.2f}x', flush=True)
signal = 'OGR能检测TE+-差异' if d_te_mn > 1.5 * d_rand else 'TE+-差异不显著'
print(f'  -> {signal}', flush=True)

# ── Step 5: Test 3 - 文献闭环 ──
print(f'\n{"="*60}', flush=True)
print('Test 3: OGR亚群分类 → 文献验证闭环', flush=True)
print('='*60, flush=True)

# Use CTB4a region to classify subspecies
# Indica anchor: MH63 random region (Chr04:2Mb)
# Japonica anchor: NIP random region (Chr04:2Mb)
ind_anchor_emb = get_embedding(chr04_seqs['MH63'][2000000:2005000])
jap_anchor_emb = get_embedding(chr04_seqs['NIP'][2000000:2005000])

# Classify each variety's CTB4a region
ctb4a_s, ctb4a_e = 1332869, 1342717
classify_results = {}
for label in ['MH63', 'NIP', 'NONA']:
    seq = chr04_seqs[label][ctb4a_s:min(ctb4a_e, len(chr04_seqs[label]))]
    emb = get_embedding(seq)
    d_ind = cosine(emb, ind_anchor_emb)
    d_jap = cosine(emb, jap_anchor_emb)
    pred = 'INDICA' if d_ind < d_jap else 'JAPONICA'
    classify_results[label] = {'d_ind': float(d_ind), 'd_jap': float(d_jap), 'pred': pred}

# Literature ground truth from CNKI matrix + breeding records
literature_truth = {
    'MH63': {'subspecies': 'indica', 'known_traits': '明恢63, indica II, 汕优63中抗稻瘟病/纹枯病/耐盐碱/耐高温[吴方喜2011]'},
    'NIP': {'subspecies': 'japonica', 'known_traits': '日本晴, temperate japonica, 冷耐受1级[RICEPTEDB], COLD1克隆品种[Ma2015]'},
    'NONA': {'subspecies': 'aus', 'known_traits': 'NONA_BOKRA, aus/cA, Sub1供体[Xu2006], 728独占冷TE, 冷适应特化'},
}

hdr = '  {:8s} {:8s} {:8s} {:10s} {:10s} {}'.format('品种','d_ind','d_jap','OGR预测','文献记录','一致性')
print(hdr, flush=True)
for label in ['MH63', 'NIP', 'NONA']:
    r = classify_results[label]
    lt = literature_truth[label]
    match = '✓' if (r['pred'] == 'INDICA' and lt['subspecies'] in ['indica','aus']) or \
                   (r['pred'] == 'JAPONICA' and lt['subspecies'] == 'japonica') else '?'
    di = r['d_ind']; dj = r['d_jap']; pr = r['pred']; ss = lt['subspecies']
    print(f'  {label:<8} {di:<8.4f} {dj:<8.4f} {pr:<10} {ss:<10} {match}', flush=True)

print(f'\n文献闭环总结:', flush=True)
correct = sum(1 for l in ['MH63','NIP','NONA']
              if (classify_results[l]['pred']=='INDICA' and literature_truth[l]['subspecies'] in ['indica','aus'])
              or (classify_results[l]['pred']=='JAPONICA' and literature_truth[l]['subspecies']=='japonica'))
print(f'  OGR分类正确率: {correct}/3', flush=True)
print(f'  文献验证来源: CNKI品种×性状矩阵 + 吴方喜2011 + 白和盛2001 + RICEPTEDB冷等级', flush=True)

# ── Save ──
output = {
    'test1_subspecies_embedding': results_t1,
    'test2_te8262_effect': {
        'MH63_NIP_dist': float(d_te_mn), 'MH63_NONA_dist': float(d_te_mnona),
        'NIP_NONA_dist': float(d_te_nipnona), 'random_control': float(d_rand),
        'ratio': float(d_te_mn/d_rand)
    },
    'test3_literature_loop': {
        'classification': classify_results,
        'ground_truth': literature_truth,
        'accuracy': f'{correct}/3'
    }
}
out_path = 'D:/project/AutoResearch-rice-T2T/data/ogr_validation_results.json'
with open(out_path, 'w', encoding='utf-8') as f:
    json.dump(output, f, indent=2, ensure_ascii=False)
print(f'\n结果已保存: {out_path}', flush=True)
print('DONE', flush=True)
