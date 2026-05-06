"""OGR滑动窗口嵌入分析: Chr04冷TE热点区域的品种间嵌入差异"""
import torch, pickle, json, os, sqlite3, gzip
import numpy as np
from transformers import AutoModelForCausalLM, AutoTokenizer
from scipy.spatial.distance import cosine

MODEL_DIR = 'D:/project/AutoResearch-rice-T2T/onegenome_rice_weights'
DB = 'D:/project/AutoResearch-rice-T2T/data/rice_integration.db'
CACHE_DIR = 'D:/project/AutoResearch-rice-T2T/onegenome_rice_weights'

# Step 1: Load Chr04 sequences
print('Loading Chr04 sequences...')
seqs = {}
for label, cache_name in [('MH63','chr04_seq.pkl'),('NIP','nip_chr04.pkl'),('NONA','nona_chr04.pkl')]:
    cache = os.path.join(CACHE_DIR, cache_name)
    if os.path.exists(cache):
        with open(cache, 'rb') as f:
            seqs[label] = pickle.load(f)
    else:
        path = {'MH63':'MH63RS3','NIP':'R','NONA':'N'}[label]
        with gzip.open(f'D:/fanjiyinzu/RICEPTEDB/{path}.genome.fasta.gz', 'rb') as f:
            raw = f.read()
        text = raw.decode('ascii', errors='ignore')
        parts = []; in_chr = False
        for line in text.split('\n'):
            if line.startswith('>') and '04' in line: in_chr = True; continue
            elif line.startswith('>') and in_chr: break
            elif in_chr: parts.append(line.strip())
        seqs[label] = ''.join(parts)
        with open(cache, 'wb') as f:
            pickle.dump(seqs[label], f)
    print(f'  {label}: {len(seqs[label]):,}bp')

# Step 2: Load OGR
print('Loading OGR model...')
tokenizer = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModelForCausalLM.from_pretrained(
    MODEL_DIR, torch_dtype=torch.float32, device_map='cuda:0',
    trust_remote_code=True, low_cpu_mem_usage=True)
model.eval()

def emb(seq):
    inputs = tokenizer(seq[:4000], return_tensors='pt', truncation=True, max_length=2048).to('cuda:0')
    with torch.no_grad():
        out = model(**inputs, output_hidden_states=True)
        return out.hidden_states[-1].mean(dim=1).squeeze().cpu().numpy()

# Step 3: Sliding window analysis on Chr04 cold hotspot
print('\nSliding window analysis (20kb window, 200kb step)...')
WINDOW = 20000; STEP = 200000
region_start = 1100000; region_end = 5500000

results = []
min_len = min(len(seqs[s]) for s in seqs)
effective_end = min(region_end, min_len - WINDOW)

for pos in range(region_start, effective_end, STEP):
    mh63_s = seqs['MH63'][pos:pos+WINDOW]
    nip_s = seqs['NIP'][pos:pos+WINDOW]
    nona_s = seqs['NONA'][pos:pos+WINDOW]
    if len(mh63_s) < 19000: continue

    e_m = emb(mh63_s); e_n = emb(nip_s); e_a = emb(nona_s)
    d_mn = cosine(e_m, e_n); d_ma = cosine(e_m, e_a); d_na = cosine(e_n, e_a)

    # Count TEs in this window
    conn = sqlite3.connect(DB); c = conn.cursor()
    te_count = c.execute(
        "SELECT COUNT(*), SUM(CASE WHEN ta.category='strong' THEN 1 ELSE 0 END) FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id WHERE ti.chromosome='Chr04' AND ti.position BETWEEN ? AND ?",
        (pos, pos+WINDOW)).fetchone()
    conn.close()

    results.append({
        'pos': pos, 'pos_mb': pos/1e6,
        'd_MH63_NIP': float(d_mn), 'd_MH63_NONA': float(d_ma), 'd_NIP_NONA': float(d_na),
        'd_max': float(max(d_mn, d_ma, d_na)),
        'te_count': te_count[0] if te_count else 0,
        'te_strong': te_count[1] if te_count else 0,
    })
    print(f'  {pos/1e6:.1f}Mb: MN={d_mn:.4f} MA={d_ma:.4f} NA={d_na:.4f} TEs={te_count[0]} strong={te_count[1]}')

# Step 4: Fine-grained analysis of top divergent regions
print('\nFine-grained analysis of top regions...')
results.sort(key=lambda x: -x['d_max'])
top_regions = results[:3]

fine_results = []
for r in top_regions:
    center = r['pos']
    for offset in range(-50000, 50001, 5000):
        pos = center + offset
        if pos < 0 or pos + WINDOW > min_len: continue
        mh63_s = seqs['MH63'][pos:pos+WINDOW]
        nip_s = seqs['NIP'][pos:pos+WINDOW]
        nona_s = seqs['NONA'][pos:pos+WINDOW]
        if len(mh63_s) < 19000: continue
        e_m = emb(mh63_s); e_n = emb(nip_s); e_a = emb(nona_s)
        fine_results.append({
            'pos': pos, 'center': center,
            'd_MH63_NIP': float(cosine(e_m, e_n)),
            'd_MH63_NONA': float(cosine(e_m, e_a)),
            'd_NIP_NONA': float(cosine(e_n, e_a)),
        })
    print(f'  Center {center/1e6:.1f}Mb: {len(fine_results)//3} fine windows done')

# Step 5: Save
output = {
    'sliding_window': results,
    'fine_grained': fine_results,
    'n_windows': len(results),
    'n_fine': len(fine_results),
}
out_path = 'D:/project/AutoResearch-rice-T2T/data/ogr_sliding_window.json'
with open(out_path, 'w') as f:
    json.dump(output, f, indent=2)
print(f'\nSaved {len(results)} coarse + {len(fine_results)} fine windows to {out_path}')
print('DONE')
