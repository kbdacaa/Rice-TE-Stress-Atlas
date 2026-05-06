"""P5: OGR基因组基础模型 × TE多态性交叉验证实验"""
import torch, sqlite3, json, os, gzip, random, pickle
import numpy as np
from transformers import AutoModelForCausalLM, AutoTokenizer
from scipy.spatial.distance import cosine, euclidean
from scipy.stats import spearmanr, mannwhitneyu

MODEL_DIR = 'D:/project/AutoResearch-rice-T2T/onegenome_rice_weights'
DB = 'D:/project/AutoResearch-rice-T2T/data/rice_integration.db'
MH63_GENOME = 'D:/fanjiyinzu/RICEPTEDB/MH63RS3.genome.fasta.gz'
OUT = 'D:/project/AutoResearch-rice-T2T/data/ogr_te_results.json'

# ── Step 1: Load Chr04 sequence ──
print('Step 1: Loading MH63 Chr04...')
cache = 'D:/project/AutoResearch-rice-T2T/onegenome_rice_weights/chr04_seq.pkl'
if os.path.exists(cache):
    with open(cache, 'rb') as f:
        chr04_seq = pickle.load(f)
else:
    with gzip.open(MH63_GENOME, 'rb') as f:
        raw = f.read()
    text = raw.decode('ascii', errors='ignore')
    seq = ''; in_chr04 = False
    for line in text.split('\n'):
        if line.startswith('>') and '04' in line:
            in_chr04 = True; continue
        elif line.startswith('>') and in_chr04: break
        elif in_chr04: seq += line.strip()
    chr04_seq = seq
    with open(cache, 'wb') as f:
        pickle.dump(seq, f)
print(f'  Chr04: {len(chr04_seq):,}bp')

# ── Step 2: Load OGR model ──
print('Step 2: Loading OGR model...')
tokenizer = AutoTokenizer.from_pretrained(MODEL_DIR, trust_remote_code=True)
model = AutoModelForCausalLM.from_pretrained(
    MODEL_DIR, torch_dtype=torch.float32, device_map='cuda:0',
    trust_remote_code=True, low_cpu_mem_usage=True)
model.eval()
print('  Model loaded')

def get_embedding(seq):
    if len(seq) < 100: return None
    inputs = tokenizer(seq[:4000], return_tensors='pt', truncation=True, max_length=1024).to('cuda:0')
    with torch.no_grad():
        out = model(**inputs, output_hidden_states=True)
        return out.hidden_states[-1].mean(dim=1).squeeze().cpu().numpy()

# ── Step 3: Extract TE flanking sequences ──
print('Step 3: Extracting TE flanking sequences...')
conn = sqlite3.connect(DB); c = conn.cursor()

# Get cold-significant TEs on Chr04
cold_tes = []
for r in c.execute("""
    SELECT ta.te_id, ti.position, ta.odds_ratio, ta.direction, ta.category
    FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id
    WHERE ti.chromosome='Chr04' AND ti.position BETWEEN 1100000 AND 5500000
    ORDER BY ti.position
"""):
    cold_tes.append({'te_id': r[0], 'pos': r[1], 'or_val': r[2], 'direction': r[3], 'category': r[4]})

# Also get neutral TEs (category='none') as controls
neutral_tes = []
for r in c.execute("""
    SELECT ta.te_id, ti.position, ta.odds_ratio
    FROM te_cold_association ta JOIN te_insertions ti ON ta.te_id=ti.id
    WHERE ti.chromosome='Chr04' AND ti.position BETWEEN 1100000 AND 5500000
    AND ta.category='none' AND ta.odds_ratio = 1.0
    ORDER BY RANDOM() LIMIT 30
"""):
    neutral_tes.append({'te_id': r[0], 'pos': r[1], 'or_val': r[2], 'direction': 'neutral', 'category': 'none'})

conn.close()

# Select TEs for experiment
test_tes = []
# All cold-significant
for t in cold_tes:
    if t['category'] == 'strong':
        test_tes.append(t)
# Random neutral
for t in neutral_tes[:15]:
    test_tes.append(t)
# Add the key diagnostic TEs
key_positions = [1334699, 1336356, 1341230, 1406577, 1406820, 1485066]  # TE#8261-8279
for pos in key_positions:
    found = [t for t in cold_tes if abs(t['pos'] - pos) < 100]
    if found and found[0] not in test_tes:
        test_tes.append(found[0])

test_tes = list({t['te_id']: t for t in test_tes}.values())  # deduplicate
ns = sum(1 for t in test_tes if t['category']=='strong')
nn = sum(1 for t in test_tes if t['category']=='none')
print(f'  Test set: {len(test_tes)} TEs ({ns} strong, {nn} neutral)')

# ── Step 4: Extract sequences and get embeddings ──
print('Step 4: Running OGR inference...')
FLANK = 500
results = []
for i, te in enumerate(test_tes):
    pos = te['pos']
    start = max(0, pos - FLANK)
    end = min(len(chr04_seq), pos + FLANK)
    seq = chr04_seq[start:end]
    if len(seq) < 400: continue

    emb = get_embedding(seq)
    if emb is not None:
        results.append({
            'te_id': te['te_id'], 'pos': pos,
            'or_val': te['or_val'], 'direction': te['direction'], 'category': te['category'],
            'seq_len': len(seq),
            'embedding': emb.tolist()
        })

    if (i+1) % 10 == 0:
        print(f'  Progress: {i+1}/{len(test_tes)}')

print(f'  Completed: {len(results)} embeddings')

# ── Step 5: Analysis ──
print('\n' + '='*60)
print('Step 5: Cross-Validation Analysis')
print('='*60)

strong_embs = [r['embedding'] for r in results if r['category'] == 'strong']
neutral_embs = [r['embedding'] for r in results if r['category'] == 'none']
prot_embs = [r['embedding'] for r in results if 'protective' in str(r['direction'])]
risk_embs = [r['embedding'] for r in results if 'associated' in str(r['direction'])]

# Test 1: Can OGR distinguish cold-significant from neutral TE regions?
if strong_embs and neutral_embs:
    # Pairwise cosine distances within each group
    strong_dists = []
    for i in range(min(len(strong_embs), 10)):
        for j in range(i+1, min(len(strong_embs), 10)):
            strong_dists.append(cosine(strong_embs[i], strong_embs[j]))
    neutral_dists = []
    for i in range(min(len(neutral_embs), 10)):
        for j in range(i+1, min(len(neutral_embs), 10)):
            neutral_dists.append(cosine(neutral_embs[i], neutral_embs[j]))

    # Cross-group distances
    cross_dists = []
    for se in strong_embs[:10]:
        for ne in neutral_embs[:10]:
            cross_dists.append(cosine(se, ne))

    print(f'\nTest 1: Embedding Distance Comparison')
    print(f'  Strong TEs (within-group): {np.mean(strong_dists):.4f} ± {np.std(strong_dists):.4f}')
    print(f'  Neutral TEs (within-group): {np.mean(neutral_dists):.4f} ± {np.std(neutral_dists):.4f}')
    print(f'  Cross-group (strong↔neutral): {np.mean(cross_dists):.4f} ± {np.std(cross_dists):.4f}')

    # Mann-Whitney test: are cross-group distances larger than within-group?
    try:
        u1, p1 = mannwhitneyu(cross_dists, strong_dists, alternative='greater')
        print(f'  Cross vs Strong-within: MW U={u1:.0f}, p={p1:.4f}')
    except: pass

# Test 2: Does OGR embedding distance correlate with GWAS OR?
ors = [r['or_val'] for r in results if r['category'] == 'strong']
or_embeddings = [r['embedding'] for r in results if r['category'] == 'strong']
if len(ors) >= 10:
    # Compute pairwise distances
    dists = []
    or_diffs = []
    for i in range(len(or_embeddings)):
        for j in range(i+1, len(or_embeddings)):
            dists.append(cosine(or_embeddings[i], or_embeddings[j]))
            or_diffs.append(abs(np.log2(max(ors[i], 0.01)) - np.log2(max(ors[j], 0.01))))
    rho, p_rho = spearmanr(dists, or_diffs)
    print(f'\nTest 2: OGR Distance vs GWAS OR Difference')
    print(f'  Spearman ρ = {rho:.4f}, p = {p_rho:.4f}')

# Test 3: Protective vs Risk TE embedding comparison
if prot_embs and risk_embs:
    prot_centroid = np.mean(prot_embs, axis=0)
    risk_centroid = np.mean(risk_embs, axis=0)
    centroid_dist = cosine(prot_centroid, risk_centroid)

    # Within-group consistency
    prot_dists = []
    for i in range(min(len(prot_embs), 8)):
        for j in range(i+1, min(len(prot_embs), 8)):
            prot_dists.append(cosine(prot_embs[i], prot_embs[j]))
    risk_dists = []
    for i in range(min(len(risk_embs), 8)):
        for j in range(i+1, min(len(risk_embs), 8)):
            risk_dists.append(cosine(risk_embs[i], risk_embs[j]))

    print(f'\nTest 3: Protective vs Risk TE Embeddings')
    print(f'  Protective within-group: {np.mean(prot_dists):.4f} ± {np.std(prot_dists):.4f}')
    print(f'  Risk within-group: {np.mean(risk_dists):.4f} ± {np.std(risk_dists):.4f}')
    print(f'  Protective↔Risk centroid distance: {centroid_dist:.4f}')
    sep = 'OGR separates protective vs risk' if centroid_dist > np.mean(prot_dists)*1.2 else 'No clear separation'
    print(f'  -> {sep}')

# Save
output = {
    'n_tes': len(results),
    'test1': {
        'strong_within_mean': float(np.mean(strong_dists)) if strong_dists else None,
        'neutral_within_mean': float(np.mean(neutral_dists)) if neutral_dists else None,
        'cross_group_mean': float(np.mean(cross_dists)) if cross_dists else None,
    },
    'test2_spearman': {'rho': float(rho) if ors else None, 'p': float(p_rho) if ors else None},
    'test3': {
        'protective_within': float(np.mean(prot_dists)) if prot_dists else None,
        'risk_within': float(np.mean(risk_dists)) if risk_dists else None,
        'centroid_distance': float(centroid_dist) if prot_embs and risk_embs else None,
    }
}
with open(OUT, 'w') as f:
    json.dump(output, f, indent=2)
print(f'\nResults saved to: {OUT}')
