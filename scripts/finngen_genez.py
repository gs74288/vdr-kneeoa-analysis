# -*- coding: utf-8 -*-
"""Download FinnGen endpoints and build gene-level signed z-scores (same rule as GCST007090)."""
import sys, io, os, gzip, re
import numpy as np
import pandas as pd
from scipy import stats
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
GTF = os.path.join(D, 'Homo_sapiens.GRCh38.112.gtf.gz')
BASE = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/'
PAD = 10000
OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))

# ---------- GTF intervals ----------
print('parsing GTF ...', flush=True)
chrom_genes = {}
pat = re.compile(r'gene_id "([^".]+)(?:\.[0-9]+)?";.*?gene_name "([^"]+)";')
with gzip.open(GTF, 'rt', errors='replace') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'gene':
            continue
        m = pat.search(p[8])
        if not m:
            continue
        chrom_genes.setdefault(p[0], []).append((int(p[3]) - PAD, int(p[4]) + PAD, m.group(2)))
for c in chrom_genes:
    chrom_genes[c].sort()
    arr = np.array([(s, e) for s, e, _ in chrom_genes[c]], dtype=np.int64)
    chrom_genes[c] = (arr[:, 0], arr[:, 1], np.array([n for _, _, n in chrom_genes[c]]))
print('chromosomes:', len(chrom_genes), flush=True)


def build(ep):
    gz = os.path.join(D, 'finngen_R12_%s.gz' % ep)
    if not os.path.exists(gz):
        print('downloading', ep, flush=True)
        with OP.open(urllib.request.Request(BASE + 'finngen_R12_%s.gz' % ep,
                                            headers={'User-Agent': 'Mozilla/5.0'}),
                     timeout=3600) as r, open(gz, 'wb') as fh:
            while True:
                b = r.read(8 << 20)
                if not b:
                    break
                fh.write(b)
        print('  %.0f MB' % (os.path.getsize(gz) / 1e6), flush=True)

    best = {}   # gene -> (p, beta)
    cols = None
    n = 0
    with gzip.open(gz, 'rt', errors='replace') as fh:
        header = fh.readline().rstrip('\n').lstrip('#').split('\t')
        cols = {k: i for i, k in enumerate(header)}
        ci, pi, qi, bi = cols['chrom'], cols['pos'], cols['pval'], cols['beta']
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(header):
                continue
            try:
                c = f[ci]
                if c not in chrom_genes:
                    continue
                pos = int(f[pi]); pval = float(f[qi]); beta = float(f[bi])
            except Exception:
                continue
            if pval <= 0 or pval > 1:
                continue
            starts, ends, names = chrom_genes[c]
            j = np.searchsorted(starts, pos, side='right') - 1
            while j >= 0 and ends[j] >= pos:
                g = names[j]
                cur = best.get(g)
                if cur is None or pval < cur[0]:
                    best[g] = (pval, beta)
                j -= 1
                if j < 0:
                    break
                if starts[j] < pos - 300000:
                    break
            n += 1
    print('%s: %d SNPs scanned, %d genes' % (ep, n, len(best)), flush=True)

    rows = []
    for g, (p, b) in best.items():
        z = stats.norm.isf(p / 2)
        rows.append((g, z if b >= 0 else -z, p))
    df = pd.DataFrame(rows, columns=['GENE', 'Z', 'P'])
    out = os.path.join(D, 'finngen_%s_gene_z.tsv' % ep)
    df.to_csv(out, sep='\t', index=False)
    print('written', out, df.shape, flush=True)
    return out


for ep in ('T2D', 'M13_ARTHROSIS'):
    build(ep)
