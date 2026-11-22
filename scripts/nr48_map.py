# -*- coding: utf-8 -*-
"""#1: cis-MR screen across the 48 human nuclear receptors against knee OA."""
import sys, io, os, json, gzip, time, re
import numpy as np
import pandas as pd
import urllib.request, urllib.error

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

NRS = ['NR0B1', 'NR0B2', 'NR1D1', 'NR1D2', 'NR1H2', 'NR1H3', 'NR1H4', 'NR1I2', 'NR1I3',
       'NR2C1', 'NR2C2', 'NR2E1', 'NR2E3', 'NR2F1', 'NR2F2', 'NR2F6', 'NR3C1', 'NR3C2',
       'NR4A1', 'NR4A2', 'NR4A3', 'NR5A1', 'NR5A2', 'NR6A1',
       'THRA', 'THRB', 'RARA', 'RARB', 'RARG', 'PPARA', 'PPARD', 'PPARG',
       'RORA', 'RORB', 'RORC', 'RXRA', 'RXRB', 'RXRG', 'ESR1', 'ESR2', 'ESRRA', 'ESRRB',
       'ESRRG', 'AR', 'PGR', 'VDR', 'HNF4A', 'HNF4G']
print('NR list:', len(NRS))

# ---- GTF lookup ----
GTF = os.path.join(D, 'Homo_sapiens.GRCh38.112.gtf.gz')
want = set(NRS)
found = {}
with gzip.open(GTF, 'rt', errors='replace') as fh:
    for ln in fh:
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene':
            continue
        attr = f[8]
        m = re.search(r'gene_name "([^"]+)"', attr)
        if not m or m.group(1) not in want:
            continue
        g = re.search(r'gene_id "([^"]+)"', attr)
        if not g:
            continue
        found[m.group(1)] = (g.group(1), f[0], int(f[3]), int(f[4]), f[6])
print('mapped:', len(found))
missing = [g for g in NRS if g not in found]
print('missing:', missing)

meta = pd.DataFrame([(k, v[0], v[1], v[2] if v[4] == '+' else v[3], v[4]) for k, v in found.items()],
                    columns=['symbol', 'ensg', 'chr', 'tss', 'strand'])
meta.to_csv(os.path.join(T, 'nr48_genes.csv'), index=False)
print(meta.head(4).to_string(index=False))
