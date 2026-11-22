# -*- coding: utf-8 -*-
"""Download GSE152805 synovium (3 donors) and compute VDR detection by cell class."""
import sys, io, os, gzip, urllib.request
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {'User-Agent': 'Mozilla/5.0'}
D = r"C:\Users\liuha\OpenClawVDR\data\GSE152805_syn"
os.makedirs(D, exist_ok=True)

SAMPLES = {'113': 'GSM4626763_SY_113', '116': 'GSM4626764_SY_116', '118': 'GSM4626765_SY_118'}
FORMS = {'GSM4626763': 'GSM4626nnn', 'GSM4626764': 'GSM4626nnn', 'GSM4626765': 'GSM4626nnn'}

rows = []
for donor, prefix in SAMPLES.items():
    gsm = prefix.split('_')[0]
    print('\n=== donor %s (%s) ===' % (donor, gsm))
    for suffix in ('barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz'):
        name = '%s.%s' % (prefix, suffix)
        p = os.path.join(D, name)
        if not os.path.exists(p):
            url = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4626nnn/%s/suppl/%s' % (gsm, name)
            with OP.open(urllib.request.Request(url, headers=UA), timeout=1200) as r, open(p, 'wb') as fh:
                fh.write(r.read())
        print('  %-34s %6.1f MB' % (name, os.path.getsize(p) / 1e6), flush=True)

    feat = pd.read_csv(os.path.join(D, '%s.features.tsv.gz' % prefix), sep='\t', header=None) \
        if os.path.exists(os.path.join(D, '%s.features.tsv.gz' % prefix)) else None
    if feat is None:
        feat = pd.read_csv(os.path.join(D, '%s.genes.tsv.gz' % prefix), sep='\t', header=None)
    with gzip.open(os.path.join(D, '%s.matrix.mtx.gz' % prefix), 'rb') as fh:
        M = scipy.io.mmread(fh).tocsc()
    print('  matrix:', M.shape)

    sym = feat[1].astype(str).values if feat.shape[1] > 1 else feat[0].astype(str).values
    s = pd.Series(sym)
    keep = ~s.duplicated()
    M = M[keep.values, :]
    sym = s[keep.values].values
    idx = {g: i for i, g in enumerate(sym)}

    tot = np.asarray(M.sum(0)).ravel()
    ng = np.asarray((M > 0).sum(0)).ravel()
    ok = (ng >= 200) & (tot > 0)
    print('  cells: %d | passing QC: %d' % (M.shape[1], ok.sum()))

    def vec(g):
        return np.asarray(M[idx[g], :].todense()).ravel() if g in idx else np.zeros(M.shape[1])

    vdr = vec('VDR')
    mac = (vec('CD68') > 0) & (vec('LYZ') > 0)
    fib = ((vec('PRG4') > 0) | (vec('THY1') > 0) | (vec('DCN') > 0)) & ~mac
    tcell = ((vec('CD3D') > 0) | (vec('CD3E') > 0)) & ~mac & ~fib
    endo = ((vec('PECAM1') > 0) | (vec('VWF') > 0)) & ~mac & ~fib & ~tcell
    classes = {'macrophage': mac, 'fibroblast': fib, 'T/NK': tcell, 'endothelial': endo}
    assigned = np.zeros(M.shape[1], dtype=bool)
    for m in classes.values():
        assigned |= m
    classes['other'] = ok & ~assigned

    for name, m in classes.items():
        mm = m & ok
        if mm.sum() >= 20:
            det = int((vdr[mm] > 0).sum())
            rows.append(dict(donor=donor, cell_class=name, n=int(mm.sum()),
                             detected=det, pct=round(100 * det / mm.sum(), 1)))
            print('    %-12s n=%-6d VDR+ %5d (%.1f%%)' % (name, mm.sum(), det, 100 * det / mm.sum()))

res = pd.DataFrame(rows)
res.to_csv(r"C:\Users\liuha\OpenClawVDR\repo\tables\GSE152805_synovium_vdr.csv", index=False)
print('\n', res.to_string(index=False))
