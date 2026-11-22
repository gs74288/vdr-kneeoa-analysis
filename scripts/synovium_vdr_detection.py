# -*- coding: utf-8 -*-
"""GSE283080 inflamed synovium: VDR detection in macrophage-like cells."""
import sys, io, os, gzip, urllib.request
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {'User-Agent': 'Mozilla/5.0'}
D = r"C:\Users\liuha\OpenClawVDR\data\synovium\GSE283080"
os.makedirs(D, exist_ok=True)
SERIES = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8655nnn/'

SAMPLES = {'Infla1': 'GSM8655601_Infla1', 'Infla2': 'GSM8655602_Infla2'}


def fetch(name, gsm):
    p = os.path.join(D, name)
    if os.path.exists(p):
        return p
    url = '%s%s/suppl/%s' % (SERIES, gsm, name)
    with OP.open(urllib.request.Request(url, headers=UA), timeout=900) as r, open(p, 'wb') as fh:
        fh.write(r.read())
    print('  %-46s %.1f MB' % (name, os.path.getsize(p) / 1e6), flush=True)
    return p


for label, prefix in SAMPLES.items():
    gsm = prefix.split('_')[0]
    print('\n===', label, '===')
    for suffix in ('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz'):
        fetch('%s_%s' % (prefix, suffix), gsm)

    feat = pd.read_csv(os.path.join(D, '%s_features.tsv.gz' % prefix), sep='\t', header=None)
    bars = pd.read_csv(os.path.join(D, '%s_barcodes.tsv.gz' % prefix), sep='\t', header=None)
    with gzip.open(os.path.join(D, '%s_matrix.mtx.gz' % prefix), 'rb') as fh:
        M = scipy.io.mmread(fh).tocsc()
    print('  matrix:', M.shape, '| genes', len(feat), '| cells', len(bars))

    symbols = feat[1].astype(str).values if feat.shape[1] > 1 else feat[0].astype(str).values
    # dedupe symbols
    s = pd.Series(symbols)
    keep = ~s.duplicated()
    M = M[keep.values, :]
    symbols = s[keep.values].values
    idx = {g: i for i, g in enumerate(symbols)}

    tot = np.asarray(M.sum(0)).ravel()
    ngene = np.asarray((M > 0).sum(0)).ravel()
    ok = (ngene >= 200) & (tot > 0)
    print('  cells passing QC (>=200 genes): %d/%d' % (ok.sum(), len(ok)))

    def vec(g):
        return np.asarray(M[idx[g], :].todense()).ravel() if g in idx else np.zeros(M.shape[1])

    vdr = vec('VDR')
    cd68, lyz, cd14 = vec('CD68'), vec('LYZ'), vec('CD14')
    prg4, thy1, col1a1 = vec('PRG4'), vec('THY1'), vec('COL1A1')

    mac = (cd68 > 0) & (lyz > 0) & ok
    fib = ((prg4 > 0) | (thy1 > 0) | (col1a1 > 0)) & ~mac & ok
    other = ok & ~mac & ~fib
    print('  macrophage-like: %d | fibroblast-like: %d | other: %d' % (mac.sum(), fib.sum(), other.sum()))
    for name, m in (('macrophage-like', mac), ('fibroblast-like', fib), ('other', other)):
        if m.sum() >= 20:
            print('    %-16s VDR detected %5d/%5d (%.1f%%)'
                  % (name, int((vdr[m] > 0).sum()), int(m.sum()), 100 * (vdr[m] > 0).mean()))
    print('  overall VDR detection: %.1f%%' % (100 * (vdr[ok] > 0).mean()))
