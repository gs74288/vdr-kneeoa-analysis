# -*- coding: utf-8 -*-
"""H3K27ac peak overlap for the 12q13 signal; also locate java.exe after winget install."""
import sys, io, os, gzip, glob
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
CH = os.path.join(D, 'chromatin')

print('=== locate java.exe ===')
cands = []
for pat in (r"C:\Program Files\Eclipse Adoptium\*\bin\java.exe",
            r"C:\Program Files\Java\*\bin\java.exe",
            r"C:\Program Files (x86)\Eclipse Adoptium\*\bin\java.exe",
            os.path.expandvars(r"%LOCALAPPDATA%\Programs\Eclipse Adoptium\*\bin\java.exe")):
    cands += glob.glob(pat)
print('  found:', cands if cands else 'none')

print('\n=== peaks ===')
bed = pd.read_csv(os.path.join(CH, 'GSE188124_H3K27ac_peaks_GRCh38.bed.gz'), sep='\t',
                  header=None, comment='#')
print('  columns:', bed.shape, '| first row:', bed.iloc[0].tolist()[:6])
bed.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'signal', 'p', 'q', 'peak'][:bed.shape[1]]
bed['chr'] = bed['chr'].astype(str)
print('  chroms:', bed.chr.value_counts().head(5).to_dict())
p12 = bed[bed.chr.isin(['chr12', '12'])].copy()
p12['chr'] = '12'
p12 = p12.sort_values('start').reset_index(drop=True)
print('  chr12 peaks:', len(p12), '| total bp: %.1f Mb' % ((p12.end - p12.start).sum() / 1e6))

# merge overlapping peaks
iv = []
for s, e in zip(p12.start.values, p12.end.values):
    if iv and s <= iv[-1][1]:
        iv[-1][1] = max(iv[-1][1], e)
    else:
        iv.append([s, e])
print('  merged intervals:', len(iv))


def hit(pos):
    lo, hi = 0, len(iv) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        if pos < iv[m][0]:
            hi = m - 1
        elif pos > iv[m][1]:
            lo = m + 1
        else:
            return True
    return False


print('\n=== GWAS file ===')
G = os.path.join(D, 'GCST90566800_chr12_46.8-48.9Mb.tsv')
g = pd.read_csv(G, sep='\t', header=None)
print('  shape:', g.shape)
# identify p-value column: the one giving ~80 values < 5e-8
best = None
for c in range(4, 12):
    try:
        v = pd.to_numeric(g[c], errors='coerce')
        n = int((v < 5e-8).sum())
        if 20 <= n <= 300:
            print('    col %d: n(p<5e-8)=%d' % (c, n))
            if best is None or abs(n - 80) < abs(best[1] - 80):
                best = (c, n)
    except Exception:
        pass
print('  chosen p column:', best)
pcol = best[0]
g['pos38'] = pd.to_numeric(g[1], errors='coerce')
g['p'] = pd.to_numeric(g[pcol], errors='coerce')
g['rsid'] = g[g.columns[-1]].astype(str)
sig = g[g.p < 5e-8].dropna(subset=['pos38'])
print('  significant variants:', len(sig), '| span %.3f-%.3f Mb'
      % (sig.pos38.min() / 1e6, sig.pos38.max() / 1e6))

in_peak_sig = int(sig.pos38.apply(lambda x: hit(int(x))).sum())
print('  significant variants in H3K27ac peaks: %d/%d (%.1f%%)'
      % (in_peak_sig, len(sig), 100 * in_peak_sig / len(sig)))
allv = g.dropna(subset=['pos38'])
in_peak_all = int(allv.pos38.apply(lambda x: hit(int(x))).sum())
print('  all variants in file in peaks: %d/%d (%.1f%%)'
      % (in_peak_all, len(allv), 100 * in_peak_all / len(allv)))

REGIONS = {'VDR body': (47841537, 47943048), 'COL2A1 body': (47972967, 48004554),
           'VDR promoter(TSS)': (47942000, 47944049),
           'COL2A1 promoter(TSS)': (47972000, 47974000)}
for k, (a, b) in REGIONS.items():
    sub = g[(g.pos38 >= a) & (g.pos38 <= b)]
    n_peak = int(sub.pos38.apply(lambda x: hit(int(x))).sum()) if len(sub) else 0
    n_sig = int((sub.p < 5e-8).sum())
    print('  %-22s variants=%-5d inPeaks=%-5d (%.0f%%) sig=%d'
          % (k, len(sub), n_peak, 100 * n_peak / max(len(sub), 1), n_sig))

for rs, pos in (('rs11168351', 48009982), ('rs7975232', 47845054)):
    print('  %s (%d) in peak: %s' % (rs, pos, hit(pos)))

# nearest peak distances for the two sentinels
for rs, pos in (('rs11168351', 48009982), ('rs7975232', 47845054)):
    dist = min(min(abs(pos - s), abs(pos - e)) if not (s <= pos <= e) else 0 for s, e in iv)
    print('  %s nearest peak distance: %d bp' % (rs, dist))
