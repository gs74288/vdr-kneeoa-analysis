# -*- coding: utf-8 -*-
"""VDR regulon activity by donor-state pseudobulk (decoupleR ULM), both atlases."""
import sys, io, os, re, inspect
import numpy as np
import pandas as pd
import anndata as ad
import decoupler as dc
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

net = pd.read_csv(os.path.join(D, 'collectri_human.csv'))
net = net[['source', 'target', 'weight']].dropna()
print('net:', net.shape, '| TFs', net.source.nunique())

print('ulm signature:', inspect.signature(dc.mt.ulm))


def run_ulm(mat, net):
    # decoupler expects samples in rows and features in columns
    out = dc.mt.ulm(data=mat.T, net=net, tmin=5, verbose=False)
    if isinstance(out, tuple):
        return out[0]
    return out


# ---------- GSE104782 ----------
print('\n=== GSE104782 pseudobulk ===')
cnt = pd.read_csv(os.path.join(D, 'GSE104782', 'GSE104782_allcells_UMI_count.txt.gz'),
                  sep='\t', index_col=0)
if cnt.shape[0] > cnt.shape[1]:
    cnt = cnt.T
cnt = cnt.groupby(level=0).sum()
lab = pd.read_excel(os.path.join(D, 'GSE104782',
                                 'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
lab = lab[['Cell', 'Cluster']].dropna()
pos = {c: i for i, c in enumerate(cnt.index)}
cells = [c for c in lab.Cell if c in pos]
X = cnt.values[[pos[c] for c in cells]].astype(np.float32)
X = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e4)
genes_b = list(cnt.columns)
states_b = lab.set_index('Cell').loc[cells, 'Cluster'].values
donors_b = np.array([c.split('_')[0] for c in cells])

units = {}
for dn in sorted(set(donors_b)):
    for st in sorted(set(states_b)):
        m = (donors_b == dn) & (states_b == st)
        if m.sum() >= 10:
            units['%s|%s' % (dn, st)] = X[m].mean(0)
mat_b = pd.DataFrame(np.vstack(list(units.values())).T, index=genes_b, columns=list(units))
print('pseudobulk:', mat_b.shape)
act_b = run_ulm(mat_b, net)
print('activity:', act_b.shape)
act_b.to_csv(os.path.join(T, 'tf_activity_gse104782.csv'))

# ---------- GSE324993 ----------
print('\n=== GSE324993 pseudobulk (chunked) ===')
a = ad.read_h5ad(os.path.join(D, 'GSE324993.h5ad'), backed='r')
genes_a = list(a.var_names)
n_cells, n_genes = a.shape
states_a = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in a.obs['cell_type'].values])
donors_a = np.array([s.split('_')[0] for s in a.obs_names])

keys = []
for dn in sorted(set(donors_a)):
    for st in sorted(set(states_a)):
        m = (donors_a == dn) & (states_a == st)
        if m.sum() >= 10:
            keys.append('%s|%s' % (dn, st))
masks = {k: (donors_a == k.split('|')[0]) & (states_a == k.split('|')[1]) for k in keys}
counts = {k: int(masks[k].sum()) for k in keys}
print('%d units' % len(keys))

acc = {k: np.zeros(len(genes_a), dtype=np.float64) for k in keys}
CH = 1500
for start in range(0, n_genes, CH):
    end = min(start + CH, n_genes)
    Xc = a[:, start:end].to_memory().X
    Xc = Xc.toarray() if hasattr(Xc, 'toarray') else np.asarray(Xc)
    Xc = Xc.astype(np.float32)
    for k in keys:
        acc[k][start:end] += Xc[masks[k]].sum(0)
    print('  chunk %d-%d' % (start, end), flush=True)

mat_a = pd.DataFrame({k: acc[k] / counts[k] for k in keys}, index=genes_a)
print('pseudobulk:', mat_a.shape)
mat_a.to_csv(os.path.join(T, 'pseudobulk_gse324993_donor_state.csv'))
act_a = run_ulm(mat_a, net)
print('activity:', act_a.shape)
act_a.to_csv(os.path.join(T, 'tf_activity_gse324993.csv'))

# ---------- analysis ----------
def summarise(act, prefix, focus, atlas):
    # act: samples x sources
    idx = pd.DataFrame([s.split('|') for s in act.index], columns=['donor', 'state'],
                       index=act.index)
    if 'VDR' not in act.columns:
        print('%s: VDR not in activity matrix' % atlas)
        return
    rows = []
    for dn in sorted(idx.donor.unique()):
        m = idx.donor == dn
        f = m & (idx.state == focus)
        o = m & (idx.state != focus)
        if f.sum() and o.sum():
            rows.append((dn, float(act.loc[f, 'VDR'].mean()), float(act.loc[o, 'VDR'].mean())))
    r = pd.DataFrame(rows, columns=['donor', 'focus', 'other'])
    r['delta'] = r.focus - r.other
    print('\n--- %s: VDR regulon activity, %s vs other states ---' % (atlas, focus))
    print(r.round(3).to_string(index=False))
    if len(r) >= 3:
        w = stats.wilcoxon(r.focus, r.other)
        print('positive %d/%d | mean delta %+.3f | Wilcoxon P=%.4g'
              % (int((r.delta > 0).sum()), len(r), r.delta.mean(), w.pvalue))
    # rank of VDR among TFs in focus units
    f = idx.state == focus
    ranks = []
    for s in act.index[f]:
        v = act.loc[s].rank(ascending=False)
        ranks.append(v['VDR'])
    print('VDR rank among %d TFs in %s units: median %.0f (range %.0f-%.0f)'
          % (act.shape[1], focus, np.median(ranks), min(ranks), max(ranks)))


summarise(act_a, 'a', 'RegC', 'GSE324993 / RegC')
summarise(act_a, 'a', 'RepC', 'GSE324993 / RepC')
summarise(act_b, 'b', 'RegC', 'GSE104782 / RegC')
