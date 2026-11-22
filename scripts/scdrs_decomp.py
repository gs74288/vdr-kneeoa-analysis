# -*- coding: utf-8 -*-
"""scDRS decomposition: standardised-expression contributions + VDR+/VDR- within state."""
import sys, io, os, re
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

z = pd.read_csv(os.path.join(D, 'GCST007090_gene_z.tsv'), sep='\t').dropna()
z = z[~z.GENE.duplicated()]
W = dict(zip(z.GENE, z.Z.astype(float)))
print('gene set:', len(W))


def build(keep, wv, me_state, gmean, gsd, n, corr):
    """contribution = weight * standardised mean expression in that state."""
    zs = (me_state - gmean) / np.maximum(gsd, 1e-9)
    con = wv * zs
    tot = np.abs(con).sum()
    sh = 100 * np.abs(con) / tot if tot > 0 else con * 0
    return pd.DataFrame({'gene': keep, 'weight': wv, 'mean_expr': me_state,
                         'z_expr': zs, 'contribution': con, 'share_pct': sh, 'corr': corr})


def report(df, atlas, focus, ncell, out):
    df.to_csv(out, index=False)
    d = df.reindex(df.contribution.abs().sort_values(ascending=False).index)
    print('\n--- %s %s: top 12 by |contribution| (n=%d cells) ---' % (atlas, focus, ncell))
    print(d.head(12)[['gene', 'z_expr', 'share_pct', 'corr']].round(3).to_string(index=False))
    for g in ('VDR', 'COL2A1'):
        r = df[df.gene == g]
        if len(r):
            rank = int((df.contribution.abs() > abs(r.contribution.values[0])).sum()) + 1
            print('  %-7s rank %d/%d | share %.2f%% | z %+.3f | corr %+.3f'
                  % (g, rank, len(df), r.share_pct.values[0],
                     r.z_expr.values[0], r['corr'].values[0]))
        else:
            print('  %s absent' % g)
    print('written', os.path.basename(out))


def donor_contrast(states, donors, scores, vdr, focus, atlas):
    fm = states == focus
    vd, sc, dn = vdr[fm], scores[fm], donors[fm]
    per = []
    for d in sorted(set(dn)):
        mm = dn == d
        p, ng = sc[mm & vd], sc[mm & ~vd]
        if len(p) >= 5 and len(ng) >= 5:
            per.append(dict(donor=d, n_pos=len(p), n_neg=len(ng),
                            mean_pos=p.mean(), mean_neg=ng.mean(), delta=p.mean() - ng.mean()))
    per = pd.DataFrame(per)
    print('\n--- %s %s: VDR+ vs VDR- within state ---' % (atlas, focus))
    print(per.round(3).to_string(index=False) if len(per) else 'not enough donors')
    if len(per) >= 3:
        wt = stats.wilcoxon(per.mean_pos, per.mean_neg)
        tt = stats.ttest_rel(per.mean_pos, per.mean_neg)
        print('donors positive: %d/%d | mean delta %+.3f | Wilcoxon P = %.3g | t P = %.3g'
              % (int((per.delta > 0).sum()), len(per), per.delta.mean(), wt.pvalue, tt.pvalue))
        per['state'] = focus
        per['atlas'] = atlas
        per.to_csv(os.path.join(T, 'scdrs_vdr_within_%s.csv' % atlas.lower()), index=False)


# ================= GSE104782 =================
print('\n================ GSE104782 ================')
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
ci = np.array([pos[c] for c in cells])
X = cnt.values[ci].astype(np.float32)
X = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e4)
genes = list(cnt.columns)
scores = pd.read_csv(os.path.join(T, 'scdrs_GSE104782_cell_scores.csv')) \
    .rename(columns={'Cell': 'cell'}).set_index('cell').loc[cells].norm_score.values
states = lab.set_index('Cell').loc[cells, 'Cluster'].values
donors = np.array([c.split('_')[0] for c in cells])
vdr = X[:, genes.index('VDR')] > 0

gi = {g: i for i, g in enumerate(genes)}
keep = [g for g in genes if g in W]
cols = np.array([gi[g] for g in keep])
wv = np.array([W[g] for g in keep])
Xk = X[:, cols]
gmean, gsd = Xk.mean(0), Xk.std(0)

rows = []
for st in sorted(set(states)):
    m = states == st
    if m.sum() < 30:
        continue
    Xs = Xk[m]
    y = scores[m]
    me = Xs.mean(0)
    ys = y - y.mean()
    xs = Xs - me
    den = np.sqrt((xs ** 2).sum(0) * (ys ** 2).sum())
    cor = np.where(den > 0, (xs * ys[:, None]).sum(0) / np.maximum(den, 1e-12), 0)
    df = build(keep, wv, me, gmean, gsd, m.sum(), cor)
    df['state'] = st
    df['atlas'] = 'GSE104782'
    rows.append(df)
dfb = pd.concat(rows, ignore_index=True)
report(dfb[dfb.state == 'RegC'], 'GSE104782', 'RegC', int((states == 'RegC').sum()),
       os.path.join(T, 'scdrs_contribution_gse104782.csv'))
dfb.to_csv(os.path.join(T, 'scdrs_contribution_gse104782_allstates.csv'), index=False)
donor_contrast(states, donors, scores, vdr, 'RegC', 'GSE104782')

# ================= GSE324993 (chunked, backed) =================
print('\n================ GSE324993 ================')
a = ad.read_h5ad(os.path.join(D, 'GSE324993.h5ad'), backed='r')
genes = list(a.var_names)
n_cells, n_genes = a.shape
print('shape:', a.shape)
states = np.array([re.sub(r'\(.*\)', '', str(s)).strip() for s in a.obs['cell_type'].values])
donors = np.array([s.split('_')[0] for s in a.obs_names])
scores = pd.read_csv(os.path.join(D, 'scdrs_results_cell_scores.csv')) \
    .set_index('cell').loc[a.obs_names].norm_score.values
ustates = sorted(set(states))
masks = {s: (states == s) for s in ustates}
nst = {s: int(masks[s].sum()) for s in ustates}
print('states:', nst)

keep = [g for g in genes if g in W]
kk = {g: i for i, g in enumerate(keep)}
wv = np.array([W[g] for g in keep])
G = len(keep)
sum_x = {s: np.zeros(G) for s in ustates}
sum_xy = {s: np.zeros(G) for s in ustates}
sum_x2 = {s: np.zeros(G) for s in ustates}
sum_y = {}; sum_y2 = {}
for s in ustates:
    y = scores[masks[s]]
    sum_y[s] = float(y.sum()); sum_y2[s] = float((y ** 2).sum())

vdr_vec = np.zeros(n_cells, dtype=bool)
vdr_idx = genes.index('VDR')
CH = 1200
for start in range(0, n_genes, CH):
    end = min(start + CH, n_genes)
    sub = a[:, start:end].to_memory()
    Xc = sub.X
    Xc = Xc.toarray() if hasattr(Xc, 'toarray') else np.asarray(Xc, dtype=np.float32)
    Xc = Xc.astype(np.float32)
    if start <= vdr_idx < end:
        vdr_vec = Xc[:, vdr_idx - start] > 0
    names = genes[start:end]
    loc = [(i, kk[g]) for i, g in enumerate(names) if g in kk]
    if loc:
        src = np.array([l[0] for l in loc])
        dst = np.array([l[1] for l in loc])
        for s in ustates:
            m = masks[s]
            if nst[s] < 30:
                continue
            Xs = Xc[m][:, src]
            sum_x[s][dst] += Xs.sum(0)
            sum_x2[s][dst] += (Xs ** 2).sum(0)
            sum_xy[s][dst] += (Xs * scores[m][:, None]).sum(0)
    print('  chunk %d-%d' % (start, end), flush=True)

N = float(n_cells)
gmean = sum(sum_x[s] for s in ustates) / N
gvar = sum(sum_x2[s] for s in ustates) / N - gmean ** 2
gsd = np.sqrt(np.maximum(gvar, 0))

rows = []
for s in ustates:
    if nst[s] < 30:
        continue
    n = nst[s]
    me = sum_x[s] / n
    sy, sy2 = sum_y[s], sum_y2[s]
    cov = sum_xy[s] - sum_x[s] * sy / n
    vx = sum_x2[s] - sum_x[s] ** 2 / n
    vy = sy2 - sy ** 2 / n
    den = np.sqrt(np.maximum(vx, 0) * np.maximum(vy, 0))
    cor = np.where(den > 0, cov / np.maximum(den, 1e-12), 0)
    df = build(keep, wv, me, gmean, gsd, n, cor)
    df['state'] = s
    df['atlas'] = 'GSE324993'
    rows.append(df)
dfa = pd.concat(rows, ignore_index=True)
report(dfa[dfa.state == 'RepC'], 'GSE324993', 'RepC', nst['RepC'],
       os.path.join(T, 'scdrs_contribution_gse324993.csv'))
dfa.to_csv(os.path.join(T, 'scdrs_contribution_gse324993_allstates.csv'), index=False)
donor_contrast(states, donors, scores, vdr_vec, 'RepC', 'GSE324993')
