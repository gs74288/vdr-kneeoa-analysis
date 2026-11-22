# -*- coding: utf-8 -*-
"""#3: core TFs of RegC across both atlases + RegC marker genes for CMap."""
import sys, io, os
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"


def prep(path):
    a = pd.read_csv(path, index_col=0)
    idx = pd.DataFrame([s.split('|') for s in a.index], columns=['donor', 'state'], index=a.index)
    a.index = a.index
    return a, idx


acts = {}
for atlas, f in (('GSE324993', 'tf_activity_gse324993.csv'),
                 ('GSE104782', 'tf_activity_gse104782.csv')):
    a, idx = prep(os.path.join(T, f))
    means = a.groupby(idx.state.values).mean()
    z = (means - means.mean()) / means.std().replace(0, np.nan)
    acts[atlas] = z
    print('%s: %d samples x %d TFs | states %s' % (atlas, a.shape[0], a.shape[1],
                                                   list(means.index)))

common = acts['GSE324993'].columns.intersection(acts['GSE104782'].columns)
print('\nshared TFs:', len(common))
if 'RegC' in acts['GSE324993'].index and 'RegC' in acts['GSE104782'].index:
    z1 = acts['GSE324993'].loc['RegC', common]
    z2 = acts['GSE104782'].loc['RegC', common]
    comb = pd.DataFrame({'z_GSE324993': z1, 'z_GSE104782': z2})
    comb['min_z'] = comb.min(axis=1)
    comb['mean_z'] = comb.mean(axis=1)
    top = comb.sort_values('min_z', ascending=False).head(25)
    print('\n=== TFs highest in RegC in BOTH atlases (ranked by worse atlas) ===')
    print(top.round(2).to_string())
    top.to_csv(os.path.join(T, 'regC_core_tfs.csv'))

# ---- RegC marker genes (pseudobulk) ----
pb = pd.read_csv(os.path.join(T, 'pseudobulk_gse324993_donor_state.csv'), index_col=0)
states = pd.Series([c.split('|')[1] for c in pb.columns], index=pb.columns)
reg = pb.loc[:, states == 'RegC'].mean(1)
oth = pb.loc[:, states != 'RegC'].mean(1)
fc = (reg - oth)
# fraction expressing
frac = (pb > 0.1).mean(1)
tab = pd.DataFrame({'mean_RegC': reg, 'mean_other': oth, 'diff': fc, 'frac_pos': frac})
tab = tab[(tab.frac_pos > 0.5)]
up = tab.sort_values('diff', ascending=False).head(300)
down = tab.sort_values('diff').head(100)
up.to_csv(os.path.join(T, 'regC_markers_up.csv'))
print('\n=== top 20 RegC-upregulated genes (pseudobulk) ===')
print(up.head(20).round(2).to_string())
print('\n138 TFs with min_z>0:', int((comb.min_z > 0).sum()),
      '| min_z>0.5:', int((comb.min_z > 0.5).sum()))
