# -*- coding: utf-8 -*-
"""Figure 7: the OA genetic signal localises to regulatory/repair chondrocytes in two atlases."""
import sys, io, os, gzip, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7,
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'axes.titlesize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7})

D = r"C:\Users\liuha\OpenClawVDR\data"
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

POS = '#B03A2E'
NEG = '#7F8C8D'
ACC = '#2C3E50'

# ---------------- GSE324993 ----------------
sc = pd.read_csv(os.path.join(D, 'scdrs_results_cell_scores.csv'))
rows, started = [], False
with gzip.open(os.path.join(r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\data",
                            'GSE324993_cell_annotation.txt.gz'), 'rt', errors='replace') as fh:
    for ln in fh:
        ln = ln.rstrip('\n')
        if ln.startswith('Cell_Label'):
            started = True
            continue
        if not started or not ln.strip():
            continue
        p = ln.split('\t')
        if len(p) < 2:
            continue
        rows.append((p[0].strip(), re.sub(r'\(.*\)', '', p[1].strip()).strip()))
ann = pd.DataFrame(rows, columns=['cell', 'state'])
ann['donor'] = ann.cell.str.split('_').str[0]
a = sc.merge(ann, on='cell', how='inner')
print('GSE324993 cells:', len(a), '| states:', a.state.nunique(), '| donors:', a.donor.nunique())

ORDER_A = ['RepC', 'RegC', 'PreHTC', 'FC', 'CPC', 'RegC2', 'MTC', 'PreFC']
LAB_A = {'RepC': 'RepC', 'RegC': 'RegC', 'PreHTC': 'preHTC', 'FC': 'FC', 'CPC': 'CPC',
         'RegC2': 'RegC2', 'MTC': 'MTC', 'PreFC': 'preFC'}
STATES_A = [s for s in ORDER_A if s in set(a.state)]

# ---------------- GSE104782 ----------------
b = pd.read_csv(os.path.join(T, 'scdrs_GSE104782_cell_scores.csv'))
b = b.rename(columns={'Cell': 'cell'})
lab = pd.read_excel(os.path.join(D, 'GSE104782',
                                 'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
lab = lab[['Cell', 'Cluster']].dropna().rename(columns={'Cell': 'cell', 'Cluster': 'state'})
lab['donor'] = lab.cell.str.split('_').str[0]
b = b.merge(lab, on='cell', how='inner')
print('GSE104782 cells:', len(b), '| states:', b.state.nunique(), '| donors:', b.donor.nunique())

ORDER_B = ['RegC', 'ProC', 'FC', 'preHTC', 'HomC', 'EC', 'HTC']
STATES_B = [s for s in ORDER_B if s in set(b.state)]

# ---------------- per-donor means ----------------
def donor_means(df, states):
    m = df.groupby(['state', 'donor']).norm_score.mean().reset_index()
    return m[m.state.isin(states)]

dm_a = donor_means(a, STATES_A)
dm_b = donor_means(b, STATES_B)
# keep only donors that contribute to that state with >=5 cells
cnt_a = a.groupby(['state', 'donor']).size().rename('n').reset_index()
cnt_b = b.groupby(['state', 'donor']).size().rename('n').reset_index()
dm_a = dm_a.merge(cnt_a, on=['state', 'donor']); dm_a = dm_a[dm_a.n >= 5]
dm_b = dm_b.merge(cnt_b, on=['state', 'donor']); dm_b = dm_b[dm_b.n >= 5]

mean_a = a.groupby('state').norm_score.mean().reindex(STATES_A)
mean_b = b.groupby('state').norm_score.mean().reindex(STATES_B)
n_a = a.groupby('state').size().reindex(STATES_A)
n_b = b.groupby('state').size().reindex(STATES_B)


def panel(ax, means, dmeans, states, hl, title, ylab):
    x = np.arange(len(states))
    cols = [POS if s == hl else NEG for s in states]
    ax.bar(x, means.values, width=0.62, color=cols, alpha=0.28, edgecolor=cols,
           linewidth=0.6, zorder=1)
    for i, s in enumerate(states):
        v = dmeans.loc[dmeans.state == s, 'norm_score'].values
        jit = np.linspace(-0.13, 0.13, len(v)) if len(v) > 1 else np.array([0.0])
        ax.scatter(np.full(len(v), i) + jit, v, s=9, color=ACC, alpha=0.85,
                   linewidths=0, zorder=3)
    ax.axhline(0, color='k', lw=0.6, zorder=2)
    lo = min(0.0, float(dmeans.norm_score.min()))
    hi = max(0.0, float(dmeans.norm_score.max()))
    span = max(hi - lo, 1e-6)
    ax.set_ylim(lo - 0.08 * span, hi + 0.12 * span)
    step = 0.5 if (hi - lo) > 1.5 else 0.25
    t0 = np.ceil(ax.get_ylim()[0] / step) * step
    t1 = np.floor(ax.get_ylim()[1] / step) * step
    ax.set_yticks(np.arange(t0, t1 + 1e-9, step))
    ax.set_xticks(x)
    ax.set_xticklabels([LAB_A.get(s, s) for s in states], rotation=30, ha='right')
    ax.set_ylabel(ylab)
    ax.set_title(title, pad=4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(0, means.values[0] + 0.03 * span, '%+.2f' % means.values[0],
            ha='center', va='bottom', fontsize=6.5, color=POS, fontweight='bold')


fig = plt.figure(figsize=(6.7, 6.5))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.25, 1.25, 1.15],
              hspace=0.62, wspace=0.28,
              left=0.085, right=0.985, top=0.945, bottom=0.075)

axA = fig.add_subplot(gs[0, :])
panel(axA, mean_a, dm_a, STATES_A, 'RepC',
      'GSE324993  (44,348 cells; 6 donors)', 'scDRS score')

axB = fig.add_subplot(gs[1, :])
panel(axB, mean_b, dm_b, STATES_B, 'RegC',
      'GSE104782  (1,464 cells; 10 donors)', 'scDRS score')

# ---------------- C: paired donor deltas ----------------
per_a = pd.read_csv(r"C:\Users\liuha\OpenClawVDR\repo\tables\scDRS_RepC_per_donor.csv", index_col=0)
per_b = pd.read_csv(os.path.join(T, 'scDRS_GSE104782_top_state_per_donor.csv'), index_col=0)

axC = fig.add_subplot(gs[2, 0])
axD = fig.add_subplot(gs[2, 1])

for ax, per, hl, lab, ptx in [
        (axC, per_a, 'RepC', 'RepC', 'P = 0.031'),
        (axD, per_b, 'RegC', 'RegC', 'P = 0.006')]:
    xs = np.arange(len(per))
    ax.plot([xs, xs], [per['other'].values, per[hl].values if hl in per else per['top'].values],
            color='#BDC3C7', lw=0.8, zorder=1)
    ax.scatter(xs, per['other'].values, s=10, color=NEG, zorder=2, linewidths=0)
    ax.scatter(xs, (per[hl] if hl in per else per['top']).values, s=10, color=POS,
               zorder=3, linewidths=0)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(per.index, rotation=45, ha='right', fontsize=6)
    ax.set_ylabel('scDRS score')
    ax.set_title('%s versus other cells, by donor' % lab, pad=4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(0.98, 0.04, ptx, transform=ax.transAxes, ha='right', fontsize=6.5)

# legend
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=POS, markersize=4, label='highlighted state'),
           Line2D([0], [0], marker='o', color='w', markerfacecolor=NEG, markersize=4, label='all other cells')]
fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False, fontsize=6.5,
           bbox_to_anchor=(0.5, 0.005))

for ax, letter in [(axA, 'A'), (axB, 'B'), (axC, 'C'), (axD, 'D')]:
    ax.text(-0.055, 1.12, letter, transform=ax.transAxes, fontsize=10, fontweight='bold',
            va='top', ha='right')

png = os.path.join(OUT, 'Figure_7.png')

# ---- layout check: no clipping, no overlapping tick labels ----
fig.canvas.draw()
rnd = fig.canvas.get_renderer()
fb = fig.bbox
issues = []
for axn, ax in (('A', axA), ('B', axB), ('C', axC), ('D', axD)):
    for t in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        if not t.get_text():
            continue
        bb = t.get_window_extent(renderer=rnd)
        if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fb.width or bb.y1 > fb.height:
            issues.append(('outside-figure', axn, t.get_text()))
    labs = [t for t in ax.get_xticklabels() if t.get_text()]
    for i in range(len(labs)):
        for j in range(i + 1, len(labs)):
            if labs[i].get_window_extent(rnd).overlaps(labs[j].get_window_extent(rnd)):
                issues.append(('overlap', axn, labs[i].get_text(), labs[j].get_text()))
print('layout issues:', issues if issues else 'none')
print('A ylim', axA.get_ylim(), '| yticks', [t.get_text() for t in axA.get_yticklabels()])
print('B ylim', axB.get_ylim(), '| yticks', [t.get_text() for t in axB.get_yticklabels()])
for t in list(axA.get_yticklabels()) + list(axB.get_yticklabels()):
    if t.get_text() == '3':
        print('  label "3" bbox', t.get_window_extent(renderer=rnd), '| figure', fb)

fig.savefig(png, dpi=300, facecolor='white')
plt.close(fig)
print('saved', png)

from PIL import Image
im = Image.open(png).convert('RGB')
print('px:', im.size)
w, h = im.size
if w != 2010:
    im = im.resize((2010, int(h * 2010 / w)), Image.LANCZOS)
tif = os.path.join(OUT, 'Figure_7.tif')
im.save(tif, format='TIFF', dpi=(300, 300), compression='tiff_lzw')
print('tiff:', im.size, '%.2f MB' % (os.path.getsize(tif) / 1e6))

# write the numbers used, for audit
with open(os.path.join(T, 'fig7_values.csv'), 'w', encoding='utf-8', newline='') as fh:
    fh.write('panel,atlas,state,n,mean\n')
    for s in STATES_A:
        fh.write('A,GSE324993,%s,%d,%.4f\n' % (LAB_A[s], n_a[s], mean_a[s]))
    for s in STATES_B:
        fh.write('B,GSE104782,%s,%d,%.4f\n' % (s, n_b[s], mean_b[s]))
print('wrote fig7_values.csv')
