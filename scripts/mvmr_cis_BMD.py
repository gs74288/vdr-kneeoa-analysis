# -*- coding: utf-8 -*-
"""Redo the cis MVMR with the 2025 outcome (GCST90566800)."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"

new = pd.read_csv(os.path.join(TAB, 'cisMR_GCST90566800_instruments.csv'))
old = pd.read_csv(os.path.join(TAB, 'mvmr_cis_input.csv'))

d = old.merge(new[['snp', 'bo', 'seo', 'p_out']], on='snp', how='left')
print('=== harmonised (new outcome) ===')
print(d[['snp', 'bx', 'bm', 'bo', 'seo', 'p_out']].to_string(index=False))
assert d.bo.notna().all()

# univariable (correct IVW weights: 1 / se_wald^2)
wald = d.bo / d.bx
se_wald = (d.seo / d.bx).abs()
w = 1 / se_wald ** 2
b_ivw = (w * wald).sum() / w.sum()
se_fixed = np.sqrt(1 / w.sum())
q = (w * (wald - b_ivw) ** 2).sum()
dfree_q = len(d) - 1
phi = max(1.0, q / dfree_q)
se_r = se_fixed * np.sqrt(phi)
print('\nunivariable IVW (new outcome):')
print('  fixed   OR %.3f (%.3f-%.3f)' % (np.exp(b_ivw),
      np.exp(b_ivw - 1.96 * se_fixed), np.exp(b_ivw + 1.96 * se_fixed)))
print('  random  OR %.3f (%.3f-%.3f)  Q=%.2f df=%d P=%.4f' % (
      np.exp(b_ivw), np.exp(b_ivw - 1.96 * se_r), np.exp(b_ivw + 1.96 * se_r),
      q, dfree_q, 1 - stats.chi2.cdf(q, dfree_q)))

# MVMR-IVW
X = np.column_stack([d.bx.values, d.bm.values])
Y = d.bo.values
W = np.diag(1 / d.seo.values ** 2)
XtWX = X.T @ W @ X
theta = np.linalg.solve(XtWX, X.T @ W @ Y)
resid = Y - X @ theta
df_ov = len(Y) - X.shape[1]
sigma2 = float(resid.T @ W @ resid) / max(df_ov, 1)
cov = sigma2 * np.linalg.inv(XtWX)
se = np.sqrt(np.diag(cov))

print('\n=== MVMR-IVW (new outcome, VDR | BMD) ===')
for i, lab in enumerate(('VDR expression', 'BMD')):
    z = theta[i] / se[i]
    print('  %-15s OR %.3f (95%% CI %.3f-%.3f)  P = %.3f' % (
        lab, np.exp(theta[i]), np.exp(theta[i] - 1.96 * se[i]),
        np.exp(theta[i] + 1.96 * se[i]), 2 * (1 - stats.norm.cdf(abs(z)))))
print('  overidentification df =', df_ov)

print('\nconditional strength (unchanged by outcome):')
for target, other, lab in ((0, 1, 'VDR'), (1, 0, 'BMD')):
    Xo = X[:, [other]]
    fit = X[:, target] - Xo @ np.linalg.lstsq(Xo, X[:, target], rcond=None)[0]
    print('  %-4s conditional F = %.2f' %
          (lab, (fit ** 2 / (d.seo.values ** 2)).sum() / max(len(d) - 2, 1)))

print('\nBMD association of each instrument (from OpenGWAS, n=583,314):')
for _, r in d.iterrows():
    print('  %-12s beta %+.5f (se %.5f)' % (r.snp, r.bm, r.sm))
