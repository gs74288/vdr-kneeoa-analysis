# -*- coding: utf-8 -*-
"""Clean summary of the conditional analysis with sensible filters."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables\cojo"
d = pd.read_csv(os.path.join(OUT, 'conditional_oa.csv'))

bk = d.b2_orig.values[84]
sek = d.se2_orig.values[84]
r = d.r_to_cond.values
se = d.se2_orig.values
b = d.b2_orig.values

b_cond = b - (r * se / sek) * bk
se_cond = se * np.sqrt(np.clip(1 - r ** 2, 1e-6, None))
z = b_cond / se_cond
p_cond = 2 * (1 - stats.norm.cdf(np.abs(z)))
d['b_cond'] = b_cond
d['se_cond'] = se_cond
d['p_cond'] = p_cond
d['abs_r'] = np.abs(r)

# keep common, well-measured, non-collinear variants
m = (d.abs_r < 0.8) & (d.se2_orig > 0.002) & (d.se2_orig < 0.05)
s = d[m].copy()

print('=== conditional analysis, %d informative variants ===' % len(s))
print('  original  min p : %.2e' % s.p2_orig.min())
print('  conditional min p: %.2e' % s.p_cond.min())
print('  still p<5e-8 : %d' % int((s.p_cond < 5e-8).sum()))
print()
print('top residual signals:')
print(s.nsmallest(5, 'p_cond')[['snp', 'pos', 'p2_orig', 'p_cond', 'abs_r']].to_string(index=False))

v = d[(d.pos >= 47_841_537) & (d.pos <= 47_943_048)]
print('\nVDR gene body (all variants):')
print('  n = %d, original min p = %.2e, conditional min p = %.2e'
      % (len(v), v.p2_orig.min(), v.p_cond.clip(upper=1).min()))
print('  variants with conditional p < 0.05: %d'
      % int((v.p_cond.clip(upper=1) < 0.05).sum()))

c = d[(d.pos >= 47_972_967) & (d.pos <= 48_004_554)]
print('\nCOL2A1 gene body:')
print('  n = %d, original min p = %.2e, conditional min p = %.2e'
      % (len(c), c.p2_orig.min(), c.p_cond.clip(upper=1).min()))

d.to_csv(os.path.join(OUT, 'conditional_oa_final.csv'), index=False)
print('\nwritten conditional_oa_final.csv')
