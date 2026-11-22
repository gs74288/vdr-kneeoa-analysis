# -*- coding: utf-8 -*-
"""Build the COJO sumstats file (IDs as chr:pos, GRCh37)."""
import sys, io, os
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\OpenClawVDR\tools\cojo"
OFFSET = 393_783
N = 1182675

d = pd.read_csv(os.path.join(TAB, 'coloc_harmonised_input_GCST90566800.csv'))
d['pos37'] = d.pos + OFFSET
d['SNP'] = d.pos37.map(lambda p: '12:%d' % p)

out = pd.DataFrame({
    'SNP': d.SNP,
    'A1': d.ea2.str.upper(),
    'A2': d.oa2.str.upper(),
    'freq': d.eaf2,
    'b': d.b2h,
    'se': d.se2,
    'p': d.p2,
    'N': N,
})
out = out.dropna()
out = out[out.freq.between(0.001, 0.999)]
# keep SNPs present in the reference panel
bim = pd.read_csv(os.path.join(OUT, 'ref.bim'), sep='\t', header=None,
                  names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
out = out[out.SNP.isin(set(bim.snp))]
out.to_csv(os.path.join(OUT, 'cojo_sumstats.txt'), sep=' ', index=False)
print('sumstats rows:', len(out))

cond_id = '12:%d' % (48_009_982 + OFFSET)
open(os.path.join(OUT, 'cond.txt'), 'w').write(cond_id + '\n')
print('conditioning SNP:', cond_id, '| present in sumstats:', cond_id in set(out.SNP))
