# -*- coding: utf-8 -*-
"""Conditional (COJO-style) analysis of the OA signal at 12q13.

Conditions the knee OA association on rs11168351 using 1000G EUR LD, then exports
the conditional summary statistics for a repeat colocalisation.
"""
import sys, io, os, gzip, struct, zlib
import numpy as np
import pandas as pd
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
DIRECT = urllib.request.build_opener(urllib.request.ProxyHandler({}))
LD = r"C:\Users\liuha\OpenClawVDR\data\ld1000g"
OUT = os.path.join(TAB, 'cojo')
os.makedirs(OUT, exist_ok=True)
VCF = ('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
       'ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')
PEAK = 48_009_982
OFFSET = 393_783
LO38, HI38 = PEAK - 100_000, PEAK + 100_000
COND_RS = 'rs11168351'

ci = pd.read_csv(os.path.join(TAB, 'coloc_harmonised_input_GCST90566800.csv'))
w = ci[(ci.pos >= LO38) & (ci.pos <= HI38)].copy()
print('variants in window:', len(w))

# GRCh38 -> GRCh37 for the 1000G fetch
lo37, hi37 = LO38 + OFFSET, HI38 + OFFSET


def reg2bins(beg, end):
    end -= 1
    out = [0]
    for sh, base in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        for k in range(1 + (beg >> sh), 1 + (end >> sh) + 1):
            out.append(base + k)
    return out


tbi = gzip.decompress(open(os.path.join(LD, 'chr12.vcf.gz.tbi'), 'rb').read())
off = 4
v = struct.unpack_from('<8i', tbi, off)
n_ref, l_nm = v[0], v[7]
off += 32
names = tbi[off:off + l_nm].rstrip(b'\x00').decode().split('\x00')
off += l_nm
my = {}
for i in range(n_ref):
    n_bin = struct.unpack_from('<i', tbi, off)[0]
    off += 4
    bins = {}
    for _ in range(n_bin):
        b, nc = struct.unpack_from('<Ii', tbi, off)
        off += 8
        bins[b] = [struct.unpack_from('<QQ', tbi, off + 16 * k) for k in range(nc)]
        off += 16 * nc
    n_intv = struct.unpack_from('<i', tbi, off)[0]
    off += 4 + 8 * n_intv
    if names[i] == '12':
        my = bins
        break

chunks = []
for b in reg2bins(lo37, hi37):
    if b in my:
        chunks.extend(my[b])
chunks.sort()
merged = []
for c in chunks:
    if merged and c[0] <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], c[1])
    else:
        merged.append(list(c))

samples = open(os.path.join(LD, 'samples.txt')).read().split('\n')
panel = pd.read_csv(os.path.join(LD, 'panel.txt'), sep='\t')
eur = panel.loc[panel['super_pop'] == 'EUR', 'sample'].tolist()
IDX = {s: i for i, s in enumerate(samples)}
EUR = [IDX[s] for s in eur if s in IDX]

geno = {}
for c in merged:
    a, b = c[0] >> 16, (c[1] >> 16) + 65536
    req = urllib.request.Request(VCF, headers={'User-Agent': 'openclaw',
                                               'Range': 'bytes=%d-%d' % (a, b - 1)})
    blob = DIRECT.open(req, timeout=300).read()
    rest, buf = blob, []
    while rest:
        d = zlib.decompressobj(31)
        try:
            buf.append(d.decompress(rest))
        except Exception:
            break
        rest = d.unused_data
        if not rest:
            break
    for ln in b''.join(buf).decode('utf-8', 'replace').split('\n'):
        if not ln or ln.startswith('#'):
            continue
        f = ln.split('\t')
        if len(f) < 10 or f[0] != '12':
            continue
        try:
            p = int(f[1])
        except ValueError:
            continue
        if p < lo37 or p > hi37 or p in geno:
            continue
        gts = f[9:]
        vec, ok = [], True
        for cc in EUR:
            g = gts[cc].split(':')[0]
            if g == '0|0':
                vec.append(0)
            elif g in ('0|1', '1|0'):
                vec.append(1)
            elif g == '1|1':
                vec.append(2)
            else:
                ok = False
                break
        if ok:
            geno[p] = np.array(vec, float)
print('genotypes:', len(geno))

w['pos37'] = w.pos + OFFSET
w = w[w.pos37.isin(geno.keys())].copy()
print('matched variants:', len(w))
pos = w.pos37.values
G = np.vstack([geno[p] for p in pos])
Gc = G - G.mean(1, keepdims=True)
sd = Gc.std(1, keepdims=True)
sd[sd == 0] = 1
Z = Gc / sd
R = (Z @ Z.T) / Z.shape[1]
np.fill_diagonal(R, 1.0)

cond_idx = list(w.snp).index(COND_RS) if COND_RS in list(w.snp) else None
print('conditioning variant index:', cond_idx)
assert cond_idx is not None

b = w.b2h.values
se = w.se2.values
b_k, se_k = b[cond_idx], se[cond_idx]
r_k = R[:, cond_idx]

b_cond = b - (r_k * se / se_k) * b_k
se_cond = se * np.sqrt(np.clip(1 - r_k ** 2, 1e-12, None))
z_cond = b_cond / se_cond
p_cond = 2 * (1 - __import__('scipy.stats', fromlist=['norm']).norm.cdf(np.abs(z_cond)))

out = pd.DataFrame({'snp': w.snp, 'pos': w.pos, 'pos37': w.pos37,
                    'b1': w.b1, 'se1': w.se1, 'p1': w.p1,
                    'b2_orig': b, 'se2_orig': se, 'p2_orig': w.p2,
                    'b2_cond': b_cond, 'se2_cond': se_cond, 'p2_cond': p_cond,
                    'r_to_cond': r_k})
out.to_csv(os.path.join(OUT, 'conditional_oa.csv'), index=False)
np.save(os.path.join(OUT, 'ld_region.npy'), R)
open(os.path.join(OUT, 'ld_snps.txt'), 'w').write('\n'.join(w.snp) + '\n')

print('\n=== after conditioning on %s ===' % COND_RS)
print('  original OA p range : %.2e .. %.2e' % (w.p2.min(), w.p2.max()))
print('  conditional p range : %.2e .. %.2e' % (p_cond.min(), p_cond.max()))
print('  variants still p<5e-8 after conditioning: %d' % int((p_cond < 5e-8).sum()))
top = out.nsmallest(5, 'p2_cond')[['snp', 'p2_orig', 'p2_cond', 'r_to_cond']]
print(top.to_string(index=False))
print('\nwritten conditional_oa.csv')
