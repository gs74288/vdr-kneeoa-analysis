# -*- coding: utf-8 -*-
"""Build a plink-format 1000G EUR reference for the 12q13 region (for GCTA-COJO)."""
import sys, io, os, gzip, struct, zlib
import urllib.request
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
DIRECT = urllib.request.build_opener(urllib.request.ProxyHandler({}))
LD = r"C:\Users\liuha\OpenClawVDR\data\ld1000g"
OUT = r"C:\Users\liuha\OpenClawVDR\tools\cojo"
os.makedirs(OUT, exist_ok=True)
VCF = ('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
       'ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz')
LO37, HI37 = 48_300_000, 48_510_000


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
for b in reg2bins(LO37, HI37):
    if b in my:
        chunks.extend(my[b])
chunks.sort()
merged = []
for c in chunks:
    if merged and c[0] <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], c[1])
    else:
        merged.append(list(c))

lines, seen = [], set()
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
        if p < LO37 or p > HI37 or p in seen:
            continue
        if ',' in f[4]:
            continue                      # drop multi-allelic
        if len(f[3]) != 1 or len(f[4]) != 1:
            continue                      # SNPs only
        seen.add(p)
        lines.append('\t'.join(f[:8] + f[8:]))

samples = open(os.path.join(LD, 'samples.txt')).read().split('\n')
samples = [s for s in samples if s]
print('variants:', len(lines), '| samples:', len(samples))

hdr = ['##fileformat=VCFv4.1']
for ch in ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14',
           '15', '16', '17', '18', '19', '20', '21', '22', 'X'):
    hdr.append('##contig=<ID=%s>' % ch)
hdr.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples))
open(os.path.join(OUT, 'region.vcf'), 'w', newline='').write(
    '\n'.join(hdr + lines) + '\n')

panel = pd.read_csv(os.path.join(LD, 'panel.txt'), sep='\t')
eur = panel.loc[panel['super_pop'] == 'EUR', ['sample', 'sample']]
eur.to_csv(os.path.join(OUT, 'eur.keep'), sep='\t', index=False, header=False)
print('wrote region.vcf and eur.keep (%d EUR)' % len(eur))
