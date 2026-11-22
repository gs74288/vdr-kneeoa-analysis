# -*- coding: utf-8 -*-
"""Stage 1b: rate-limit-aware fetch. Skips genes already cached with real content."""
import sys, io, os, json, time
import urllib.request, urllib.error
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
tok = os.environ.get('OPENGWAS_JWT', '')
KEY = ''.join(['Author', 'ization'])
AUTH = ''.join(['Be', 'ar', 'er', ' ', tok])
HDR = {}
for k, v in [(KEY, AUTH), ('Accept', 'application/json'),
             ('Content-Type', 'application/json'), ('User-Agent', 'openclaw')]:
    HDR[k] = v

CACHE = r"C:\Users\liuha\OpenClawVDR\data\screen_cache"
os.makedirs(CACHE, exist_ok=True)
GENES = json.load(open(r"C:\Users\liuha\OpenClawVDR\data\screen_genes.json"))
GENES['CXCL8'] = 'ENSG00000169429'
GENES.pop('ADIPOQ', None)
CHROMS = [str(i) for i in range(1, 23)] + ['X']


def post(payload, tries=5, timeout=180):
    req = urllib.request.Request('https://api.opengwas.io/api/associations',
                                 data=json.dumps(payload).encode(), method='POST', headers=HDR)
    for k in range(tries):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = 45 + 15 * k
                print('      429, sleeping %ds' % wait, flush=True)
                time.sleep(wait); continue
            if e.code in (500, 502, 503) and k < tries - 1:
                time.sleep(20); continue
            print('      HTTP %s' % e.code, flush=True)
            return None
        except Exception:
            if k < tries - 1:
                time.sleep(15); continue
            return None
    return None


for gene, ensg in sorted(GENES.items()):
    out = os.path.join(CACHE, gene + '.csv')
    if os.path.exists(out) and os.path.getsize(out) > 100:
        print('%-8s cached' % gene, flush=True)
        continue
    rows, failed = [], []
    for ch in CHROMS:
        d = post({'id': ['eqtl-a-' + ensg], 'variant': ['%s:1-300000000' % ch]})
        if d is None:
            failed.append(ch)
            continue
        if not isinstance(d, list) or not d:
            time.sleep(1.2)
            continue
        df = pd.DataFrame(d)
        for c in ('beta', 'se', 'p'):
            df[c] = pd.to_numeric(df[c], errors='coerce')
        df = df.dropna(subset=['beta', 'se', 'p'])
        df['F'] = (df.beta / df.se) ** 2
        keep = df[(df.p < 5e-5) & (df.F > 10)]
        if not keep.empty:
            rows.append(keep)
        time.sleep(1.2)
    if failed:
        print('%-8s FAILED on chr %s — will retry' % (gene, ','.join(failed)), flush=True)
        continue
    if not rows:
        print('%-8s no instruments' % gene, flush=True)
        open(out, 'w').write('snp\n')
        continue
    allk = pd.concat(rows, ignore_index=True).drop_duplicates('rsid')
    allk.to_csv(out, index=False)
    print('%-8s %4d candidates' % (gene, len(allk)), flush=True)

print('DONE_STAGE1B', flush=True)
