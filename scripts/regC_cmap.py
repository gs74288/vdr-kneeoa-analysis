# -*- coding: utf-8 -*-
"""CMap via Enrichr (LINCS L1000) + synovium dataset sizes."""
import sys, io, os, re, json, urllib.request, urllib.parse, urllib.error
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
OP = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {'User-Agent': 'Mozilla/5.0'}
T = r"C:\Users\liuha\OpenClawVDR\repo\tables"

genes = pd.read_csv(os.path.join(T, 'regC_markers_up.csv'), index_col=0).head(150).index.tolist()
print('genes to query:', len(genes))

print('\n=== Enrichr addList ===')
payload = urllib.parse.urlencode({'list': '\n'.join(genes),
                                  'description': 'RegC markers GSE324993'}).encode()
try:
    r = OP.open(urllib.request.Request('https://maayanlab.cloud/Enrichr/addList',
                                       data=payload, headers=UA), timeout=90)
    d = json.loads(r.read())
    uid = d['userListId']
    print('  userListId', uid)
    for lib in ('LINCS_L1000_Chem_Pert_Consensus_Sigs', 'DSigDB', 'GO_Biological_Process_2023'):
        q = urllib.parse.urlencode({'userListId': uid, 'backgroundType': lib})
        try:
            e = json.loads(OP.open(urllib.request.Request(
                'https://maayanlab.cloud/Enrichr/enrich?' + q, headers=UA), timeout=120).read())
            rows = e.get(lib, [])[:8]
            print('\n  --- %s (%d hits) ---' % (lib, len(e.get(lib, []))))
            for row in rows:
                term, p, z, comb, genes_hit = row[1], row[2], row[3], row[4], row[5]
                print('    %-58s P=%.2e  %s' % (term[:58], p, ','.join(genes_hit[:4])))
        except Exception as ex:
            print('  %s failed: %s' % (lib, type(ex).__name__))
except Exception as e:
    print('  Enrichr unavailable:', type(e).__name__, str(e)[:200])

print('\n=== synovium dataset sizes ===')
for acc in ('GSE216651', 'GSE294277', 'GSE283080', 'GSE279915'):
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/' % (acc[:-3], acc)
    try:
        html = OP.open(urllib.request.Request(url, headers=UA), timeout=60).read().decode('utf8', 'replace')
        files = re.findall(r'<a href="([^"]+)">[^<]+</a>\s+(\d{4}-\d\d-\d\d \d\d:\d\d)\s+([\d.]+[KMG]?)', html)
        print('\n  %s' % acc)
        for f in files:
            print('     %-62s %8s' % (f[0], f[2]))
        if not files:
            print('     (no files)')
    except Exception as e:
        print('  %s failed: %s' % (acc, type(e).__name__))
