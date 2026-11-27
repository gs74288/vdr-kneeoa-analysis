# -*- coding: utf-8 -*-
"""Pilot embedding extraction with Geneformer V2-104M on 1,000 chondrocytes."""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, r"C:\Users\liuha\OpenClawVDR\tools\geneformer\src")
import geneformer

MODEL = r"C:\Users\liuha\OpenClawVDR\tools\geneformer\model_v2_104M\Geneformer-V2-104M"
DATA = r"C:\Users\liuha\OpenClawVDR\data\gf_104782\tokenized\pilot.dataset"
OUT = r"C:\Users\liuha\OpenClawVDR\data\gf_104782\embs"
os.makedirs(OUT, exist_ok=True)

ex = geneformer.EmbExtractor(model_version="V2", emb_mode="cls", max_ncells=1000,
                             forward_batch_size=4, nproc=1, emb_label=["state", "donor"])
ex.extract_embs(MODEL, DATA, OUT, "pilot")
print('done')
for f in os.listdir(OUT):
    print('  ', f, os.path.getsize(os.path.join(OUT, f)) / 1e6, 'MB')

