# -*- coding: utf-8 -*-
"""Pilot: convert symbols to Ensembl and tokenize 1,000 cells with Geneformer V2."""
import sys, io, os, pickle
import numpy as np
import pandas as pd
import anndata as ad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
SRC = r"C:\Users\liuha\OpenClawVDR\data\gf_pilot\pilot.h5ad"
RAW_DIR = r"C:\Users\liuha\OpenClawVDR\data\gf_pilot\raw"
OUT_DIR = r"C:\Users\liuha\OpenClawVDR\data\gf_pilot\tokenized"
DICT = r"C:\Users\liuha\OpenClawVDR\tools\geneformer\src\geneformer\gene_name_id_dict_gc104M.pkl"
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

with open(DICT, 'rb') as fh:
    name2id = pickle.load(fh)
print('dict size:', len(name2id))

a = ad.read_h5ad(SRC)
mapped = [name2id.get(str(g)) for g in a.var_names]
keep = [i for i, m in enumerate(mapped) if m]
print('genes mapped: %d of %d' % (len(keep), a.n_vars))
b = a[:, keep].copy()
ids = [mapped[i] for i in keep]
b.var['ensembl_id'] = ids
b.var_names = ids
b.var_names_make_unique()
b.obs['n_counts'] = np.asarray(b.X.sum(axis=1)).ravel()
print('subset:', b.shape)

# Geneformer expects raw counts in X for V2 tokenizer with keep_counts
b.write(os.path.join(RAW_DIR, 'pilot.h5ad'))
print('written raw/pilot.h5ad')

from geneformer.tokenizer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer({"state": "state", "donor": "donor"}, nproc=1)
tk.tokenize_data(RAW_DIR, OUT_DIR, "pilot", file_format="h5ad")
print('tokenized files:', os.listdir(OUT_DIR))
for f in os.listdir(OUT_DIR):
    print('  ', f, os.path.getsize(os.path.join(OUT_DIR, f)) / 1e6, 'MB')

