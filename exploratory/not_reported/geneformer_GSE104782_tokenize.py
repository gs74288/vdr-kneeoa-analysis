# -*- coding: utf-8 -*-
"""GSE104782: symbol->Ensembl and Geneformer V2 tokenization."""
import sys, io, os, pickle
import numpy as np
import anndata as ad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
BASE = r"C:\Users\liuha\OpenClawVDR\data\gf_104782"
SRC = os.path.join(BASE, 'pilot.h5ad')
RAW_DIR = os.path.join(BASE, 'raw')
OUT_DIR = os.path.join(BASE, 'tokenized')
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
print('subset:', b.shape, '| n_counts median %.0f' % np.median(b.obs['n_counts']))
b.write(os.path.join(RAW_DIR, 'pilot.h5ad'))
print('written raw/pilot.h5ad')

import geneformer  # noqa: F401
from geneformer.tokenizer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer({"state": "state", "donor": "donor"}, nproc=1)
tk.tokenize_data(RAW_DIR, OUT_DIR, "pilot", file_format="h5ad")
print('tokenized:', os.listdir(OUT_DIR))
