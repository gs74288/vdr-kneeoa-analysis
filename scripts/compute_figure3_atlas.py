# -*- coding: utf-8 -*-
"""Rebuild Figure 3 from the GSE104782 processed matrix + published cluster labels."""
import sys, io, os, gzip, re
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
DATA = r"C:\Users\liuha\OpenClawVDR\data\GSE104782"
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
os.makedirs(TAB, exist_ok=True)

# ---- load ----
mat = pd.read_csv(os.path.join(DATA, 'GSE104782_allcells_UMI_count.txt.gz'),
                  sep='\t', index_col=0)
print('matrix (genes x cells):', mat.shape)
X = mat.T.astype(np.float32)

meta = pd.read_excel(os.path.join(DATA,
       'GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx'))
meta = meta.dropna(subset=['Cluster']).set_index('Cell')
common = [c for c in X.index if c in meta.index]
print('cells with labels:', len(common))
X = X.loc[common]
meta = meta.loc[common]
meta['donor'] = [re.match(r'([A-Za-z]+\d+)', c).group(1) if re.match(r'([A-Za-z]+\d+)', c) else 'NA'
                 for c in common]

adata = ad.AnnData(X.values, obs=meta[['Cluster', 'donor']].copy(),
                   var=pd.DataFrame(index=mat.index))
adata.obs['Cluster'] = adata.obs['Cluster'].astype('category')
adata.var_names_make_unique()
print(adata)

# ---- standard pipeline ----
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
sc.tl.umap(adata, random_state=1)
print('UMAP done')

# ---- VDR detection ----
vdr_idx = list(adata.var_names).index('VDR') if 'VDR' in adata.var_names else None
if vdr_idx is None:
    print('VDR not in the filtered matrix')
    adata.obs['VDR_detected'] = 0
else:
    raw = mat.loc['VDR', common].values.astype(float)
    adata.obs['VDR_detected'] = (raw > 0).astype(int)
    adata.obs['VDR_umi'] = raw
print('VDR detected: %d / %d (%.2f%%)' % (adata.obs.VDR_detected.sum(),
      adata.n_obs, 100 * adata.obs.VDR_detected.mean()))

frac = adata.obs.groupby('Cluster', observed=True).agg(
    n=('VDR_detected', 'size'), det=('VDR_detected', 'sum'))
frac['frac'] = frac.det / frac.n
print(frac.sort_values('frac', ascending=False))

out = pd.DataFrame({'cell': adata.obs_names,
                    'UMAP1': adata.obsm['X_umap'][:, 0],
                    'UMAP2': adata.obsm['X_umap'][:, 1],
                    'Cluster': adata.obs['Cluster'].astype(str).values,
                    'donor': adata.obs['donor'].values,
                    'VDR_detected': adata.obs['VDR_detected'].values})
out.to_csv(os.path.join(TAB, 'fig3_umap.csv'), index=False)
frac.to_csv(os.path.join(TAB, 'fig3_vdr_by_state.csv'))
print('written fig3_umap.csv and fig3_vdr_by_state.csv')
