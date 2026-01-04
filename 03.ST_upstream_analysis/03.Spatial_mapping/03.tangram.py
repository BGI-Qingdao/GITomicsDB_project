import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
import tangram as tg

ad_sc=sc.read_h5ad(sys.argv[1])
ad_sc = ad_sc[~ad_sc.obs['group'].isin(["Discard"]), :]
ad_sp=sc.read_h5ad(sys.argv[2])

xs = ad_sp.obs.new_x.values
ys = ad_sp.obs.new_y.values
plt.axis('off')
plt.scatter(xs, ys, s=.7);
plt.gca().invert_yaxis()

sc.pp.normalize_total(ad_sc)
tg.pp_adatas(ad_sc, ad_sp, genes=None)
assert ad_sc.uns['training_genes'] == ad_sp.uns['training_genes']

ad_map = tg.map_cells_to_space(ad_sc, ad_sp, mode='clusters', cluster_label='celltype')

tg.project_cell_annotations(ad_map, ad_sp, annotation="celltype")

annotation_list = list(pd.unique(ad_sc.obs['celltype']))
tg.plot_cell_annotation_sc(ad_sp, annotation_list, x='new_x', y='new_y',spot_size= 50, scale_factor=0.1, perc=0.001)

plt.savefig(sys.argv[3] + ".map.pdf", format='pdf', bbox_inches='tight')
plt.close()

ad_sp.write_h5ad(sys.argv[3]+".spatial.h5ad")
ad_map.write_h5ad(sys.argv[3]+".mapped.h5ad")

#pred
df = ad_sp.obsm['tangram_ct_pred']
df.to_csv(sys.argv[3]+'.pred.txt', sep='\t', header=True, index=True)

#max
max_values = df.max(axis=1)
max_columns = df.idxmax(axis=1)
ad_sp.obs['max_celltype']=max_columns
ad_sp.obs['max_celltype_values']=max_values
meta=ad_sp.obs
meta.to_csv(sys.argv[3]+'.meta.txt', sep='\t', header=True, index=True)

plt.figure(figsize=(12,12))
sns.scatterplot(data=meta,x='new_x',y='new_y',hue='max_celltype',palette='tab20',s=2,linewidth=0)
plt.gca().set_aspect('equal')
plt.legend(bbox_to_anchor=(1.01, 0), loc=3, borderaxespad=0)
plt.savefig(sys.argv[3]+'.max_celltype.pdf', dpi=300, bbox_inches='tight', format='pdf')

