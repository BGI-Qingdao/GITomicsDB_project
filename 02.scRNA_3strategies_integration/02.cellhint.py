import scanpy as sc
import numpy as np
import pandas as pd
import anndata 
import cellhint
import seaborn as sns
import matplotlib.pyplot as plt
import sys
adata = sc.read_h5ad(sys.argv[1])
prefix= sys.argv[2]

sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key = 'organ', subset = True)

sc.pp.scale(adata, max_value = 10)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

alignment = cellhint.harmonize(adata, 'organ', 'celltype') #,dataset_order=order)
alignment.write(f'{prefix}.alignment.pkl')
#corr = alignment.base_distance.to_meta()
member_mat = alignment.base_distance.to_meta(turn_binary = True)
member_mat.to_csv(f'{prefix}.corr.csv', sep = '\t', header = True, index=True)
alignment.relation.to_csv(f'{prefix}.HT.csv', sep = ',', index = False)
sns.clustermap(member_mat)
plt.savefig(f'{prefix}.heatmap.pdf')
cellhint.treeplot(alignment,save=f'{prefix}.tree.pdf')
