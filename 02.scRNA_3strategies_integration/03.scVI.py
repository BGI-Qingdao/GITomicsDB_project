import anndata as ad
import numpy as np
import scvi
import scanpy as sc

adata = ad.read_h5ad('/ldfsqd1/ST_OCEAN/USER/guolidong/FishGut_V4/21.anno.H5ad/fish9hm.gene7700.h5ad')
adata = adata[adata.obs['group']!='Discard'].copy()
adata.obs['so'] = adata.obs.apply(lambda row:  f'{row["species"]}.{row["organ"]}',axis=1)
scvi.model.SCVI.setup_anndata(adata, batch_key="so")
vae = scvi.model.SCVI(adata)
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
np.savetxt('scVI.csv',adata.obsm["X_scVI"],delimiter=",")
adata.write('fish9hm.scvi.h5ad',compression='gzip')
sc.pp.neighbors(adata, use_rep="X_scVI")
for x in [0.1,0.2,0.3,0.5,0.8,1]:
    sc.tl.leiden(adata,resolution=x,key_added=f'leiden_{x}')
sc.tl.umap(adata)
adata.write('fish9hm.scvi.leiden.h5ad',compression='gzip')
