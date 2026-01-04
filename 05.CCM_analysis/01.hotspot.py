# -*- coding: utf-8 -*-
import scanpy as sc
import hotspot
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy
import pandas as pd
def load_gene_list(path_allmarkers):
    All_markers = pd.read_csv(path_allmarkers,sep = "\t",usecols = ['gene','count'])
    # All_markers = All_markers[All_markers['avg_log2FC'] > 0.2]
    All_markers = All_markers["gene"].tolist()
    return All_markers

def performHotSpot(adata,path_allmarkers,prefix):
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs["total_counts"] = np.asarray(adata.X.sum(1)).ravel()
    adata.layers["csc_counts"] = scipy.sparse.csc_matrix(adata.X)#.tocsc()
    #sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    custom_genes = load_gene_list(path_allmarkers)
    valid_genes = [g for g in custom_genes if g in adata.var_names]
    adata.var['custom_selected']=adata.var_names.isin(valid_genes)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    subadata = adata[:,adata.var['custom_selected']].copy()

    sc.pp.neighbors(subadata)
    sc.tl.umap(subadata)
    #Only use top 5000 genses to generate gene modules
    #subadata = adata[:,adata.var['highly_variable']]
    #Create the Hotspot object and the neighborhood graph
    hs = hotspot.Hotspot(
        subadata,
        layer_key="csc_counts",
        model='normal',
        latent_obsm_key="X_umap",
        umi_counts_obs_key="total_counts",
    )

    hs.create_knn_graph(
        weighted_graph=False,
        n_neighbors=10,
    )
    hs_results = hs.compute_autocorrelations(jobs=30)
    # Select the genes with significant spatial autocorrelation
    hs_genes = hs_results.index[hs_results.FDR < 0.05]
    # Compute pair-wise local correlations between these genes
    lcz = hs.compute_local_correlations(hs_genes, jobs=20)
    modules = hs.create_modules(
        min_gene_threshold=40, core_only=True, fdr_threshold=0.05
    )
    plt.rcParams.update({
        'pdf.fonttype':42,
        'ps.fonttype':42,
        'font.family':'Arial'
    })
    plt.figure()
    hs.plot_local_correlations(vmin=-5, vmax=5)
    plt.savefig(f'{prefix}.hotspot_correlations4.0.pdf')
    results = hs.results.join(hs.modules).sort_values('Z', ascending=False).to_csv(f'{prefix}.modules4.0.csv')
    module_scores=hs.calculate_module_scores()
    module_scores.to_csv(f'{prefix}.modulescore4.0.csv')

adata = sc.read(sys.argv[1])
#adata = adata[adata.obs['group'] == 'Epithelial'].copy()
performHotSpot(adata,sys.argv[2],sys.argv[3])
