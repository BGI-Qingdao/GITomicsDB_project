# -*- coding: utf-8 -*-

import argparse as arg
from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores, GenePairFinder
import scanpy as sc
import pandas as pd
parser = arg.ArgumentParser()
parser.add_argument('-mappath')
parser.add_argument('-species1')
parser.add_argument('-species2')
parser.add_argument('-fn1')
parser.add_argument('-fn2')
parser.add_argument('-celltype1')
parser.add_argument('-celltype2')
parser.add_argument('-NUMITERS', type=int, default=3)
parser.add_argument('-cpu', type=int, default=5)
parser.add_argument('-mappingtable')
parser.add_argument('-enrichedGene')
args = parser.parse_args()
# read h5ad
adata1 = sc.read_h5ad(args.fn1)
adata2 = sc.read_h5ad(args.fn2)
# ensure unique gene names
adata1.var_names_make_unique()
adata2.var_names_make_unique()
# species files
filenames = {args.species1: args.fn1,args.species2: args.fn2}
# celltype keys
keys = {args.species1: args.celltype1,args.species2: args.celltype2}
# run SAMAP
sm = SAMAP(filenames,keys=keys,f_maps=args.mappath,save_processed=False)
sm.run(NUMITERS=args.NUMITERS,ncpus=args.cpu)
# mapping scores
D, MappingTable = get_mapping_scores(sm,keys,n_top=0)
MappingTable.to_csv(args.mappingtable,sep='\t')
# save UMAP
def save_umap(fname, data):
  df=pd.DataFrame({'samap_1': data[:,0],'samap_2': data[:,1]})
  df.to_csv(fname,sep='\t',index=False)
save_umap(f'{args.species1}.samap.csv',sm.sams[args.species1].adata.obsm['X_umap_samap'])
save_umap(f'{args.species2}.samap.csv',sm.sams[args.species2].adata.obsm['X_umap_samap'])
# enriched genes
gpf = GenePairFinder(sm,keys=keys)

gene_pairs = gpf.find_all(align_thr=0.1,thr=0.05)

gene_pairs.to_csv(args.enrichedGene,sep='\t',index=False)
