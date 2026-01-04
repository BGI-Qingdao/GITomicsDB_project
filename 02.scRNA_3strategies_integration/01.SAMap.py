# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 08:06:34 2022

@author: Ray
"""

import argparse as arg
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                                sankey_plot, chord_plot, CellTypeTriangles,
                                ParalogSubstitutions, FunctionalEnrichment,
                                convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pickle
import os

parser = arg.ArgumentParser(formatter_class=arg.RawDescriptionHelpFormatter, description = '''This programe is used to perform SAMap.''')
parser.add_argument('-mappath', action='store', help='Input blast path.')
parser.add_argument('-species1', action='store', help='Input expression list.')
parser.add_argument('-species2', action='store', help='Output expression matrix.')
parser.add_argument('-fn1', action='store', help='Species h5ad')
parser.add_argument('-fn2', action='store', help='Species h5ad')
parser.add_argument('-celltype1', action='store', help='Celltype1.')
parser.add_argument('-celltype2', action='store', help='Celltype2.')
parser.add_argument('-NUMITERS', type=int, default=3, action='store', help='NUMITERS.')
parser.add_argument('-cpu', type=int, default=20, help='Number of cpu used.')
parser.add_argument('-mappingtable', action='store', help='Mappingtable path.')
parser.add_argument('-enrichedGene', action='store', help='Enriched genes table')
args = parser.parse_args()

fname = args.species1+args.species2+'.sm'

filenames = {args.species1:args.fn1, args.species2:args.fn2}
keys={args.species1:args.celltype1, args.species2:args.celltype2}
if not os.path.exists(fname):
    print(args.mappath)
    sm = SAMAP(
        filenames,
        keys=keys,
        f_maps = args.mappath,
        save_processed=False,
       )
    sm.run(NUMITERS=args.NUMITERS, ncpus=args.cpu)
    #with open(args.species1+args.species2+'.sm', 'wb') as f:
    #    pickle.dump(sm, f)
else:
    f = open(fname, 'rb')
    sm = pickle.load(f)

D,MappingTable=get_mapping_scores(sm, keys, n_top=0)
MappingTable.to_csv(args.mappingtable, sep='\t')
