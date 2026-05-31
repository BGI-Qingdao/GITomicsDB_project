#!/bin/bash

A=$1
B=$2

source /dellfsqd2/ST_OCEAN/USER/guolidong/anaconda3_new/install/bin/activate saturn
python3 01.performSAMap.py \
    -mappath /dellfsqd2/ST_OCEAN/USER/guolidong/FishGut/11.sammap/ALL_BLAST/maps/ \
    -species1 $A \
    -species2 $B \
    -fn1 ./inh5ad/${A}.IS.lognorm.h5ad \
    -fn2 ./inh5ad/${B}.IS.lognorm.h5ad \
    -celltype1 celltypeIS \
    -celltype2 celltypeIS \
    -NUMITERS 3 \
    -cpu 30 \
    -mappingtable 02.${A}_${B}.mappingtable.txt\
    -enrichedGene 02.${A}_${B}.table.txt >${A}_${B}.log 2>&1
