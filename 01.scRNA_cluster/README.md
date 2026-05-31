## For each dataset ( species+organ ), we annotated cell type via :

* QC each library using SOUPX and DoubletFinder

See 01.scRNA_qc.r for detailed scripts

* Removing batch effect by Harmony and clustering by Seurat FindCluster

See 02.harmony.R for detailed scripts  