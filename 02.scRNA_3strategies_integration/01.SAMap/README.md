## Scripts for run SAMap between 10 species

* step1:  build gene to gene network by BLASTP

See 00.build_gene2gene_of_LP.sh as example.

* step2: build SAMap cell type simularity matrix by run SAMap between each two species

See 01.batch_10_species_VS_10_species_SAMap.sh and the other two scripts called by it.

Why I run SAMap between each two species?  

Because all 10 species with near 1 millon cells will consume > 1 TB level memeory.

Plus, the result different between run 10 species at one VS run one by one is acceptable.