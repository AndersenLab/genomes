# first, start with classic elegans
nextflow run ../main.nf --projects c_elegans/PRJNA13758 --wb_version WS276 --output elegans_test

# briggsae from wormbase
nextflow run ../main.nf --projects c_briggsae/PRJNA10731 --wb_version WS276 --output briggsae_WB_test

# tropicalis manual
nextflow run ../main.nf --genome /projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa.gz \
--gff /projects/b1059/projects/Nicolas/c.tropicalis/NIC58/predictions/merger/NIC58.final_annotation.renamed.csq.gff --species c_tropicalis --projects NIC58_nanopre \
--wb_version June2021 --output tropicalis_man_test