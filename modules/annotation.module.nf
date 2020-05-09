

process decompress {

    input:
        tuple val(row), path("file.gz")
    
    output:
        tuple val(row), path("file_decompressed")

    """
        gzip -dc file.gz > file_decompressed
    """
}

process snpeff_db {

    tag { "${row.name}" }

    publishDir "${params.output}/${row.out_dir}", mode: 'copy'

    input:
        tuple val(name), \
              val(row), \
              path("${name}/sequences.fa"), \
              path("${name}/protein.fa"), \
              path("${name}/cds.fa"), \
              path("${name}/genes.gff.gz"), \
              path("snpeff_config.txt")

    """
        #
        {
            cat snpeff_config.txt;
            echo "${name}.genome : ${row.species}";
            has_mtdna=\$(gzip -dc ${name}/genome.fa.gz | grep 'MtDNA' | wc -l)
            if (( \${has_mtdna} > 0 )); then
                echo "${name}.MtDNA.codonTable : Invertebrate_Mitochondrial"
            fi;
        } > snpEff.config
        snpEff build -c snpEff.config \
                     -dataDir . \
                     -gff3 \
                     -v ${name}
    """

}
