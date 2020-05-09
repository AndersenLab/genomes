

process gzip_to_bgzip {

    publishDir "${params.output}/${row.out_dir}", mode: 'copy'

    input:
        tuple val(row), path("genome_in.fa.gz")

    output:
        tuple val(row), path("${row.genome}.fa.gz")

    """
        gzip -dc genome_in.fa.gz > genome.fa
        bgzip genome.fa
        mv genome.fa.gz ${row.genome}.fa.gz
    """

}

process bwa_index {

    publishDir "${params.output}/${row.out_dir}", mode: 'copy'

    input:
        tuple val(row), path("${row.genome}.fa.gz")

    output:
        path("${row.genome}.fa.gz.amb")
        path("${row.genome}.fa.gz.ann")
        path("${row.genome}.fa.gz.bwt")
        path("${row.genome}.fa.gz.pac")
        path("${row.genome}.fa.gz.sa")

    """
        bwa index ${row.genome}.fa.gz
    """
}

process samtools_faidx {

    publishDir "${params.output}/${row.out_dir}", mode: 'copy'

    input:
        tuple val(row), path("${row.genome}.fa.gz")

    output:
        path("${row.genome}.fa.gz.fai")
        path("${row.genome}.fa.gz.gzi")

    """
        samtools faidx ${row.genome}.fa.gz
    """
}

process create_sequence_dictionary {

    publishDir "${params.output}/${row.out_dir}", mode: 'copy'

    input:
        tuple val(row), path("${row.genome}.fa.gz")

    output:
        path("${row.genome}.dict")
        
    """
    picard CreateSequenceDictionary \
        R=${row.genome}.fa.gz \
        O=${row.genome}.dict
    """
}