
process download_url {

    executor 'local'
    container null

    tag { "${row.species}/${row.project}" }

    input:
        tuple val(row), val(url)
    
    output:
        tuple val(row), path("genome.fa.gz")

    """
    wget ${url} -O genome.fa.gz
    """
}

process gzip_to_bgzip {

    label = "sm"
    publishDir "${row.out_dir}", mode: 'copy', pattern: "*fa.gz"

    input:
        tuple val(row), path(genome_in)

    output:
        tuple val(row.name), path("${row.genome}.fa"), emit: uncompressed
        tuple val(row), path("${row.genome}.fa.gz"), emit: compressed

    """
    gzip -dc ${genome_in} > ${row.genome}.fa
    bgzip ${row.genome}.fa -c > ${row.genome}.fa.gz # part of samtools
    """
}

process bwa_index {

    label = "sm"
    publishDir "${row.out_dir}", mode: 'copy'

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

    label = "sm"
    publishDir "${row.out_dir}", mode: 'copy'

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

    label = "sm"
    publishDir "${row.out_dir}", mode: 'copy'

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