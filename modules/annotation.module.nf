

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

    publishDir "${params.output}/${row.out_dir}/snpeff", mode: 'copy'

    input:
        tuple val(name), \
              val(row), \
              path("${name}/sequences.fa"), \
              path("${name}/genes.gtf.gz"), \
              path("snpeff_config_base.txt")

    output:
        path("snpEff.config")
        path("${name}/snpEffectPredictor.bin")
        path("${name}/genes.gtf.gz")
        path("${name}/sequences.fa")

    """
        {
            cat snpeff_config_base.txt;
            echo "${name}.genome : ${row.species}";
            has_mtdna=\$(grep 'MtDNA' ${name}/sequences.fa | wc -l)
            if (( \${has_mtdna} > 0 )); then
                echo "${name}.MtDNA.codonTable : Invertebrate_Mitochondrial"
            fi;
        } > snpEff.config
        snpEff build -c snpEff.config \
                     -dataDir . \
                     -gtf22 \
                     -v ${name}
    """

}


process format_csq {
    /*
        Generate a GFF file for CSQ annotation by BCFTools
    */

    publishDir "${params.output}/${row.out_dir}/csq", mode: 'copy'

    input:
        tuple val(row), path("in.gff.gz")

    output:
        path("${row.name}.csq.gff3.gz")
        path("${row.name}.csq.gff3.gz.tbi")

    shell:
    '''
        # to prep the gff3 for bcftools csq
        gzip -dc in.gff.gz | \
            awk '$2 ~ "WormBase.*"' | \
            sed -e 's/ID=Transcript:/ID=transcript:/g' \
                -e 's/ID=Gene:/ID=gene:/g' \
                -e 's/Parent=Transcript:/Parent=transcript:/g' \
                -e 's/Parent=Gene:/Parent=gene:/g' \
                -e 's/Parent=Pseudogene:/Parent=transcript:/g' > prep.gff
        
        format_csq.R

        {
            gzip -dc in.gff.gz | grep '^##';
            bedtools sort -i out.gff3;
        } > !{row.name}.csq.gff3
        
        bgzip !{row.name}.csq.gff3
        tabix -p gff !{row.name}.csq.gff3.gz
    '''

}