

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
        path("${row.name}.AA_Scores.tsv")
        path("${row.name}.AA_Length.tsv")

    """
        # Ryan's fix for formatting gff3 for bcsq
        zcat in.gff.gz | grep -P "\tWormBase\t" > simple.wormbase.gff3

        Rscript --vanilla ${workflow.projectDir}/bin/bcsq_gff_format.R

        cp fixed.wormbase.gff3 ${row.name}.csq.gff3
        bgzip ${row.name}.csq.gff3
        tabix -p gff ${row.name}.csq.gff3.gz

        # also run AA_scores and AA_length
        Rscript --vanilla ${workflow.projectDir}/bin/AA_Scores_Table.R ${workflow.projectDir}/bin/BLOSUM62
        mv AA_Scores.tsv ${row.name}.AA_Scores.tsv

        Rscript --vanilla ${workflow.projectDir}/bin/AA_Length.R ${row.name}.csq.gff3.gz
        mv gff_AA_Length.tsv ${row.name}.AA_Length.tsv

    """

}

process extract_lcrs {
    /*
        Extract low complexity regions for:

        dust
        RepeatMasker
    */

    publishDir "${params.output}/${row.out_dir}/lcr", mode: 'copy'

    input:
        tuple val(row), path("in.gff.gz")

    output:
        path("${row.name}.dust.bed.gz")
        path("${row.name}.dust.bed.gz.tbi")
        path("${row.name}.repeat_masker.bed.gz")
        path("${row.name}.repeat_masker.bed.gz.tbi")

    shell:
    '''
        gzip -dc in.gff.gz | \
        awk -v OFS="\t" '$2 == "dust" { print $1, $4, $5, $2 > "!{row.name}.dust.bed" } 
                         $2 == "RepeatMasker" {  print $1, $4, $5, $2 > "!{row.name}.repeat_masker.bed" }'
        bgzip !{row.name}.dust.bed
        bgzip !{row.name}.repeat_masker.bed
        tabix -f -p bed !{row.name}.dust.bed.gz
        tabix -f -p bed !{row.name}.repeat_masker.bed.gz
    '''
}