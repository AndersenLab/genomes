

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

// aa length currently doesn't work unless c.e. need to potentially change gff for other species
process format_csq {

    tag { "${row.name}" }

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
        path("*.gff")

    """
    if [[ ${row.species} = "c_elegans" ]]
    then
        # Ryan's fix for formatting gff3 for bcsq
        zcat in.gff.gz | grep -P "\tWormBase\t" > simple.wormbase.gff3

        Rscript --vanilla ${workflow.projectDir}/bin/bcsq_gff_format.R

        cp fixed.wormbase.gff3 ${row.name}.csq.gff3
        bgzip ${row.name}.csq.gff3
        tabix -p gff ${row.name}.csq.gff3.gz

    else
        # to prep the gff3 for bcftools csq
        gzip -dc in.gff.gz | \
            awk '\$2 ~ "WormBase.*"' | \
            sed -e 's/ID=Transcript:/ID=transcript:/g' \
                -e 's/ID=Gene:/ID=gene:/g' \
                -e 's/Parent=Transcript:/Parent=transcript:/g' \
                -e 's/Parent=Gene:/Parent=gene:/g' \
                -e 's/Parent=Pseudogene:/Parent=transcript:/g' > prep.gff
        
        format_csq.R
        {
            # gzip -dc in.gff.gz | grep '^##'; # I don't understand this line
            bedtools sort -i out.gff3;
        } > ${row.name}.csq.gff3
        
        bgzip ${row.name}.csq.gff3
        tabix -p gff ${row.name}.csq.gff3.gz

    fi

    # also run AA_scores and AA_length
    Rscript --vanilla ${workflow.projectDir}/bin/AA_Length.R ${row.name}.csq.gff3.gz
    mv gff_AA_Length.tsv ${row.name}.AA_Length.tsv

    Rscript --vanilla ${workflow.projectDir}/bin/AA_Scores_Table.R ${workflow.projectDir}/bin/BLOSUM62
    mv AA_Scores.tsv ${row.name}.AA_Scores.tsv

    # also get gene file for nemascan
    Rscript --vanilla ${workflow.projectDir}/bin/gene_file_nemascan.R ${row.name}.csq.gff3.gz ${row.species}

    """

}

// aa length currently doesn't work unless c.e. need to potentially change gff for other species
process format_csq_manual {

    /*
        Generate a GFF file for CSQ annotation by BCFTools
    */

    publishDir "${params.output}/${row.out_dir}/csq", mode: 'copy'
    // publishDir "${params.output}/csq", mode: 'copy' // use for debug test

    input:
        tuple val("row"), path("gff_file")

    output:
        path("${row.name}.csq.gff3.gz")
        // path("${row.name}.csq.gff3.gz.tbi")
        // path("${row.name}.csq.gff3")
        // path("${name}.csq.gff3.tbi")
        path("${row.name}.AA_Scores.tsv")
        path("${row.name}.AA_Length.tsv")
        path("*.gff")

    """
    bedtools sort -i ${gff_file} > ${row.name}.csq.gff3
    bgzip ${row.name}.csq.gff3
    tabix -p gff ${row.name}.csq.gff3.gz

    Rscript --vanilla ${workflow.projectDir}/bin/AA_Length.R ${row.name}.csq.gff3.gz
    mv gff_AA_Length.tsv ${row.name}.AA_Length.tsv
    
    # also run AA_scores and AA_length
    Rscript --vanilla ${workflow.projectDir}/bin/AA_Scores_Table.R ${workflow.projectDir}/bin/BLOSUM62
    mv AA_Scores.tsv ${row.name}.AA_Scores.tsv

    # also get gene file for nemascan
    Rscript --vanilla ${workflow.projectDir}/bin/gene_file_nemascan.R ${row.name}.csq.gff3.gz ${row.species}

    """

}


// use gff instead of gtf to create snpeff config manually
process snpeff_db_manual {

    publishDir "${params.output}/${row.out_dir}/snpeff", mode: 'copy'
    // publishDir "${params.output}/snpeff", mode: 'copy' // use for debug test

    input:
        tuple val(name), val(row), \
              path("${name}/sequences.fa"), \
              path("${name}/genes.gff"), \
              path("snpeff_config_base.txt")

    output:
        path("snpEff.config")
        path("${name}/snpEffectPredictor.bin")
        path("${name}/genes.gff")
        path("${name}/sequences.fa")

    """
    sp=`echo ${name} | cut -d '.' -f 1`

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
                     -gff3 \
                     -v ${name}
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

