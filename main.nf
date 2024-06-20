#!/usr/bin/env nextflow 
/*
Andersen Lab Genome preparation pipeline

This pipeline is designed to download genomic datasets
and from wormbase and process them for downstream analysis.

Authors:
    Daniel E. Cook
    Mike Sauria
*/

nextflow.enable.dsl=2
//assert System.getenv("NXF_VER") >= "23.0"
assert nextflow.version.matches('23.0+')

/* Includes */
// Downloads
include { download_url as download_genome;
          download_url as download_gtf;
          download_url as download_gff3; } from './modules/genome.module.nf'

// Genome
include { gzip_to_bgzip } from './modules/genome.module.nf'
include { bwa_index } from './modules/genome.module.nf'
include { samtools_faidx } from './modules/genome.module.nf'
include { create_sequence_dictionary } from './modules/genome.module.nf'

// Annotation
include { snpeff_db } from './modules/annotation.module.nf'
include { snpeff_db_manual } from './modules/annotation.module.nf'
include { format_csq } from './modules/annotation.module.nf'
include { format_csq_manual } from './modules/annotation.module.nf'
include { extract_lcrs } from './modules/annotation.module.nf'

// Constants
WORMBASE_PREFIX = "https://downloads.wormbase.org"

/*
    Params
*/
params.snpeff_config = "${workflow.projectDir}/data/snpeff_config_base.txt"
params.genome = null // set for manual genome not from wormbase
params.gff = null // set for manual genome not from wormbase

if(!params.output) {
    params.outputDir = "${workflow.launchDir}/genomes_test"
} else {
    params.outputDir = params.output
}


if(!params.genome) {
    params.wb_version="WS290"
    params.projects="""c_elegans/PRJNA13758,c_briggsae/PRJNA10731,c_tropicalis/PRJNA53597"""
    project_list = params.projects.split(",")
    params.project = null
} else {
    params.wb_version=null
    params.project=null
    project_list = null
    if (!params.project || !params.wb_version) {
        log.info("If a genome parameter is specified, project and wb_version parameters must also be specified")
        exit 1
    }
    params.species=params.genome.split('/')[-1]
    params.projects="${params.species}/${params.project}"
}

if (params.help) {
    log.info('''
>AAGACGACTAGAGGGGGCTATCGACTACGAAACTCGACTAGCTCAGCGGGATCAGCATCACGATGGGGGCCTATCTACGACAAAATCAGCTACGAAA
AGACCATCTATCATAAAAAATATATATCTCTTTCTAGCGACGATAAACTCTCTTTCATAAATCTCGGGATCTAGCTATCGCTATATATATATATATGC
GAAATA      CGCG       GA ATATA AAAA    TCG TCGAT GC       GGGC     CGATCGA TAGAT GA      TATATCGC
TTAAC ACTAGAGGGG CTATCGAC  CGAA CT GACTA CT  GCG  AT AGCATCACG TGGGGGCCTATC  CGAC AA TCAGCTACGAAAT
AGCCC TCTATCATAA    TATAT T TCT TC AGCGA GA A A T TC    ATAAAT TCGGGATCTAGC A CGC AT    ATATATATGC
GCGAT TCTAC   AG GCGGGGGA AT TA AA AAGAC CG TC AT GC AGCTGGGGGC    ACG   GA TA AT GA CTATATATATCGC
AATGC ACTAGAG GG CTATCGAC ACG A CT GACTA CT AGCGG AT AGCATCACGATGGG GCCTATC ACG C AA TCAGCTACGAAAT
ACTCC TCTATCA AA AAATATAT TCTC  TC AGCGA GA AAACT TC TTCATAAATCTCGG ATCTAGC ATCG  AT TATATATATATGC
TTAATA       FCG       GA ATATA AAA     TCG TCGAT GC        GG     ACGATCGA TAGAT GA CTATATATATCGC
AACACGACTAGAGGGGGCTATCGACTACGAAACTCGACTAGCTCAGCGGGATCAGCATCACGATGGGGGCCTATCTACGACAAAATCAGCTACGAAAT
CTACCATCTATCATAAAAAATATATATCTCTTTCTAGCGACGATAAACTCTCTTTCATAAATCTCGGGATCTAGCTATCGCTATATATATATATATGC

''' + """
To run the pipeline:

nextflow main.nf --projects c_elegans/PRJNA13758 --wb_version WS276

    parameters              description                                    Set/Default
    ==========              ===========                                    ========================
    --wb_version            wormbase version to build                      ${params.wb_version}
    --projects              comma-delimited list of `species/project_id`   ${params.projects}
    --genome                Path to manually curated genome                ${params.genome}
    --output                Path of output folder                          ${params.outputDir}

    username                                                            ${"whoami".execute().in.text}

""")
    exit 1
} else {
    log.info("""
G E N O M E S - N F   P I P E L I N E
========================================
wb_version =            ${params.wb_version}
projects   =            ${params.projects}
genome     =            ${params.genome}
output     =            ${params.outputDir}
    """)
}

log.info "Genome Management"


workflow {
    // Download genome if wsbuild OR use provided genome and gff to build annotation files
    if(!params.genome) {
        println "Downloading ${params.wb_version} --> ${project_list}"

        Channel.of("${WORMBASE_PREFIX}") | fetch_projects

        genome_set = fetch_projects.out.splitCsv(header: true, sep: "\t")
                    .filter { project_list.contains("${it.species}/${it.project}") }
                    .map { row ->
                        // Create output directory stub
                        row.name = "${row.species}.${row.project}.${params.wb_version}"
                        row.genome = "${row.name}.genomic";
                        row.out_dir = "${params.outputDir}/${row.species}/genomes/${row.project}/${params.wb_version}";
                        row;
                    }

        // Download genome and index
        genome_set.map { row -> [row, "${WORMBASE_PREFIX}/releases/${params.wb_version}/species/${row.species}/${row.project}/${row.species}.${row.project}.${params.wb_version}.genomic.fa.gz"] }
            | download_genome | gzip_to_bgzip
        gzip_to_bgzip.out.compressed | (samtools_faidx & bwa_index & create_sequence_dictionary)

        // Download
        genome_set.map { row -> [row, "${WORMBASE_PREFIX}/releases/${params.wb_version}/species/${row.species}/${row.project}/${row.species}.${row.project}.${params.wb_version}.canonical_geneset.gtf.gz"] }
            | download_gtf
        genome_set.map { row -> [row, "${WORMBASE_PREFIX}/releases/${params.wb_version}/species/${row.species}/${row.project}/${row.species}.${row.project}.${params.wb_version}.annotations.gff3.gz"] }
            | download_gff3

        /* SnpEff */
        genome_eff = download_genome.out
            .map { row, genome -> [row.name, row] }
            .join(gzip_to_bgzip.out.uncompressed)
            .join(download_gtf.out
                .map { row, gtf -> [row.name, gtf] })
            .combine(Channel.fromPath(params.snpeff_config)) | snpeff_db

        /* CSQ Annotations */
        download_gff3.out | format_csq

        /* Extract LCRs and other annotations */
        download_gff3.out | extract_lcrs

       // get gene and transcript tracks for cendr (elegans only right now)
        download_gff3.out | cendr_browser_tracks
    } else {
        myFile = file("${params.genome}")
        println("Managing genome: ${myFile.getBaseName()}")

        // Setup genome files
       Channel.of("${params.species}")
            .combine(Channel.of("${params.project}")) | manual_setup

        genome_set = manual_setup.out.splitCsv(header: true, sep: "\t")
            .map { row ->
                // Create output directory stub
                row.name = "${row.species}.${row.project}.${params.wb_version}"
                row.genome = "${row.name}.genomic";
                row.out_dir = "${params.outputDir}/${row.species}/genomes/${row.project}/${params.wb_version}";
                row; }

        genome_set.combine(Channel.fromPath("${params.genome}")) | gzip_to_bgzip
        gzip_to_bgzip.out.compressed | (bwa_index & samtools_faidx & create_sequence_dictionary)


        if (params.gff) {
            /* SnpEff */
            genome_set
                .map{ row -> [row.name, row] }
                .join(gzip_to_bgzip.out.uncompressed)
                .combine(Channel.fromPath("${params.gff}"))// gff 
                .combine(Channel.fromPath(params.snpeff_config)) | snpeff_db_manual 
            
            /* CSQ Annotations */
            genome_set
                .combine(Channel.fromPath("${params.gff}")) | format_csq_manual
        }

        /* Extract LCRs and other annotations */
        // we don't have this!
    }

    


}

process fetch_projects {

    executor 'local'
    container null

    publishDir "${params.outputDir}/", mode: 'copy'

    input:
        val(WORMBASE_PREFIX)

    output:
        path("project_species.tsv")

    shell:
    """
    set -e
    species=(`wget -q -O - ${WORMBASE_PREFIX}/species/ | grep folder | \\
             awk '{split(\$0,A,"href=\\""); split(A[2], B, "/"); print B[1]}'`)

    echo -e "species\\tproject" > project_species.tsv
    for S in \${species[*]}; do
        wget -q -O - ${WORMBASE_PREFIX}/species/\${S}/ | grep folder | \\
            awk -v OFS="\\t" -v species=\${S} '{split(\$0,A,"href=\\""); split(A[2], B, "/"); \\
                                                if (( B[1] ~ /^PR/ )) print species,B[1]}' >> project_species.tsv
    done
    """
}

process manual_setup {

    executor "local"
    container null

    input:
        tuple val(species), val(project)

    output:
        file("setup_file.txt")

    """
    echo -e "species\tproject" > setup_file.txt
    echo "${species}\t${project}" >> setup_file.txt
    """
}



// download big bed file for gene and transcript tracks for CeNDR
// looks like this file isonly currently supported for elegans and briggsae but not tropicalis.
process cendr_browser_tracks {

    label "sm"

    publishDir "${row.out_dir}/browser_tracks/", mode: 'copy'

    input:
        tuple val(row), path("gff")

    output:
        tuple path("*.bed.gz"), path("*.bed.gz.tbi"), optional: true

    """
    if [[ ${row.species} != 'c_tropicalis' ]];
    then

        function zip_index {
            bgzip -f \${1}
            tabix \${1}.gz
        }

        species=`echo "${row.species}" | sed "s/c_//"`

        # Generate the transcripts track; Confusingly, this track is derived from one called elegans_genes on wormbase.
        # Add parenthetical gene name for transcripts.

        wget -O gene_file.bb ${WORMBASE_PREFIX}/releases/${params.wb_version}/MULTI_SPECIES/hub/\${species}/\${species}_genes_${params.wb_version}.bb
        bigBedToBed gene_file.bb tmp.bed
        sortBed -i tmp.bed > \${species}_transcripts_${params.wb_version}.bed
        bgzip -f \${species}_transcripts_${params.wb_version}.bed
        tabix \${species}_transcripts_${params.wb_version}.bed.gz
        rm tmp.bed

        # Generate Gene Track BED File 
        gzip -dc ${gff} | \
        grep 'locus' | \
        awk '\$2 == "WormBase" && \$3 == "gene"' > temp.gff
        sortBed -i temp.gff > sorted.gff

        # Install with conda install gawk
        convert2bed --sort-tmpdir=./ -i gff - < sorted.gff > final.bed
        gawk -v OFS='\\t' '{ match(\$0, "locus=([^;\\t]+)", f); \$4=f[1]; print \$1, \$2, \$3, \$4, 100, \$6  }' final.bed | \
        uniq > \${species}_gene.${params.wb_version}.bed
        zip_index \${species}_gene.${params.wb_version}.bed

    fi
    """
}



