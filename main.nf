#!/usr/bin/env nextflow 
/*
Andersen Lab Genome preparation pipeline

This pipeline is designed to download genomic datasets
and from wormbase and process them for downstream analysis.

Author: Daniel E. Cook
*/
nextflow.preview.dsl=2
//assert System.getenv("NXF_VER") == "20.01.0-rc1"
assert nextflow.version.matches('20.0+')

// Constants
WORMBASE_PREFIX = "ftp://ftp.wormbase.org/pub/wormbase/releases"

/*
    Params
*/
//params.output="genomes"
params.snpeff_config = "${workflow.projectDir}/data/snpeff_config_base.txt"
params.genome = null // set for manual genome not from wormbase
params.gff = null // set for manual genome not from wormbase
params.output = "${workflow.projectDir}/genomes-test/"

if(!params.genome) {
    params.wb_version="WS276"
    params.projects="""c_elegans/PRJNA13758,c_briggsae/PRJNA10731,c_tropicalis/PRJNA53597"""
    project_list = params.projects.split(",")
} else {
    params.wb_version=null
    params.projects=null
    project_list = null
}


/* Includes */
// Downloads
include { download_url as download_genome;
          download_url as download_gtf;
          download_url as download_gff3; } from './modules/download.module.nf'

// Genome
include gzip_to_bgzip from './modules/genome.module.nf'
include bwa_index from './modules/genome.module.nf'
include samtools_faidx from './modules/genome.module.nf'
include create_sequence_dictionary from './modules/genome.module.nf'

// Annotation
include snpeff_db from './modules/annotation.module.nf'
include snpeff_db_manual from './modules/annotation.module.nf'
include { decompress as decompress_genome; } from './modules/annotation.module.nf'
include format_csq from './modules/annotation.module.nf'
include format_csq_manual from './modules/annotation.module.nf'
include extract_lcrs from './modules/annotation.module.nf'
include convert_gff_to_gtf from './modules/annotation.module.nf'

def log_summary() {
/*
    Generates a log
*/

out = '''
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
    --output                Path of output folder                          ${params.output}

    username                                                            ${"whoami".execute().in.text}

"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}

log.info "Genome Management"


workflow {
    // Download genome if wsbuild OR use provided genome and gff to build annotation files
    if(!params.genome) {
        println "Downloading ${params.wb_version} --> ${project_list}"

        fetch_projects()
        genome_set = fetch_projects.out.splitCsv(header: true, sep: "\t")
                    .filter { project_list.contains("${it.species}/${it.project}") }
                    .map { row ->
                        // Create output directory stub
                        row.name = "${row.species}.${row.project}.${params.wb_version}"
                        row.genome = "${row.name}.genome";
                        row.out_dir = "${row.species}/genomes/${row.project}/${params.wb_version}";
                        row;
                    }

        // Download genome and index
        format_dl(genome_set, "genomic.fa.gz") | download_genome | gzip_to_bgzip | (bwa_index & samtools_faidx & create_sequence_dictionary)

        // Download
        format_dl(genome_set, "canonical_geneset.gtf.gz") | download_gtf
        format_dl(genome_set, "annotations.gff3.gz") | download_gff3

        /* SnpEff */
        genome_eff = decompress_genome(download_genome.out).map { row, genome -> [row.name, row, genome] }
        gtf_eff = download_gtf.out.map { row, gtf -> [row.name, gtf] }
        genome_eff.join(gtf_eff)
                .combine(Channel.fromPath(params.snpeff_config)) | snpeff_db

        /* CSQ Annotations */
        download_gff3.out | format_csq

        /* Extract LCRs and other annotations */
        download_gff3.out | extract_lcrs
    } else {
        myFile = file("${params.genome}")
        println("Managing genome: ${myFile.getBaseName()}")

        /* SnpEff */
        Channel.from("${myFile.getBaseName()}".replaceFirst(/.genome.fa/, "")) // name
            .combine(Channel.from("${myFile.getParent()}")) //outdir
            .combine(Channel.fromPath("${params.genome}")) // genome
            .combine(Channel.fromPath("${params.gff}"))// gff 
            .combine(Channel.fromPath(params.snpeff_config)) | convert_gff_to_gtf | snpeff_db_manual 

        /* CSQ Annotations */
        Channel.from("${myFile.getBaseName()}".replaceFirst(/.genome.fa/, ""))
            .combine(Channel.from("${myFile.getParent()}"))
            .combine(Channel.fromPath("${params.gff}")) | format_csq_manual

        /* Extract LCRs and other annotations */
        //
    }

    


}

process fetch_projects {

    executor 'local'

    publishDir "${params.output}/", mode: 'copy'

    output:
        path("project_species.tsv")

    shell:
    '''
    set -e
    wormbase_url="ftp://ftp.wormbase.org/pub/wormbase"
    species=$(curl ${wormbase_url}/species/ | awk '{ print $9 }' | grep -v "README")

    function fetch_projects {
        species=${1}
        curl "${wormbase_url}/species/${species}/" | \
            awk -v OFS="\t" -v species=${species} '$9 ~ "^PR" { print species, $9 }'
    }
    export wormbase_url
    export -f fetch_projects
    {
        echo -e "species\tproject";
        parallel --retry-failed -j 8 --verbose fetch_projects ::: ${species};
    } > project_species.tsv
    '''
}

def format_dl(ch, fname) {
    ch.map { row ->
    [row,
     "${WORMBASE_PREFIX}/${params.wb_version}/species/${row.species}/${row.project}/${row.species}.${row.project}.${params.wb_version}.${fname}"]
    }
}

def download(ch, fname) {
    ch.map { row ->
        [row, 
         format_url(row, fname)]
    }
}

