#!/usr/bin/env nextflow 
/*
Andersen Lab Genome preparation pipeline

This pipeline is designed to download genomic datasets
and from wormbase and process them for downstream analysis.

Author: Daniel E. Cook
*/
nextflow.preview.dsl=2
assert System.getenv("NXF_VER") == "20.01.0-rc1"

log.info "Genome Management"

/*
    Params
*/
params.output="genome"
params.build="WS276"
params.projects="PRJNA13758,PRJNA10731,PRJNA53597"
project_list = params.projects.split(",")

PROJECT_ALIAS = [PRJNA13758: "N2"]

process fetch_projects {

    publishDir "${params.output}/"

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



workflow {
    fetch_projects()
    fetch_projects.out.splitCsv(header: true, sep: "\t")
                      .filter { species, project -> project_list.Contains(project) }.view()
}