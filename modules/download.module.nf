

process download_url {

    tag { "${row.species}/${row.project}" }

    input:
        tuple val(row), val(url)
    
    output:
        tuple val(row), path("out.file")

    """
        wget -O out.file ${url}
    """

}