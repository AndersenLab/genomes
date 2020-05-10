# genomes

Scripts and tools for managing reference genomes

## Overview

This repo contains a nextflow pipeline that downloads, indexes, and builds annotation databases for reference genomes from wormbase.


1. Downloads genomes from wormbase specified as `species/project` for a given `build`

## Usage

The pipeline can be run locally or on Quest. For example:

```bash
nextflow run main.nf -resume -profile local --wb_version=WS276 --projects=c_elegans/PRJNA13758
```

### `-profile`

Can be set to `local` or `quest`. The pipeline uses the `andersenlab/genomes` docker image built from [`env/genome.Dockerfile`](env/genome.Dockerfile). The image is automatically built on github using actions.

### `-wb_version`

The wormbase version to build.

### `--projects`

A comma-delimited list of `species/project_id` identifiers. To see a list of available projects that can be downloaded run the pipeline without specifying any projects. The first step of the pipeline downloads the list of available projects for all species on wormbase.

```bash
nextflow run main.nf --projects=""
```

This will generate a file called `project_species.tsv` in the `params.output` folder (`./genomes` if working locally).

