# genomes

Scripts and tools for managing reference genomes

## Overview

This repo contains a nextflow pipeline that downloads, indexes, and builds annotation databases for reference genomes from wormbase. The following outputs are created:

1. A BWA Index
2. SNPeff annotation database
3. CSQ annotation database
4. Samtools faidx index
5. A GATK Sequence dictionary file

## Usage

The pipeline can be run locally or on Quest. For example:

```bash
nextflow run main.nf -resume -profile local --wb_version=WS276 --projects=c_elegans/PRJNA13758
```

By default, the pipeline will generate reference genome indices and annotations for:

* `c_elegans/PRJNA13758` - N2 based reference genome
* `c_briggsae/PRJNA10731`
* `c_tropicalis/PRJNA53597`

## Options

### `-profile`

Can be set to `local` or `quest`. The pipeline uses the `andersenlab/genomes` docker image built from [`env/genome.Dockerfile`](env/genome.Dockerfile). The image is automatically built using github actions. See [`.github/workflows/build.yml`](.github/workflows/build.yml) for details.

### `-wb_version`

The wormbase version to build. For example, `WS276`.

### `--projects`

A comma-delimited list of `species/project_id` identifiers. A table below lists the current projects that can be downloaded. This table is regenerated as the first step of the pipeline, and stored as a file called `project_species.tsv` in the `params.output` folder (`./genomes` if working locally).

The current set of available species/projects that can be built are:

| species         | project     |
|:----------------|:------------|
| b_xylophilus    | PRJEA64437  |
| c_briggsae      | PRJNA10731  |
| c_angaria       | PRJNA51225  |
| a_ceylanicum    | PRJNA231479 |
| a_suum          | PRJNA62057  |
| a_suum          | PRJNA80881  |
| b_malayi        | PRJNA10729  |
| c_brenneri      | PRJNA20035  |
| c_elegans       | PRJEB28388  |
| c_elegans       | PRJNA13758  |
| c_elegans       | PRJNA275000 |
| c_latens        | PRJNA248912 |
| c_remanei       | PRJNA248909 |
| c_remanei       | PRJNA248911 |
| c_remanei       | PRJNA53967  |
| c_inopinata     | PRJDB5687   |
| c_japonica      | PRJNA12591  |
| c_sp11          | PRJNA53597  |
| c_sp5           | PRJNA194557 |
| c_nigoni        | PRJNA384657 |
| c_sinica        | PRJNA194557 |
| c_tropicalis    | PRJNA53597  |
| d_immitis       | PRJEB1797   |
| h_bacteriophora | PRJNA13977  |
| l_loa           | PRJNA60051  |
| m_hapla         | PRJNA29083  |
| m_incognita     | PRJEA28837  |
| h_contortus     | PRJEB506    |
| h_contortus     | PRJNA205202 |
| n_americanus    | PRJNA72135  |
| p_exspectatus   | PRJEB6009   |
| o_tipulae       | PRJEB15512  |
| p_redivivus     | PRJNA186477 |
| s_ratti         | PRJEA62033  |
| s_ratti         | PRJEB125    |
| o_volvulus      | PRJEB513    |
| p_pacificus     | PRJNA12644  |
| t_muris         | PRJEB126    |
| t_spiralis      | PRJNA12603  |
| t_suis          | PRJNA208415 |
| t_suis          | PRJNA208416 |

