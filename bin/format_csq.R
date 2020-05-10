#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(vroom)

# CSQ Annotations require a GFF where the biotype
# of each gene *and* its transcripts are labeled 
# in the 9th column of the GFF [attributes]
# Wormbase provided GFF files only label the
# parent (gene) record, and not child transcript 
# records. This script copies the gene biotype 
# and appends it to the child transcript records.

# change upper case to lower case in ID and Parent
WS_gtf <- vroom::vroom("prep.gff", comment="#", delim="\t", col_names=FALSE) 

# For new genomes, a biotype may not have been assigned.
# In these instances, use the one assigned in the notes field (e.g. Note=PREDICTED protein_coding)
# And assign this as the biotype as long as no other biotype field is detected.
WS_gtf <- WS_gtf %>% 
            dplyr::mutate(biotype_from_note = stringr::str_match(X9, "Note=PREDICTED ([^;]+)")[,2]) %>%
            dplyr::mutate(X9 = ifelse(!is.na(biotype_from_note) & !grepl("biotype", X9),
                                      glue::glue("{X9};biotype={biotype_from_note}"),
                                      X9)) %>%
            dplyr::select(-biotype_from_note)

gene_biotype_table <- WS_gtf %>% 
                        dplyr::filter(X3 == "gene") %>% 
                        dplyr::mutate(Name = stringr::str_match(X9, "ID=gene:([^;]+)")[,2]) %>%
                        dplyr::mutate(biotype = stringr::str_match(X9, "biotype=([^;]+)")[,2]) %>%
                        dplyr::select(Name, biotype)

WS_gtf <- WS_gtf %>%
    dplyr::mutate(Name = stringr::str_match(X9, "Parent=gene:([^;]+)")[,2])

# Merge biotypes on transcript rows
WS_gtf <- dplyr::left_join(WS_gtf, gene_biotype_table) %>%
          dplyr::mutate(X9 = ifelse(!is.na(biotype), glue::glue("{X9};biotype={biotype}"), X9)) %>%
          dplyr::select(-Name, -biotype)

# Sort
WS_gtf <- WS_gtf %>% dplyr::arrange(X1, X4, X5)

vroom::vroom_write(WS_gtf, path = "out.gff3", col_names = FALSE)
