#!/usr/bin/env Rscript

#Libraries
library(dplyr)
library(stringr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

#Read in GFF
gff <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE)

# gff = read.delim("QX1410.csq.gff3", header=FALSE, stringsAsFactors=FALSE)
# gff = read.delim("c_elegans.PRJNA13758.WS276.csq.gff3.gz", header=FALSE, stringsAsFactors=FALSE)

names(gff) <- c("chrm_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

#Filter, Parse, Calculate DNA length
CDS <- dplyr::filter(gff, type=="CDS")
parsedCDS <- CDS %>%
  tidyr::separate("attributes", into = c("ID", "Parent", "Name", "Prediction_Status", "Wormprep", "Protein_ID", "Locus"), sep = ';') %>%
  dplyr::mutate("length"= (end - start) + 1)


#Formating to get the final table
Reformat <- parsedCDS %>%
  dplyr::group_by(Parent) %>%
  dplyr::summarise(CDNA = sum(length)) %>%
  dplyr::mutate(AA_Length = ((CDNA/3) - 1)) %>% #Calculating AA length
  dplyr::mutate(Parent = stringr::str_replace(Parent, ":", "=")) %>% # to allow for elegans and briggsae
  dplyr::rowwise() %>%
  dplyr::mutate(TRANSCRIPT = tail(stringr::str_split(Parent, "=")[[1]], n = 1)) %>%
  tidyr::separate_rows(TRANSCRIPT, sep = ",")

#Making final table
AA_Length <- Reformat %>%
  dplyr::select("TRANSCRIPT", "AA_Length")

readr::write_tsv(AA_Length, "gff_AA_Length.tsv")

















