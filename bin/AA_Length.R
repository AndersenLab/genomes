#!/usr/bin/env Rscript

#Libraries
library(dplyr)
library(stringr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

#Read in GFF
gff <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE)

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
  dplyr::mutate(AA_Length = ((CDNA/3) - 1)) #Calculating AA length


Reformat <- Reformat %>%
  dplyr::mutate("Transcript" = str_sub(Reformat$Parent, 19)) %>%
  tidyr::separate(Transcript, into = c("Transcript_1", "Rest"), sep = ",") #If not including transcript isoforms, this can be named TRANSCRIPT and just continue on to making the table.

#Optional steps to include isoforms with alternatively spliced 5' and 3' UTRS
Reformat <- Reformat %>%
  dplyr::mutate("Transcript_2" = str_sub(Reformat$Rest, 12)) %>%
  tidyr::unite("TRANSCRIPT", "Transcript_1", "Transcript_2", sep = ",", na.rm = TRUE)

Reformat <- Reformat %>%
  tidyr::separate_rows(TRANSCRIPT, sep = ",")


#Making final table
AA_Length <- Reformat %>%
  dplyr::select("TRANSCRIPT", "AA_Length")


readr::write_tsv(AA_Length, "gff_AA_Length.tsv")

















