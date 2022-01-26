#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)

# arguments:
# 1: gff file
# 2: species  
args <- commandArgs(trailingOnly = TRUE)

# load gff
gff <- read.delim(args[1], header=FALSE, stringsAsFactors=FALSE)

# change column names
colnames(gff) <- c('CHROM', "source", "type", "start", "stop", "v6", "v7", "v8", "info")

# re-arrange file for NemaScan
gene_ref <- gff %>%
  dplyr::filter(type %in% c("mRNA", "transcript")) %>%
  dplyr::mutate(gene = stringr::str_match(info, "transcript:\\s*(.*?)\\s*;")[,2],
                wbgene = stringr::str_match(info, "gene:\\s*(.*?)\\s*;")[,2],
                biotype = stringr::str_split_fixed(info, "biotype=", 2)[,2],
                biotype = stringr::str_split_fixed(biotype, ";", 2)[,1]) %>%
  dplyr::mutate(gene_name = ifelse(grepl("locus", info), stringr::str_match(info, "locus=\\s*(.*?)\\s*;")[,2], gene)) %>%
  dplyr::select(gene, chr = CHROM, strand = v7, txstart = start, txend = stop, wbgene, gene_name, biotype) %>%
  dplyr::mutate(type = "Transcript")

# save
readr::write_tsv(gene_ref, paste0(args[2], ".gff"))