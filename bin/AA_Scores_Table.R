#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

#Grantham Formatting
grantham <- readr::read_tsv("https://gist.githubusercontent.com/danielecook/501f03650bca6a3db31ff3af2d413d2a/raw/5583a134b36b60762be6cd54002a0f4044338cd0/grantham.tsv")

grantham_score <- grantham %>%
  tidyr::gather(SECOND,GSCORE, -FIRST) %>%
  dplyr::filter(GSCORE > 0)

grantham_score_a <- grantham_score %>%
  dplyr::rename(AA = FIRST, ALT_AA = SECOND)

#Need AA pairs to be in both options (AB and BA)
grantham_score_b <- grantham_score %>%
  dplyr::rename(AA = SECOND, ALT_AA = FIRST)

g_score <- dplyr::bind_rows(grantham_score_a, grantham_score_b)


#BLOSUM Formatting
BLOSUM62 <- read.table(args[1], header=TRUE, quote="\"", stringsAsFactors=FALSE, row.names = NULL)

blosum_score <- BLOSUM62 %>%
  tidyr::gather(AA,BSCORE, -row.names) %>%
  dplyr::rename(ALT_AA = row.names)

#Unite Scores into 1 Table
scores <- dplyr::left_join(blosum_score, g_score)

scores <- scores %>%
  dplyr::select("AA","ALT_AA", "GSCORE", "BSCORE")

#Combine REF and ALT
united_scored <- scores %>%
  tidyr::unite("REF_ALT_AA", "AA", "ALT_AA", sep= "|")

#Write out table
readr::write_tsv(united_scored, "AA_Scores.tsv")






















