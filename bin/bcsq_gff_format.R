#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)



#Read in GFF and Parse Columns
raw_gff <- read.delim("simple.wormbase.gff3", header=FALSE, stringsAsFactors=FALSE) 


format_gff <- function(df){ 
  #Case Change 
  case_change <- df %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "ID=Transcript:", "ID=transcript:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "ID=Gene:", "ID=gene:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "Parent=Gene:", "Parent=gene:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "Parent=Transcript:", "Parent=transcript:")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "Parent=Pseudogene:", "Parent=transcript:")) 
    
    #Gene Line Fix
  
  gtf_gene <- dplyr::filter(case_change, V3== "gene")
  
  gene_name <- tidyr::separate(gtf_gene, V9 , c("before", "after"), sep="Name=") %>% 
    tidyr::separate(after, c("Name","others"), sep=";") %>% 
    dplyr::select(pos1 = V4, Name)
  
  gene_biotype <- tidyr::separate(gtf_gene, V9, c("before","after"), sep="biotype=") %>% 
    tidyr::separate(after, c("biotype","others"), sep=";") %>% 
    dplyr::select(pos2 = V4, biotype) %>% 
    dplyr::mutate(biotype=paste0(";biotype=", biotype))
  
  name_biotype_table <- dplyr::bind_cols(gene_name, gene_biotype) %>%  
    #dplyr::filter(V4...1 == V4...3) %>%  
    dplyr::filter(pos1 == pos2) %>%
    dplyr::select(Name, biotype)
  
  tr <- dplyr::filter(case_change, stringr::str_detect(V9, "ID=transcript")) %>% 
    tidyr::separate(V9 ,c("before","after"), sep="Parent=gene:", remove=F) %>% 
    tidyr::separate(after, c("parent_gene","others"), sep=";") %>%  
    dplyr::left_join(name_biotype_table, by=c("parent_gene"="Name")) %>% 
    dplyr::select(V9, biotype)
  
  rebuilt <- dplyr::left_join(case_change, tr, by="V9") %>% 
    tidyr::unite(V9, c(V9, biotype), sep="", remove=T, na.rm=T)
  
    #CDS .1 Fix - Seperate and remove second transcript from ID=transcript string
  cds_1 <- rebuilt %>% 
    dplyr::filter(V3 == "CDS") %>%  
    tidyr::separate(V9, sep=";", into = c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8","a9")) %>% 
    tidyr::separate(a2, sep= "," , into= c("PARENT1", "PARENT2")) %>% 
    dplyr::select(-PARENT2) %>% 
    tidyr::unite(V9, a1:a9, sep= ";", na.rm= TRUE)
    
    
    #CDS .2 Fix - Build DF with second transcript and insert ID=transcript             
  cds_2 <- rebuilt %>% 
    dplyr::filter(V3 == "CDS")  %>% 
    tidyr::separate(V9, sep=";", into = c("a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8","a9")) %>% 
    tidyr::separate(a2, sep= "," , into= c("PARENT1", "PARENT2")) %>% 
    dplyr::mutate(PARENT2 = stringr::str_replace(PARENT2, "Transcript", "Parent=transcript")) %>% 
    dplyr::select(-PARENT1) %>% 
    tidyr::drop_na(PARENT2) %>% 
    tidyr::unite(V9, a1:a9, sep= ";", na.rm= TRUE)
    
    #Join .1 fix and .2 Fixed
  cds_1_2 <- dplyr::bind_rows(cds_1, cds_2)
    
  cds_fixed <- rebuilt %>% 
    dplyr::filter(V3 != "CDS") %>% 
    dplyr::bind_rows(cds_1_2)
    
   
  #Pseudogenic Transcript Fix - Add in biotype for transcript lines & Replace ID= field
  
  pseudogene_edit_1 <- cds_fixed %>% 
    dplyr::filter(V3 == "pseudogenic_transcript")  %>% 
    tidyr::separate(V9, sep=";", into = c("a1", "a2", "a3", "a4")) %>% 
    dplyr::mutate(a5 = "biotype=polymorphic_pseudogene") %>%  #Add in biotype
    tidyr::unite(V9, a1:a5, sep= ";", na.rm= TRUE) %>%
    dplyr::mutate(V9 = stringr::str_replace(V9, "ID=Pseudogene:", "ID=transcript:")) #Replace ID= field
  
  pseudogene_edit_2 <- cds_fixed %>% 
    dplyr::filter(V3 != "pseudogenic_transcript") 
  
  
  pseudogene_fixed <- dplyr::bind_rows(pseudogene_edit_1, pseudogene_edit_2 )
  
  
  
  
   #Convert to Supported Biotypes 
  
  corrected_biotypes <- pseudogene_fixed %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "biotype=asRNA", "biotype=antisense_RNA")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "biotype=pseudogene", "biotype=polymorphic_pseudogene")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "biotype=piRNA", "biotype=misc_RNA")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "biotype=ncRNA", "biotype=misc_RNA")) %>% 
    dplyr::mutate(V9 = stringr::str_replace(V9, "biotype=tRNA", "biotype=misc_RNA"))
    
    
  
      #Final GFF
  
  reformatted_gff <- corrected_biotypes[with(corrected_biotypes, order(V1, V4)), ]
  
  return(reformatted_gff)
}

BCSQ_gff <- format_gff(raw_gff)

#Write out the GFF file

write.table(BCSQ_gff, "fixed.wormbase.gff3", quote=F, sep="\t", col.names = F, row.names = F)



