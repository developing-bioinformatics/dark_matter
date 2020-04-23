#!/usr/bin/R
library(rBLAST)
library(Biostrings)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)
library(DECIPHER)
library(taxonomizr)
library(dendextend)
library(parallel)
options(scipen=999, digits=1) 

function(srr, dend_list, subtrees=50, cluster_cutoff=0.01, stag_threshold=5,consensus_threshold=0.5,consensus_min=0.5,complexity="high",cores=12 ){

  
  if (file.exists(paste('data/',srr,'.fastq',sep=''))==FALSE) {
    #system(paste('fastq-dump', srr, sep=' '))
    SRAfile <- read_xml(entrez_fetch(db="sra",id=srr,rettype='xml')) %>%
      xml_find_all("//RUN_SET") %>% 
      xml_find_all("//SRAFile") %>% 
      xml_attr('url')
    download.file(SRAfile[1], paste(srr,'.fastq',sep=""))
  }
  
  dna = readFastq(paste('data/',srr,'.fastq',sep=''))
  if (complexity=="high"){
    cluster_sreads = sread(dna[cutoff_df$cluster==4])
  } else if (complexity=="low") {
    cluster_sreads = sread(dna[cutoff_df$cluster==3])
  } else if (complexity=="all") {
    cluster_sreads = sread(dna[cutoff_df$cluster==3|cutoff_df$cluster==4])
  }
  
  #initialize data frames for for loop
  dend_seqs_df <- data.frame()
  dend_alignment_df <- data.frame()
  dend_stag_df <- data.frame()
  cluster_consensus_gap <- data.frame()
  cluster_consensus <- DNAStringSet()
  
  
  for (i in seq(1:subtrees)) {
    cat(paste('Subcluster ',i,'\n',sep=''))
    idx_test <- as.numeric(dend_list[[i]] %>% 
                             labels)
    dend_seqs_df <- cluster_sreads[idx_test]
    idx_test <- data.frame()
    dend_alignment_df <-AlignSeqs(dend_seqs_df, processors=cores)
    dend_stag_df <- StaggerAlignment(dend_alignment_df,threshold=stag_threshold, processors=cores)
    cluster_consensus_gap <- c(cluster_consensus_gap,ConsensusSequence(dend_stag_df,
                                           threshold = consensus_threshold,
                                           ambiguity=FALSE,
                                           noConsensusChar = "-",
                                           minInformation = consensus_min,
                                           includeTerminalGaps = FALSE))
    cluster_consensus = c(cluster_consensus,RemoveGaps(cluster_consensus_gap[[i]], removeGaps = "all", processors=cores))
  }
    return(cluster_consensus)
}