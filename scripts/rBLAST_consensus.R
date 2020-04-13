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
options(scipen=999, digits=1) 

function(srr, cutoff_df, complexity="high", subtrees=50, cluster_cutoff=0.01, escore=1e-10, consensus_threshold=0.8, output='taxa', cores=12 ){
  # srr is an SRA accession number
  # cutoff_df is the output from dustyCutoffs
  # complexity can be "high", "low", or "all"
  # subtrees is the number of subclusters used to create consensus sequences
  # cores is the number of CPUs to use
  # consenus_threshold is the fraction of reads required to establish a consensus base.
  # output=taxa returns the taxonomy of each read
  # output=seqs returns the consensus sequences
  # output=dend returns the list of subdendrograms
  
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
    cluster_sreads = sread(dna[cutoff_df$cluster==3])
  } else if (complexity=="low") {
    cluster_sreads = sread(dna[cutoff_df$cluster==4])
  } else if (complexity=="all") {
    cluster_sreads = sread(dna[cutoff_df$cluster==3|cutoff_df$cluster==4])
  }
  
  distanceMat <- DistanceMatrix(cluster_sreads, processors=cores) #compute distance matrix
  clustersId <- IdClusters(distanceMat, 
                           method="complete",
                           cutoff=cluster_cutoff,
                           showPlot = FALSE,
                           myXStringSet = cluster_sreads, 
                           type = "dendrogram", 
                           processors=cores)
  
  dend_list <- get_subdendrograms(clustersId, subtrees)
  #initialize data frames for for loop
  dend_seqs_df <- data.frame()
  dend_alignment_df <- data.frame()
  dend_stag_df <- data.frame()
  cluster_consensus_gap <- data.frame()
  cluster_consensus <- DNAStringSet()
  
  
  for (i in seq(1:subtrees)) {
    cat(paste('Subcluster ',i,sep=''))
    idx_test <- as.numeric(dend_list[[i]] %>% 
                             labels)
    dend_seqs_df <- cluster_sreads[idx_test]
    idx_test <- data.frame()
    dend_alignment_df <-AlignSeqs(dend_seqs_df, processors=cores)
    dend_stag_df <- StaggerAlignment(dend_alignment_df, processors=cores)
    cluster_consensus_gap <- c(cluster_consensus_gap,ConsensusSequence(dend_stag_df,
                                           threshold = consensus_threshold,
                                           ambiguity = FALSE,
                                           noConsensusChar = "+",
                                           minInformation = 1 - consensus_threshold,
                                           includeTerminalGaps = TRUE))
    cluster_consensus = c(cluster_consensus,RemoveGaps(cluster_consensus_gap[[i]], removeGaps = "all", processors=cores))
  }


  taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
  taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")
  
  bl_ncbi <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
  b_args <- paste('-num_threads ',cores,' -evalue ',escore, sep="")
  blasthits <- predict(bl_ncbi, cluster_consensus, BLAST_args = b_args)
  
  accid = as.character(blasthits$SubjectID)
  ids<-accessionToTaxa(accid, '/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
  taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
  consensus_tax=cbind(blasthits,taxlist)
  
  if (output=='taxa') {
    return(consensus_tax)
  } else if (output=='seq') {
    return(cluster_consensus)
  } else if (output=='dend') {
    return(dend_list)
  }
}