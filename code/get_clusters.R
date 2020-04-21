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

function(srr, cutoff_df, complexity="high", cluster_cutoff, escore=1e-10, cores=12 ){
  # srr is an SRA accession number
  # cutoff_df is the output from dustyCutoffs
  # complexity can be "high", "low", or "all"
  # subtrees is the number of subclusters used to create consensus sequences
  # cores is the number of CPUs to use
  # consenus_threshold is the fraction of reads required to establish a consensus base.
  # output=taxa returns the taxonomy of each read
  # output=seq returns the consensus sequences
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
    return(clustersId)
}