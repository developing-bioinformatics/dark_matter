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

function(consensus_seqs, escore=1e-10, cores=12){
  
  taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
  taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")
  
  bl_ncbi <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
  nclus = 2
  b_args <- paste('-num_threads ',cores/nclus,' -evalue ',escore, sep="")
  
  blastit = function(consensus_seqs) {
    return(predict(bl_ncbi, consensus_seqs, BLAST_args= b_args))
  }
  
  #p = proc.time()
  cl = makeCluster(nclus, type = 'FORK')
  splits = clusterSplit(cl, consensus_seqs)
  p_works = parLapply(cl, splits, blastit)
  stopCluster(cl)
  #proc.time() - p
  blasthits_pre <- rbind(p_works[[1]],p_works[[2]])
  blasthits <- blasthits_pre %>% arrange()
  accid = as.character(blasthits$SubjectID)
  ids<-accessionToTaxa(accid, '/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
  taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
  consensus_tax=cbind(blasthits,taxlist)
  return(consensus_tax)
}