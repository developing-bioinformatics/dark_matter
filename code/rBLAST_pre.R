#!/usr/bin/R
library(rBLAST)
library(parallel)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)
library(Biostrings)

function(srr_list, idx_srr, escore, cores=12 ){
  srr <- srr_list[idx_srr]
  if (file.exists(paste('data/',srr,'.fastq',sep=''))==FALSE) {
    #system(paste('fastq-dump', srr, sep=' '))
    SRAfile <- read_xml(entrez_fetch(db="sra",id=srr,rettype='xml')) %>%
      xml_find_all("//RUN_SET") %>% 
      xml_find_all("//SRAFile") %>% 
      xml_attr('url')
    download.file(SRAfile[1], paste(srr,'.fastq',sep=""))
  }
#srr = srr_list[idx_srr]

dna = readFastq('./data/', pattern=srr_list[idx_srr])
reads = sread(dna)


bl_ncbi <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
nclus = 2
b_args <- paste('-num_threads ',cores/nclus,' -evalue ',escore, sep="")

blastit = function(reads) {
  return(predict(bl_ncbi, reads, BLAST_args= b_args))
}

p = proc.time()
cl = makeCluster(nclus, type = 'FORK')
splits = clusterSplit(cl, reads)
p_works = parLapply(cl, splits, blastit)
stopCluster(cl)
proc.time() - p
blasthits_ncbi_pre <- rbind(p_works[[1]],p_works[[2]])
blasthits_ncbi = blasthits_ncbi_pre %>% 
  arrange()

return(blasthits_ncbi)
}