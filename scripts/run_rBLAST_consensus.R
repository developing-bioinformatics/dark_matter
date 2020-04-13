#!/usr/bin/R
library(rBLAST)
library(parallel)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(Biostrings)
library(DECIPHER)

setwd('/projectnb/ct-shbioinf/epederso/dark_matter/')
dustyCutoffs <- dget('scripts/dustyCutoffs.R')
rBLAST_consensus <- dget('scripts/rBLAST_consensus.R')
              #1            #2            #3          #4              #5
srr_list = c("SRR11043467","SRR11043468","SRR11043469","SRR11043470","SRR11043471", 
             "SRR11043472","SRR11043473","SRR11043474","SRR11043475","SRR11043476", #10
             "SRR11043477","SRR11043478","SRR11043479","SRR11043480","SRR11043481",
             "SRR11043482","SRR11043483","SRR11043484","SRR11043486","SRR11043487", #20
             "SRR11043488","SRR11043489","SRR11043490","SRR11043491","SRR11043493",
             "SRR11043494","SRR11043495","SRR11043496","SRR11043497","SRR11043498", #30
             "SRR11043499","SRR11043500","SRR11043501","SRR11043502","SRR11043503",
             "SRR11043504","SRR11043505","SRR11043506","SRR11043507","SRR11043508", #40
             "SRR11043509","SRR11043510","SRR11043511","SRR11043512","SRR11043513",
             "SRR11043514","SRR11043515","SRR11043516","SRR11043517","SRR11043492", #50
             "SRR11206993","SRR11206994","SRR11206995","SRR11206996","SRR11206997",
             "SRR11206998","SRR11206999","SRR11043485")                             #58



cores = c(16)
escore = c('1e-10')
idx_srr = 51

b_args <- paste('-num_threads ',cores,' -evalue ',escore, sep="")
srr = srr_list[idx_srr]

dna = readFastq('./data/', pattern=srr_list[idx_srr])
reads = sread(dna)

#BLAST search
bl_ncbi <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
blasthits_ncbi <- predict(bl_ncbi, reads, BLAST_args = b_args)

cutoff_plot <- dustyCutoffs(srr_list[idx_srr],
                            blasthits_ncbi,
                            cutoff=5000,
                            data=FALSE,
                            verbose=TRUE)
cutoff_plot #plot the object
cutoff_df <- dustyCutoffs(srr_list[idx_srr], #collect data frame
                          blasthits_ncbi,
                          cutoff=5000,
                          data=TRUE,
                          verbose=TRUE)
#Consensus sequence data
consensus_tax_high <- rBLAST_consensus(srr_list[idx_srr],
                                       cutoff_df,
                                       complexity="high",
                                       subtrees=30,
                                       cluster_cutoff=0.01,
                                       escore=escore,
                                       output='taxa',
                                       cores=cores)
consensus_tax_all <- rBLAST_consensus(srr_list[idx_srr],
                                      cutoff_df,
                                      complexity="all",
                                      subtrees=70,
                                      cluster_cutoff=0.01,
                                      escore=escore,
                                      output='taxa',
                                      cores=cores)

cluster_consensus <- rBLAST_consensus(srr_list[idx_srr],
                                      cutoff_df,
                                      complexity="high",
                                      subtrees=70,
                                      cluster_cutoff=0.01,
                                      escore=escore,
                                      output='seqs', 
                                      cores=cores)
dend_list <- rBLAST_consensus(srr_list[idx_srr],
                              cutoff_df,
                              complexity="high",
                              subtrees=40,
                              cluster_cutoff=0.01,
                              escore=escore,
                              output='dend',
                              cores=cores)

# save data image
save.image(paste('outputs/save',srr_list[idx_srr],'_consensus.RData',sep=''))

# write FASTA of consensus sequences
writeFasta(cluster_consensus,paste('outputs/consensus',srr_list[idx_srr],'.fa',sep=''))

