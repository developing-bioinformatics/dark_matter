#!/usr/bin/R
library(rBLAST)
library(parallel)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)
library(Biostrings)
library(taxonomizr)
library(DECIPHER)
library(dendextend)
library(ggdendro)
library(circlize)


setwd('/projectnb/ct-shbioinf/epederso/dark_matter/')
rBLAST_pre <- dget('code/rBLAST_pre.R') # run BLAST
dustyCutoffs <- dget('code/dustyCutoffs.R') # identify high-complexity reads that do not return BLAST results
get_clusters <- dget('code/get_clusters.R') # calculate distance matrix and return a dendrogram
# In script as tryCatch: calculate the optimal number of subtrees so that each does not contain singletons
cluster_split <- dget('code/cluster_split.R') # split clusters based on the optimal number of subtrees
get_consensus <- dget('code/get_consensus.R') # align sequences based on each subtrees to get consensus sequences
consensus_taxa <- dget('code/consensus_taxa.R') # BLAST consensus sequences and identify taxa
fracIdentity <- dget('code/fracIdentity.R') # calculate the fraction of high-complexity reads that are now represented by BLAST hits
lca_core <- dget('code/lca_core.R') # for coloring dendrograms

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
p = proc.time()
idx_srr = 9
srr = srr_list[idx_srr]
complexity_group <- "high"


blasthits_ncbi <- rBLAST_pre(srr_list,idx_srr,escore,cores)


cutoff_df <- dustyCutoffs(srr_list[idx_srr], #collect data frame
                          blasthits_ncbi,
                          cutoff=3000,
                          data=TRUE,
                          verbose=TRUE)

cutoff_plot <- dustyCutoffs(srr_list[idx_srr],
                            blasthits_ncbi,
                            cutoff=3000,
                            data=FALSE,
                            verbose=TRUE)
cutoff_plot #plot the object
try({
  clusterMax <- get_clusters(srr_list[idx_srr],
                             cutoff_df,
                             cluster_cutoff=0.01,
                             complexity=complexity_group, 
                             escore, 
                             cores)
  },silent=TRUE)

#Maximize number of subtrees with tryCatch loop. Returns highest possible number of subclusters.
num_reads <-length(cutoff_df$cluster)
#subtrees <- as.integer(num_reads/450)
subtrees <- 400
cycles <- (subtrees+seq(1:num_reads))
i <- subtrees
tryCatch(
  expr = {
    for (i in cycles){
      
    dend_list <- cluster_split(clusterMax,i)
    consensus_seqs <- get_consensus(srr=srr_list[idx_srr],
                                    dend_list,
                                    subtrees=i, 
                                    cluster_cutoff=0.01, 
                                    stag_threshold=5,
                                    consensus_threshold=0.5,
                                    consensus_min=0.5,
                                    complexity=complexity_group,
                                    cores=cores)
    }
  },
  error = function(i){
    NA
  },
  finally = {
    optSubtrees <- i-1
    dend_list <- cluster_split(clusterMax,optSubtrees)
    consensus_seqs <- get_consensus(srr=srr_list[idx_srr],
                                      dend_list,
                                      subtrees=optSubtrees, 
                                      cluster_cutoff=0.01, 
                                      stag_threshold=5,
                                      consensus_threshold=0.5,
                                      consensus_min=0.5,
                                      complexity=complexity_group,
                                      cores=cores)
  }
)

BrowseSeqs(consensus_seqs) #visual inspection
 
consensus_tax_high <- consensus_taxa(consensus_seqs, 
                                     escore, 
                                     cores)
 
#unique(consensus_tax_high$QueryID)

recovered <- fracIdentity(consensus_tax_high, 
                          cutoff_df, 
                          dend_list, 
                          complexity=complexity_group)

# save cluster plot
out1 <- paste("outputs/dusty_clusters_",srr,".png", sep="")
ggsave(cutoff_plot, file=out1, height = 9, width = 7, dpi=500)

#proc.time() - p


# write FASTA of consensus sequences
writeFasta(consensus_seqs,paste('outputs/consensus',srr_list[idx_srr],'_',complexity_group,'.fa',sep=''))


mean(cutoff_df$Dusty[cutoff_df$cluster==1])
sd(cutoff_df$Dusty[cutoff_df$cluster==1])
length(cutoff_df$Dusty[cutoff_df$cluster==1])

mean(cutoff_df$Dusty[cutoff_df$cluster==3])
sd(cutoff_df$Dusty[cutoff_df$cluster==3])
length(cutoff_df$Dusty[cutoff_df$cluster==3])

mean(cutoff_df$Dusty[cutoff_df$cluster==2])
sd(cutoff_df$Dusty[cutoff_df$cluster==2])
length(cutoff_df$Dusty[cutoff_df$cluster==2])

mean(cutoff_df$Dusty[cutoff_df$cluster==4])
sd(cutoff_df$Dusty[cutoff_df$cluster==4])
length(cutoff_df$Dusty[cutoff_df$cluster==4])

mean(cutoff_df$meanQ[cutoff_df$cluster==1])
sd(cutoff_df$meanQ[cutoff_df$cluster==1])
length(cutoff_df$meanQ[cutoff_df$cluster==1])

mean(cutoff_df$meanQ[cutoff_df$cluster==3])
sd(cutoff_df$meanQ[cutoff_df$cluster==3])
length(cutoff_df$meanQ[cutoff_df$cluster==3])

mean(cutoff_df$meanQ[cutoff_df$cluster==2])
sd(cutoff_df$meanQ[cutoff_df$cluster==2])
length(cutoff_df$meanQ[cutoff_df$cluster==2])

mean(cutoff_df$meanQ[cutoff_df$cluster==4])
sd(cutoff_df$meanQ[cutoff_df$cluster==4])
length(cutoff_df$meanQ[cutoff_df$cluster==4])

##### dendrograms

lca_df = consensus_tax_high %>%
  group_modify(~ lca_core(.x)) %>% # run lca
  slice(1) %>% #keep only first row in each group (one per read)
  summarize(last_common) 

lca_df[is.na(lca_df)] <- "Unknown"

regexp <- "[[:digit:]]+"
taxa_idx = as.numeric(str_extract(lca_df$QueryID, regexp))
taxa_labels = lca_df$last_common
taxa <- data.frame()
taxa <- cbind(taxa_idx,taxa_labels)
names <- unique(lca_df$last_common)

idx_dend = as.numeric(unique(str_extract(consensus_tax_high$QueryID, regexp))) %>% sort()
color_vec <- matrix(1,optSubtrees,1)
color_vec <- as.numeric(color_vec)
n= length(taxa_idx)
for (i in seq_len(n)){
  for (j in seq(1:length(names))){
    if (taxa[n+i] == names[j]){
      color_vec[as.numeric(taxa[i])] = j+1
    }
  }
}

circos.par(cell.padding = c(0, 0, 0, 0))
labels <- as.numeric(clusterMax %>% 
                       labels)
m = length(labels)
leg_names = c(names,'No BLAST Hits')
if (length(names) > 1){
  leg_col = c(1+seq(2:(length(names)+1)),1)
} else {
  leg_col = c(2,1)
}

if (complexity_group == "high"){
  main_title <- paste("High Complexity Dark Matter ",srr_list[idx_srr],sep="")
} else if (complexity_group == "low"){
  main_title <- paste("Low Complexity Dark Matter ",srr_list[idx_srr],sep="")
} else if (complexity_group == "all"){
  main_title <- paste("Dark Matter ",srr_list[idx_srr],sep="")
}

circos.initialize(factors = "a", xlim = c(0, m)) # only one sector
dend = color_branches(clusterMax, k = optSubtrees, col = color_vec)
dend_height = attr(dend, "height")
circos.track(ylim = c(0, dend_height), bg.border = NA, 
             track.height = 0.4, panel.fun = function(x, y) {
               circos.dendrogram(dend)
             })
title(main=main_title)
legend(x=-0.3, y=0.5, legend=leg_names, bty='n', text.col=leg_col, cex=0.6)
circos.clear()

# save data image
save.image(paste('outputs/save',srr_list[idx_srr],'_consensus_',complexity_group,'.RData',sep=''))
