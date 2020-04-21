#!/usr/bin/R
library(dplyr)
library(stringr)

function(consensus_tax, cutoff_df, dend_list, group=3) {
  regexp <- "[[:digit:]]+"
  idx_dend = as.numeric(unique(str_extract(consensus_tax$QueryID, regexp)))
  
  blast_labels <- data.frame()
  for (i in seq(1:length(dend_list))) {
    blast_labels <- c(blast_labels, as.numeric(dend_list[[idx_dend[i]]] %>% labels))
  }
  
  high_complexity_count <- cutoff_df %>% count(cluster==group)
  recovered = length(blast_labels)/high_complexity_count$n[[2]]
  return(recovered)
}