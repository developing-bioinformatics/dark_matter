#!/usr/bin/R
library(Biostrings)
library(taxonomizr)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)
library(stringi)
library(forcats)
library(rBLAST)
library(multidplyr)
library(cowplot)
library(ape)
library(DECIPHER)
#library(kmer)
library(msa)
#install.packages("kmer")

srr=c('SRR11206999') 

cluster_4_idx = (meanQscores_dustyNA$cluster==4)
NA_reads_cluster_4 = dnaMerge[cluster_4_idx]
cluster_4_sreads = NA_reads_cluster_4@sread

distanceMat <- DistanceMatrix(cluster_4_sreads, processors=NULL) #compute distance matrix

clustersId <- IdClusters(distanceMat, method="complete", cutoff=0.01,showPlot = TRUE, myXStringSet = cluster_4_sreads, type = "dendrogram", processors=NULL)
#clustersId <- IdClusters(method="inexact", cutoff=0.01, myXStringSet = cluster_4_sreads, processors=NULL)


cluster_4_alignment <- AlignSeqs(cluster_4_sreads, iterations=4, gapOpening=c(-30,-10),gapExtension=c(-5,-1), perfectMatch = 15,misMatch = -3, processors=NULL)
BrowseSeqs(cluster_4_alignment)

writeXStringSet(cluster_4_alignment, file=paste("./outputs/OHara/Alignments_cluster_4_",srr,".fasta",sep=""))

#cluster_4_stagger_alignment <- StaggerAlignment(cluster_4_sreads, tree = NULL, threshold = 3, fullLength = FALSE, processors = 1, verbose = TRUE)

IUPAC_CODE_MAP # list of ambiguity codes for reference
threshold = 0.4
cluster_4_consensus <- ConsensusSequence(cluster_4_alignment,
                                     threshold = 0.4,
                                     ambiguity = TRUE,
                                     noConsensusChar = "+",
                                     minInformation = 1 - threshold,
                                     includeTerminalGaps = TRUE)
BrowseSeqs(cluster_4_alignment)
BrowseSeqs(cluster_4_consensus)

#####

#correlate Dusty and Lengths

widths_dustyNA = cbind.data.frame(widths=widths$`reads@ranges@width`[idxNA],dusty=complexNA$complexNA)
ggplot(widths_dustyNA) +
  geom_point(aes(x=widths,y=dusty)) +
  scale_y_log10() +
  #geom_point(data=clusterCentersHits, aes(x=clusterCentersHits$meanQ,y=clusterCentersHits$Dusty)) +
  theme_minimal()

widths_dusty = cbind.data.frame(widths=widths$`reads@ranges@width`,dusty=dustyScore(dnaMerge))

# create multiple linear model
lm_fit <- lm(widths ~ dusty, data=widths_dusty)
summary(lm_fit)

# save predictions of the model in the new data frame 
# together with variable you want to plot against
corr_df <- data.frame(widths_pred = predict(lm_fit, widths_dusty),dusty=widths_dusty$dusty)

ggplot(widths_dusty) +
  geom_point(aes(x=widths,y=dusty)) +
  #scale_y_log10() +
  geom_line(data=corr_df,aes(x=widths_pred,y=dusty)) +
  theme_minimal()
