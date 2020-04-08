#!/usr/bin/R
library(rentrez)
library(Biostrings)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)


function(srr,blasthits,cutoff=5000,amplicon_length=550){
  #srr is an SRA accession number
  #blasthits is a BLAST result data frame
  #cutoff is a number from 1 to maximum dusty score of srr. Default cutoff=5000 
  system(paste('fastq-dump', srr, sep=' ')) 
  dna = readFastq('.', pattern=srr)

  reads = sread(dna)
  qscores = quality(dna) 

  print('Splitting into reads with and without BLAST hits')

  HitsQID_group = blasthits %>%
    group_by(QueryID) 


  #extract read indices
  regexp <- "[[:digit:]]+"
  idxHits = as.numeric(unique(str_extract(HitsQID_group$QueryID, regexp)))
  idxNA = unique(setdiff((seq(1:length(dna))),idxHits))
  NA_reads = dna[idxNA]

############ Quality scores of Hits
  print('Calculating quality scores for each read')
##NA
  QscoresNA = quality(dna[idxNA]) # parse quality scores
  numQscoresNA = as(QscoresNA, "matrix") # converts to numeric scores automatically
  meanQscoresNA_temp = rowMeans(numQscoresNA,na.rm=TRUE)
  meanQscoresNA = as.data.frame(meanQscoresNA_temp)
  remove(meanQscoresNA_temp)
  
  ##Hits
  QscoresHits = quality(dna[idxHits])
  numQscoresHits = as(QscoresHits, "matrix")
  #meanQscoresHits_temp = rowMeans(numQscoresHits,na.rm=TRUE)
  meanQscoresHits_temp = apply(numQscoresHits, 1, mean, na.rm=TRUE) #apply the function mean() across 
  meanQscoresHits = as.data.frame(meanQscoresHits_temp,col.names=c('QmeansHits'))
  remove(meanQscoresHits_temp)
  
  print('Calculating dusty scores for each read')
  complexNA = dustyScore(dna[idxNA])
  complexNA = as.data.frame(complexNA)
  meanQscores_dustyNA = cbind(meanQ=(meanQscoresNA$meanQscoresNA_temp),Dusty=(complexNA$complexNA))
  meanQscores_dustyNA = as.data.frame(meanQscores_dustyNA)
  
  complexHits = dustyScore(dna[idxHits])
  complexHits = as.data.frame(complexHits)
  meanQscores_dustyHits = cbind(meanQ=(meanQscoresHits$meanQscoresHits_temp),Dusty=(complexHits$complexHits))
  meanQscores_dustyHits = as.data.frame(meanQscores_dustyHits)


  print(paste('Grouping into high and low complexity populations with cutoff = ',cutoff,sep=""))
  

  clusterHits <- as.numeric((meanQscores_dustyHits$Dusty >= cutoff)+1)
  meanQscores_dustyHits$cluster <- as.factor(clusterHits)
  #clusterCentersHits <- as.data.frame(clusterHits$centers,clusterHits$cluster)
  clusterCentersHits <- as.data.frame(rbind(c(Q=mean(meanQscores_dustyHits$meanQ[meanQscores_dustyHits$cluster==1]), D=mean(meanQscores_dustyHits$Dusty[meanQscores_dustyHits$cluster==1])),c(Q=mean(meanQscores_dustyHits$meanQ[meanQscores_dustyHits$cluster==2]), D=mean(meanQscores_dustyHits$Dusty[meanQscores_dustyHits$cluster==2]))))
  
  #clusterNA = kmeans(meanQscores_dustyNA,2)
  clusterNA <- as.numeric((meanQscores_dustyNA$Dusty < cutoff)+3)
  meanQscores_dustyNA$cluster <- as.factor((clusterNA))
  #clusterCentersNA <- as.data.frame(clusterNA$centers,clusterNA$cluster)
  clusterCentersNA <- as.data.frame(rbind(c(Q=mean(meanQscores_dustyNA$meanQ[meanQscores_dustyNA$cluster==3]), D=mean(meanQscores_dustyNA$Dusty[meanQscores_dustyNA$cluster==3])),c(Q=mean(meanQscores_dustyNA$meanQ[meanQscores_dustyNA$cluster==4]), D=mean(meanQscores_dustyNA$Dusty[meanQscores_dustyNA$cluster==4]))))
  
  
  #Overlayed plot
  print('Plotting all clusters')
  
  histcomb_hits_label  <- paste("BLAST Hits (",length(idxHits),")", sep="")
  histcomb_NA_label  <- paste("No BLAST Hits (",length(idxNA),")", sep="")
  histcomb = (ggplot(meanQscores_dustyNA) + 
                geom_point(data=meanQscores_dustyHits,
                          aes(x=meanQ,y=Dusty,color=as.factor(cluster)),size=0.02) +
                geom_point(data=meanQscores_dustyNA,
                           aes(x=meanQ,y=Dusty,color=as.factor(cluster)),size=0.02) +
                #
                stat_density_2d(data=meanQscores_dustyHits,
                                aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
                stat_density_2d(data=meanQscores_dustyNA,
                                aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
                #
                geom_point(data=clusterCentersHits,
                           aes(x=clusterCentersHits$Q,y=clusterCentersHits$D,shape=as.factor(histcomb_hits_label),size=3)) +
                geom_point(data=clusterCentersNA, 
                           aes(x=clusterCentersNA$Q,y=clusterCentersNA$D,shape=as.factor(histcomb_NA_label),size=3)) +
                
                scale_y_log10() +
                xlab('Mean Q-Score (Nanopore)') +
                ylab(expression("Dusty Complexity Score (log"[10]*")")) +
                labs(color="Cluster",shape="Dataset") +
                ggtitle(paste('Grouping of BLAST Results of a ', amplicon_length,' (bp) \n Amplicon by Complexity', sep='')) +
                guides(size=FALSE) +
                theme_minimal())
  
  return(histcomb)
}

#srr = c('SRR11206993') 

#dusty_Q_plot

