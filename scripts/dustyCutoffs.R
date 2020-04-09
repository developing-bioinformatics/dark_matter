#!/usr/bin/R
library(rentrez)
library(xml2)
library(Biostrings)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringr)
options(scipen=999, digits=1) 

function(srr, blasthits, cutoff=5000, amplicon_length=550, data=FALSE, verbose=TRUE){
  #srr is an SRA accession number (string)
  #blasthits is a BLAST result data frame (before or after taxonomizr)
  #cutoff is a number from 1 to maximum dusty score of srr. 
  #Default cutoff = 5000 
  #If cutoff is NULL, kmeans clustering is used
  #If data = TRUE, the function returns a data frame with quality scores, complexity scores and cluster ids with length of dna .
  #If data = FALSE, the function returns a ggplot object.
  
  if (file.exists(paste(srr,'.fastq',sep=''))==FALSE) {
    #system(paste('fastq-dump', srr, sep=' '))
    SRAfile <- read_xml(entrez_fetch(db="sra",id=srr,rettype='xml')) %>%
      xml_find_all("//RUN_SET") %>% 
      xml_find_all("//SRAFile") %>% 
      xml_attr('url')
    download.file(SRAfile[1], paste(srr,'.fastq',sep=""))
  }
  
  dna = readFastq('.', pattern=srr)
  if (verbose==TRUE){
    cat('Extracting metadata\n')
  }
  
  reads = sread(dna)
  qscores = quality(dna) 
  escore = max(cl$E)
  
  if (verbose==TRUE){
    cat('Splitting into reads with and without BLAST hits\n')
  }
  
  #extract read indices
  regexp <- "[[:digit:]]+"
  idxHits = as.numeric(unique(str_extract(blasthits$QueryID, regexp)))
  idxNA = unique(setdiff((seq(1:length(dna))),idxHits))
  NA_reads = dna[idxNA]

## Quality scores of Hits
  if (verbose==TRUE){
    cat('Calculating quality scores for each read\n')
  }
##NA
  QscoresNA = quality(dna[idxNA]) # parse quality scores
  numQscoresNA = as(QscoresNA, "matrix") # converts to numeric scores automatically
  meanQscoresNA_temp = rowMeans(numQscoresNA,na.rm=TRUE)
  meanQscoresNA = as.data.frame(meanQscoresNA_temp)
  remove(meanQscoresNA_temp)
  
  ##Hits
  QscoresHits = quality(dna[idxHits])
  numQscoresHits = as(QscoresHits, "matrix")
  meanQscoresHits_temp = apply(numQscoresHits, 1, mean, na.rm=TRUE) #apply the function mean() across 
  meanQscoresHits = as.data.frame(meanQscoresHits_temp,col.names=c('QmeansHits'))
  remove(meanQscoresHits_temp)
  
  if (verbose==TRUE){
    cat('Calculating dusty scores for each read\n')
  }
  complexNA = dustyScore(dna[idxNA])
  complexNA = as.data.frame(complexNA)
  meanQscores_dustyNA = cbind(meanQ=(meanQscoresNA$meanQscoresNA_temp),Dusty=(complexNA$complexNA))
  meanQscores_dustyNA = as.data.frame(meanQscores_dustyNA)
  
  complexHits = dustyScore(dna[idxHits])
  complexHits = as.data.frame(complexHits)
  meanQscores_dustyHits = cbind(meanQ=(meanQscoresHits$meanQscoresHits_temp),Dusty=(complexHits$complexHits))
  meanQscores_dustyHits = as.data.frame(meanQscores_dustyHits)

  if (is.null(cutoff)) { #cluster if cutoff=NULL
    if (verbose==TRUE){
      cat('K-means clustering into high and low complexity populations\n')
    }
    
    clusterNA = kmeans(meanQscores_dustyNA,2)
    meanQscores_dustyNA$cluster <- as.factor((clusterNA$cluster+2))
    clusterCentersNA <- as.data.frame(clusterNA$centers,clusterNA$cluster)
  
    clusterHits = kmeans(meanQscores_dustyHits,2)
    meanQscores_dustyHits$cluster <- as.factor(clusterHits$cluster)
    clusterCentersHits <- as.data.frame(clusterHits$centers,clusterHits$cluster)

  } else {
    if (verbose==TRUE){
      cat(paste('Grouping into high and low complexity populations with cutoff = ',cutoff,"\n",sep=""))
    }
    clusterHits <- as.numeric((meanQscores_dustyHits$Dusty < cutoff)+1)
    meanQscores_dustyHits$cluster <- as.factor(clusterHits)
    clusterCentersHits <- as.data.frame(rbind(c(Q=mean(meanQscores_dustyHits$meanQ[meanQscores_dustyHits$cluster==1]), D=mean(meanQscores_dustyHits$Dusty[meanQscores_dustyHits$cluster==1])),c(Q=mean(meanQscores_dustyHits$meanQ[meanQscores_dustyHits$cluster==2]), D=mean(meanQscores_dustyHits$Dusty[meanQscores_dustyHits$cluster==2]))))
  
    clusterNA <- as.numeric((meanQscores_dustyNA$Dusty < cutoff)+3)
    meanQscores_dustyNA$cluster <- as.factor((clusterNA))
    clusterCentersNA <- as.data.frame(rbind(c(Q=mean(meanQscores_dustyNA$meanQ[meanQscores_dustyNA$cluster==3]), D=mean(meanQscores_dustyNA$Dusty[meanQscores_dustyNA$cluster==3])),c(Q=mean(meanQscores_dustyNA$meanQ[meanQscores_dustyNA$cluster==4]), D=mean(meanQscores_dustyNA$Dusty[meanQscores_dustyNA$cluster==4]))))
  }

  #Overlayed plot
  if (verbose==TRUE){
    cat('Building plot of all clusters\n')
  }
  
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
                labs(color="Cluster",shape=paste("Dataset\n(E-score = ",format(escore,scientific=T),")",sep="")) +
                ggtitle(paste('Grouping of BLAST Results of a ', amplicon_length,' (bp) \n Amplicon by Complexity', sep='')) +
                guides(size=FALSE) +
                theme_minimal())
  if (data==FALSE) {
    if (verbose==TRUE){
      cat('Returning ggplot object\n')
    }
    return(histcomb)
  }
  if (data==TRUE) {
    histcomb #print ggplot
    #rearrange new data frame to output all calculated values
    obj_out = setNames(data.frame(matrix(ncol=3,nrow=length(dna))),c("meanQ","Dusty","cluster"))
    obj_out[,] = as.numeric(0)

    obj_out$meanQ[idxNA] <- meanQscores_dustyNA$meanQ
    obj_out$Dusty[idxNA] <- meanQscores_dustyNA$Dusty
    obj_out$cluster[idxNA] <- as.numeric_version(meanQscores_dustyNA$cluster)
    
    obj_out$meanQ[idxHits] <- meanQscores_dustyHits$meanQ
    obj_out$Dusty[idxHits] <- meanQscores_dustyHits$Dusty
    obj_out$cluster[idxHits] <- as.numeric_version(meanQscores_dustyHits$cluster)

    if (verbose==TRUE){
      cat('Returning data frame\n')
    }
    return(obj_out)
  }
}


