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

#devtools::install_github("tidyverse/multidplyr")
#devtools::install_github("mhahsler/rBLAST")

#USER VARIABLES
srr=c('SRR11206993') 
cores = c(12)
escore = c('1e-50')

print(paste('Analyzing Accession code ',srr, sep=''))
print(paste('Multicore processing on ',cores,' cores', sep=''))
print('NCBI database')
print(paste('BLAST threshold (E-score = ',escore,')' ,sep=''))


# Read taxonomy database
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")

#SRR repository: PRJNA605442

#load fastqs
system(paste('fastq-dump', srr, sep=' ')) 
#read fastqs
dna = readFastq('.', pattern=srr)

# todo: move files to data

#Full123 = c(dna1,dna2,dna3)
#dnaMerge <- do.call('append', Full123)
#remove(Full123)
dnaMerge = dna

readMerge = sread(dnaMerge, id=id(dnaMerge))
reads = sread(dnaMerge)
qscores = quality(dnaMerge) 

# plot readlength
widths = as.data.frame(reads@ranges@width)
(widthplot <- ggplot(widths) +
    geom_histogram(aes(x=reads@ranges@width), binwidth = 10) + 
    theme_linedraw() + 
    xlab('Read Length (bp)') +
    xlim(0,2000) +
    ggtitle('Read length distribution for 550bp amplicon'))
out1 <- paste("outputs/Swamp/readlengths_",srr,".png", sep="")
ggsave(widthplot, file=out1)

# plot qscores
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE) #apply the function mean() across 
avgscores = as.data.frame(avgscores)
(qscores = ggplot(avgscores) +
    geom_histogram(aes(x=avgscores), binwidth=0.2) +
    theme_linedraw() +
    xlab('Quality Score') +
    ggtitle('Per Read Average Quality'))
out2 <- paste("outputs/Swamp/quality_",srr,".png", sep="")
ggsave(qscores, file=out2)


## PLASTID DB
## 1. make database
#makeblastdb("./data/plastid.3.1.genomic.fna", dbtype = "nucl")

## 3. open database
#plastid_genomes_db <- blast("./data/plastid.3.1.genomic.fna")

## 4. perform search (first sequence in the db should be a perfect match)
## blast
bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
print('Begin BLAST...')
b_args <- paste('-num_threads ',cores,' -evalue ',escore, sep="")
cl <- predict(bl, reads, BLAST_args = b_args)
accid = as.character(cl$SubjectID) # accession IDs of BLAST hits

# Plot results

#takes accession number and gets the taxonomic ID
ids<-accessionToTaxa(accid, '/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
#taxlist displays the taxonomic names from each ID #
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)

cltax=cbind(cl,taxlist)


lca2 = function(x) {
  require(dplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  x = x %>% filter(!is.na(superkingdom)) %>% filter(superkingdom != 'Bacteria') # need to deal with this more generally
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  numcount = countshnames==1
  lastuni = tail(names(shortnames[numcount==T]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  newtax <- as.list(ifelse(countshnames==1,shortnames,NA))
  
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni) %>%
    mutate(superkingdom = newtax$superkingdom) %>%
    mutate(phylum = newtax$phylum) %>%
    mutate(class = newtax$class) %>%
    mutate(order = newtax$order) %>%
    mutate(family = newtax$family) %>%
    mutate(genus = newtax$genus) %>%
    mutate(species = newtax$species)
  return(ret)
}

blasthits = cltax

cluster <- new_cluster(cores)
cluster_library(cluster, 'dplyr')
cluster_copy(cluster, 'lca2')

#split groups across multiple CPU cores
prepare = blasthits %>%
  group_by(QueryID) %>%
  partition(cluster) %>%  #split groups into cluster units
  do({lca2(.)}) %>%
  collect()

#count reads matching species:
spcount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(species) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(species))


#count reads matching to genera:
gencount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(genus) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(genus))

#count reads matching to families
famcount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(family) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(family))


(p1 = ggplot(spcount[1:50,]) +
    geom_col(aes(x=fct_reorder(species, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Species')
)

(p2 = ggplot(gencount[1:50,]) +
    geom_col(aes(x=fct_reorder(genus, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Genus')
)

(p3 = ggplot(famcount[1:50,]) +
    geom_col(aes(x=fct_reorder(family, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Family')
)

gr = plot_grid(p1, p2, p3, ncol=1, nrow=3)

out3 <- paste("outputs/Swamp/grid_taxonomy_",srr,".png", sep="")
ggsave(gr, file=out3, height = 9, width = 7, dpi=500)

out6 <- paste("outputs/Swamp/save_",srr,".RData", sep="")
save.image(file=out6, compress="gzip")

## Section 2
print('Splitting into reads with and without BLAST hits')

HitsQID_group = cltax %>%
  #filter(!is.na(genus)) %>% # get rid of NA genera
  group_by(QueryID) 


#extract read indices
regexp <- "[[:digit:]]+"
idxHits = as.numeric(unique(str_extract(HitsQID_group$QueryID, regexp)))
#idxNA = as.numeric(str_extract(NAHitsDiff, regexp))

idxNA = unique(setdiff((seq(1:length(dnaMerge))),idxHits))

NA_reads = dnaMerge[idxNA]

############ Quality scores of Hits

#idxNA= missingint # indexes of reads that have no genus matches (script_main)
#idxHits = sort(as.integer(seq2)) # indexes of reads that have at least one genus match

##NA
QscoresNA = quality(dnaMerge[idxNA]) # parse quality scores
numQscoresNA = as(QscoresNA, "matrix") # converts to numeric scores automatically
meanQscoresNA_temp = rowMeans(numQscoresNA,na.rm=TRUE)
meanQscoresNA = as.data.frame(meanQscoresNA_temp)
remove(meanQscoresNA_temp)

##Hits
QscoresHits = quality(dnaMerge[idxHits])
numQscoresHits = as(QscoresHits, "matrix")
#meanQscoresHits_temp = rowMeans(numQscoresHits,na.rm=TRUE)
meanQscoresHits_temp = apply(numQscoresHits, 1, mean, na.rm=TRUE) #apply the function mean() across 
meanQscoresHits = as.data.frame(meanQscoresHits_temp,col.names=c('QmeansHits'))
remove(meanQscoresHits_temp)

#meanQscoresAll = as.data.frame(meanQscoresNA,meanQscoresHits)


# ggplot(meanQscoresNA, aes(x=meanQscoresNA$meanQscoresNA_temp)) +
#   geom_histogram(bins = 100, 
#                  color="darkblue", 
#                  fill="lightblue",
#                  position="dodge")
# 
# ggplot(meanQscoresHits, aes(x=meanQscoresHits$meanQscoresHits_temp)) +
#   geom_histogram(bins = 100, 
#                  color="darkred", 
#                  fill="pink",
#                  position="dodge")
# 
# ggplot(meanQscoresNA, aes(x=meanQscoresNA$meanQscoresNA_temp)) +
#   geom_histogram(bins = 100, 
#                  color="darkblue", 
#                  fill="lightblue",
#                  position="dodge")

print('Plotting Qscores')
histcomb1 = (ggplot(meanQscoresHits) +
  geom_histogram(data = meanQscoresHits,aes(x=meanQscoresHits$meanQscoresHits_temp, color="BLAST Hits"),
                 bins = 100, 
                 fill="pink",
                 position="dodge"
  ) +
  geom_histogram(data = meanQscoresNA, aes(x=meanQscoresNA$meanQscoresNA_temp, color="No BLAST Hits"),
                 bins = 100, 
                 fill="lightblue",
                 position="dodge", 
                 alpha=0.5
  )+
  xlab('meanQscores')) +
  labs(color='Dataset')
print('...')
complexNA = dustyScore(dnaMerge[idxNA])
complexNA = as.data.frame(complexNA)
meanQscores_dustyNA = cbind(meanQ=(meanQscoresNA$meanQscoresNA_temp),Dusty=(complexNA$complexNA))
meanQscores_dustyNA = as.data.frame(meanQscores_dustyNA)

complexHits = dustyScore(dnaMerge[idxHits])
complexHits = as.data.frame(complexHits)
meanQscores_dustyHits = cbind(meanQ=(meanQscoresHits$meanQscoresHits_temp),Dusty=(complexHits$complexHits))
meanQscores_dustyHits = as.data.frame(meanQscores_dustyHits)

print('Plotting Q and Dusty score density plots')
histcomb2 = (ggplot(meanQscores_dustyNA) + 
  stat_density_2d(data=meanQscores_dustyNA,aes(x=meanQ,y=Dusty), color= 'lightblue') +
  stat_density_2d(data=meanQscores_dustyHits,aes(x=meanQ,y=Dusty), color='pink') +
  scale_y_log10() +
  theme_minimal())

clusterHits = kmeans(meanQscores_dustyHits,2)
meanQscores_dustyHits$cluster <- as.factor(clusterHits$cluster)
clusterCentersHits <- as.data.frame(clusterHits$centers,clusterHits$cluster)

# ggplot(meanQscores_dustyHits) + 
#   stat_density_2d(data=meanQscores_dustyHits, aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
#   scale_y_log10() +
#   geom_point(data=clusterCentersHits, aes(x=clusterCentersHits$meanQ,y=clusterCentersHits$Dusty)) +
#   theme_minimal()

print('...')

clusterNA = kmeans(meanQscores_dustyNA,2)
meanQscores_dustyNA$cluster <- as.factor((clusterNA$cluster+2))
clusterCentersNA <- as.data.frame(clusterNA$centers,clusterNA$cluster)

# ggplot(meanQscores_dustyNA) + 
#   stat_density_2d(data=meanQscores_dustyNA, aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
#   scale_y_log10() +
#   geom_point(data=clusterCentersNA, aes(x=clusterCentersNA$meanQ,y=clusterCentersNA$Dusty)) +
#   theme_minimal()

#Overlayed plot
print('Plotting all clusters')

histcomb3_hits_label  <- paste("BLAST Hits (",length(idxHits),")", sep="")
histcomb3_NA_label  <- paste("No BLAST Hits (",length(idxNA),")", sep="")
histcomb3 = (ggplot(meanQscores_dustyNA) + 
  stat_density_2d(data=meanQscores_dustyHits,
                  aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
  #geom_point(data=meanQscores_dustyHits,
                  #aes(x=meanQ,y=Dusty,color=as.factor(cluster)),size=0.5) +
  geom_point(data=clusterCentersHits,
             aes(x=clusterCentersHits$meanQ,y=clusterCentersHits$Dusty,shape=as.factor(histcomb3_hits_label),size=3)) +
    #
    stat_density_2d(data=meanQscores_dustyNA,
                    aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
    #geom_point(data=meanQscores_dustyNA,
               #aes(x=meanQ,y=Dusty,color=as.factor(cluster)),size=0.5) +
    geom_point(data=clusterCentersNA, 
               aes(x=clusterCentersNA$meanQ,y=clusterCentersNA$Dusty,shape=as.factor(histcomb3_NA_label),size=3)) +
    
  scale_y_log10() +
  labs(color="Cluster",shape="Dataset") +
  guides(size=FALSE) +
  theme_minimal())


gr2 = plot_grid(histcomb1, histcomb2, ncol=1, nrow=2)

out4 <- paste("outputs/Swamp/grid_Qmean_Dusty_",srr,".png", sep="")
ggsave(gr2, file=out4, height = 9, width = 7, dpi=500)

out5 <- paste("outputs/Swamp/clusters_Qmean_Dusty_",srr,".png", sep="")
ggsave(histcomb3, file=out5, height = 9, width = 7, dpi=500)

out6 <- paste("outputs/Swamp/save_",srr,".RData", sep="")
save.image(file=out6, compress="gzip")


