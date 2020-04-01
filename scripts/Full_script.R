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

# Download SRA File:
#srr=c('SRR11043481')
#system(paste('fastq-dump', srr, sep=' '))


# Read taxonomy database
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")


# read fastq
#dna = readFastq(paste(srr, '.fastq',sep=''))

#read O'Hara SRRs:
srr1=c('SRX7820009') #OH
srr2=c('SRX7820008') #OH
srr3=c('SRX7820007') #OH
#load fastqs
system(paste('fastq-dump', srr1, sep=' ')) 
system(paste('fastq-dump', srr2, sep=' '))
system(paste('fastq-dump', srr3, sep=' '))
#read fastqs
dna1 = readFastq('./data/', pattern=srr1)
dna2 = readFastq('./data/', pattern=srr2)
dna3 = readFastq('./data/', pattern=srr3)

Full123 = c(dna1,dna2,dna3)
dnaMerge <- do.call('append', Full123)
remove(Full123)

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
ggsave(widthplot, file='outputs/readlengths.png')

# plot qscores
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE) #apply the function mean() across 
avgscores = as.data.frame(avgscores)
(qscores = ggplot(avgscores) +
    geom_histogram(aes(x=avgscores), binwidth=0.2) +
    theme_linedraw() +
    xlab('Quality Score') +
    ggtitle('Per Read Average Quality'))
ggsave(qscores, file='outputs/quality.png')


## blast
bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa")
cl <- predict(bl, reads, BLAST_args = '-num_threads 12 -evalue 1e-100')
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

cluster <- new_cluster(11)
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

ggsave(gr, file='outputs/grid_taxonomy.png', height = 9, width = 7, dpi=500)

## Section 2
print('Splitting into reads with and without genus hits')

HitsQID_group = cltax %>%
  filter(!is.na(genus)) %>% # get rid of NA genera
  group_by(QueryID) %>%#group
  slice(1)
  #summarize(count=n())

NAQID_group = cltax %>%
  filter(is.na(genus)) %>%
  group_by(QueryID) %>%
  slice(1)

#NAHitsDiff = setdiff(HitsQID_group,NAQID_group)
NAHitsDiff = setdiff(NAQID_group$QueryID,HitsQID_group$QueryID) #QueryIDs in NA_group that have no other matches



len = lengths(HitsQID_group)
seq1 = seq(1:len[1])

#make artificial list with length of HitsQID_group_S
regexp <- "[[:digit:]]+"
seq2 = str_extract(HitsQID_group[[1,]], regexp)

missingint = setdiff(seq1,seq2) #QueryIDs that have no matches
missing = paste("Query_", sep="", missingint)

############ Quality scores of Hits

idxNA= missingint # indexes of reads that have no genus matches (script_main)
idxHits = sort(as.integer(seq2)) # indexes of reads that have at least one genus match

##NA
QscoresNA = quality(dnaMerge[idxNA]) # parse quality scores
numQscoresNA = as(QscoresNA, "matrix") # converts to numeric scores automatically
meanQscoresNA_temp = rowMeans(numQscoresNA,na.rm=TRUE)
meanQscoresNA = as.data.frame(meanQscoresNA_temp)
remove(meanQscoresNA_temp)

##Hits
QscoresHits = quality(dnaMerge[idxHits])
numQscoresHits = as(QscoresHits, "matrix")
meanQscoresHits_temp = rowMeans(numQscoresHits,na.rm=TRUE)
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
  geom_histogram(data = meanQscoresHits,aes(x=meanQscoresHits$meanQscoresHits_temp, color="Genus Hits"),
                 bins = 100, 
                 fill="pink",
                 position="dodge"
  ) +
  geom_histogram(data = meanQscoresNA, aes(x=meanQscoresNA$meanQscoresNA_temp, color="No Genus Hits"),
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

histcomb3 = (ggplot(meanQscores_dustyNA) + 
  stat_density_2d(data=meanQscores_dustyNA,
                  aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
  geom_point(data=clusterCentersNA, 
             aes(x=clusterCentersNA$meanQ,y=clusterCentersNA$Dusty,shape=as.factor("No Genus Hits"),size=3)) +
  #
  stat_density_2d(data=meanQscores_dustyHits,
                  aes(x=meanQ,y=Dusty,color=as.factor(cluster))) +
  geom_point(data=clusterCentersHits,
             aes(x=clusterCentersHits$meanQ,y=clusterCentersHits$Dusty,shape=as.factor("Genus Hits"),size=3)) +
  scale_y_log10() +
  labs(color="Cluster",shape="Dataset") +
  guides(size=FALSE) +
  theme_minimal())


gr2 = plot_grid(histcomb1, histcomb2, ncol=1, nrow=2)

ggsave(gr2, file='outputs/grid_q_dusty.png', height = 9, width = 7, dpi=500)

ggsave(histcomb3, file='outputs/q_dusty_clusters.png', height = 9, width = 7, dpi=500)


save.image(file='outputs/save.RData',compress="gzip")

